/*
 * SamSegmentMapper.cpp
 *
 *  Created on: 16 Apr 2016
 *      Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical
 * Bioinformatics, Bioinformatics and Omics Data Analytics and currently
 * Neuroblastoma Genomics) Copyright (C) 2018 Umut H. Toprak, Matthias
 * Schlesner, Roland Eils and DKFZ Heidelberg
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *      LICENSE: GPL
 */

#include "SamSegmentMapper.h"
#include <cmath>
#include <iostream>
#include <limits>
#include "GlobalAppConfig.h"

namespace sophia {

    SamSegmentMapper::SamSegmentMapper(ChrSize defaultReadLengthIn)
        : STARTTIME{time(nullptr)},
          PROPER_PAIR_COMPENSATION_MODE{Breakpoint::PROPER_PAIR_COMPENSATION_MODE},
          DISCORDANT_LEFT_RANGE{static_cast<ChrSize>(round(defaultReadLengthIn * 3))},
          DISCORDANT_RIGHT_RANGE{static_cast<ChrSize>(round(defaultReadLengthIn * 2.51))},
          printedBps{0u},
          chrIndexCurrent{0},
          minPos{std::numeric_limits<ChrSize>::max()},
          maxPos{std::numeric_limits<ChrSize>::min()},
          breakpointsCurrent{},
          discordantAlignmentsPool{},
          discordantAlignmentCandidatesPool{},
          discordantLowQualAlignmentsPool{} {}

    void
    SamSegmentMapper::parseSamStream() {
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        while (true) {
            auto alignment = make_shared<Alignment>();

            if (!chrConverter.isCompressedMref(alignment->getChrIndex())) {
                continue;
            }

            if (alignment->isValidLine()) {
                if (alignment->getChrIndex() != chrIndexCurrent) {
                    switchChromosome(*alignment);
                }
                alignment->continueConstruction();
                printBps(alignment->getStartPos());
                incrementCoverages(*alignment);
                assignBps(alignment);
            } else {
                break;
            }
        }
        // EOF event for the samtools pipe. printing the end of the very last
        // chromosome,
        printBps(numeric_limits<int>::max());
    }

    void
    SamSegmentMapper::switchChromosome(const Alignment &alignment) {
        // As we entered a new chromosome here, now print the previous chromosome's
        // unprinted regions
        if (chrIndexCurrent != 0) {
            printBps(numeric_limits<int>::max());
        }
        chrIndexCurrent = alignment.getChrIndex();
        breakpointsCurrent.clear();
        coverageProfiles.clear();
        discordantAlignmentsPool.clear();
        if (PROPER_PAIR_COMPENSATION_MODE) {
            discordantAlignmentCandidatesPool.clear();
        }
        discordantLowQualAlignmentsPool.clear();
        minPos = std::numeric_limits<ChrSize>::max();
        maxPos = std::numeric_limits<ChrSize>::min();
    }

    /** Does a lot of stuff, and -- just by the way -- also prints the results to stdout :(.
      **/
    void
    SamSegmentMapper::printBps(ChrSize alignmentStart) {
        for (auto &bp : breakpointsCurrent) {
            if (!bp.second.isCovFinalized() && ((bp.first) + 1 < alignmentStart)) {
                unsigned long posDiff = static_cast<unsigned long>(bp.first - minPos);
                if (bp.first != minPos) {
                    bp.second.setLeftCoverage(
                        coverageProfiles[static_cast<unsigned long>(bp.first - minPos - 1)].getCoverage());
                } else {
                    bp.second.setLeftCoverage(0);
                }
                bp.second.setRightCoverage(
                    coverageProfiles[posDiff].getCoverage());
                bp.second.setNormalSpans(
                    coverageProfiles[posDiff].getNormalSpans());
                bp.second.setLowQualSpansSoft(
                    coverageProfiles[posDiff].getLowQualSpansSoft());
                bp.second.setLowQualSpansHard(
                    coverageProfiles[posDiff].getLowQualSpansHard());
                bp.second.setUnpairedBreaksSoft(
                    coverageProfiles[posDiff].getNormalBpsSoft());
                bp.second.setUnpairedBreaksHard(
                    coverageProfiles[posDiff].getNormalBpsHard());
                bp.second.setBreaksShortIndel(
                    coverageProfiles[posDiff].getNormalBpsShortIndel());
                bp.second.setLowQualBreaksSoft(
                    coverageProfiles[posDiff].getLowQualBpsSoft());
                bp.second.setLowQualBreaksHard(
                    coverageProfiles[posDiff].getLowQualBpsHard());
                bp.second.setCovFinalized(true);
            }
        }
        if (minPos != std::numeric_limits<ChrSize>::max()) {
            while (minPos + 2 + DISCORDANT_LEFT_RANGE < alignmentStart) {
                if (minPos != maxPos) {
                    coverageProfiles.pop_front();
                    ++minPos;
                } else {
                    coverageProfiles.clear();
                    minPos = std::numeric_limits<ChrSize>::max();
                    maxPos = std::numeric_limits<ChrSize>::min();
                    break;
                }
            }
        }
        for (auto bpIt = breakpointsCurrent.begin();
             bpIt != breakpointsCurrent.end();) {
            if ((bpIt->first) + DISCORDANT_RIGHT_RANGE < alignmentStart) {
                if (bpIt->second.finalizeBreakpoint(   // Side effect: prints the breakpoint!
                        discordantAlignmentsPool, discordantLowQualAlignmentsPool,
                        discordantAlignmentCandidatesPool)) {
                    ++printedBps;
                }
                bpIt = breakpointsCurrent.erase(bpIt);
            } else {
                break;
            }
        }
        while (!discordantAlignmentsPool.empty() &&
               (discordantAlignmentsPool.front().readStartPos +
                    DISCORDANT_LEFT_RANGE + DISCORDANT_RIGHT_RANGE <
                alignmentStart)) {
            discordantAlignmentsPool.pop_front();
        }
        if (PROPER_PAIR_COMPENSATION_MODE) {
            while (!discordantAlignmentCandidatesPool.empty() &&
                   (discordantAlignmentCandidatesPool.front().readStartPos +
                        DISCORDANT_LEFT_RANGE + DISCORDANT_RIGHT_RANGE <
                    alignmentStart)) {
                discordantAlignmentCandidatesPool.pop_front();
            }
        }
        while (!discordantLowQualAlignmentsPool.empty() &&
               (discordantLowQualAlignmentsPool.front().readStartPos +
                    DISCORDANT_LEFT_RANGE + DISCORDANT_RIGHT_RANGE <
                alignmentStart)) {
            discordantLowQualAlignmentsPool.pop_front();
        }
    }

    void
    SamSegmentMapper::incrementCoverages(const Alignment &alignment) {
        if (minPos == std::numeric_limits<ChrSize>::max()) {
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                coverageProfiles.emplace_back();
            }
            minPos = alignment.getStartPos();
            maxPos = alignment.getEndPos() - 1;
        } else {
            while (alignment.getStartPos() > maxPos) {
                coverageProfiles.emplace_back();
                ++maxPos;
            }
        }
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        switch (alignment.getReadType()) {
        case 0:
        case 3:
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                if (i > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementCoverage();
                coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementNormalSpans();
            }
            if (PROPER_PAIR_COMPENSATION_MODE) {
                discordantAlignmentCandidatesPool.emplace_back(
                    alignment.getStartPos(),
                    alignment.getEndPos(),
                    std::numeric_limits<ChrIndex>::max(), // mateChrIndexIn
                    std::numeric_limits<ChrSize>::max(), // mateStartPosIn
                    -1, // sourceType
                    false);
            }
            break;
        case 1:
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                if (i > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementCoverage();
                coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementNormalSpans();
            }
            break;
        case 4:
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                if (i > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementCoverage();
                coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementNormalSpans();
            }
            if (!chrConverter.isTechnical(alignment.getMateChrIndex())
                && !chrConverter.isInBlockedRegion(alignment.getMateChrIndex(),
                                                   alignment.getMatePos())) {

                if (PROPER_PAIR_COMPENSATION_MODE) {
                    discordantAlignmentCandidatesPool.emplace_back(
                        alignment.getStartPos(),
                        alignment.getEndPos(),
                        alignment.getMateChrIndex(),
                        alignment.getMatePos(),
                        2,   // TODO is this a chromosome index
                        alignment.isInvertedMate());
                }
                if (!alignment.isNullMapq()) {
                    discordantAlignmentsPool.emplace_back(
                        alignment.getStartPos(),
                        alignment.getEndPos(),
                        alignment.getMateChrIndex(),
                        alignment.getMatePos(),
                        2,   // TODO is this a chromosome index
                        alignment.isInvertedMate());
                } else {
                    discordantLowQualAlignmentsPool.emplace_back(
                        alignment.getStartPos(), alignment.getEndPos(),
                        alignment.getMateChrIndex(), alignment.getMatePos(), 2,
                        alignment.isInvertedMate());
                }
            }
            break;
        case 2:
            if (!(alignment.isLowMapq() || alignment.isNullMapq())) {
                for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                    if (i > maxPos) {
                        coverageProfiles.emplace_back();
                        ++maxPos;
                    }
                }
            } else {
                for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                    if (i > maxPos) {
                        coverageProfiles.emplace_back();
                        ++maxPos;
                    }
                    coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementLowQualSpansHard();
                }
            }
            break;
        case 5:
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                if (i > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                coverageProfiles[static_cast<unsigned long>(i - minPos)].incrementLowQualSpansSoft();
            }
            if (!alignment.isSupplementary() &&
                !chrConverter.isTechnical(alignment.getMateChrIndex()) &&
                alignment.isDistantMate()) {
                if (!chrConverter.isInBlockedRegion(alignment.getMateChrIndex(),
                                                    alignment.getMatePos())) {
                    discordantLowQualAlignmentsPool.emplace_back(
                        alignment.getStartPos(),
                        alignment.getEndPos(),
                        alignment.getMateChrIndex(),
                        alignment.getMatePos(),
                        2,   // TODO Is this a chromosome index?
                        alignment.isInvertedMate(),
                        alignment.getReadBreakpoints());
                }
            }
            break;
        default:
            break;
        }

        switch (alignment.getReadType()) {
        case 1:
            for (auto j = 0u; j < alignment.getReadBreakpoints().size(); ++j) {
                ChrSize bpPos = alignment.getReadBreakpoints()[j];
                if (bpPos > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                switch (alignment.getReadBreakpointTypes()[j]) {
                case 'S':
                    if (bpPos == alignment.getStartPos()) {
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].decrementNormalSpans();
                    }
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementNormalBpsSoft();
                    break;
                case 'I':
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementNormalBpsShortIndel();
                    break;
                case 'D':
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementNormalBpsShortIndel();
                    for (signed int k = 0; k != alignment.getReadBreakpointsSizes()[j]; ++k) {
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos + k)].decrementNormalSpans();
                    }
                    break;
                default:
                    break;
                }
            }
            break;
        case 3:
            for (auto j = 0u; j < alignment.getReadBreakpoints().size(); ++j) {
                auto bpPos = alignment.getReadBreakpoints()[j];
                if (bpPos > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                switch (alignment.getReadBreakpointTypes()[j]) {
                case 'I':
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementNormalBpsShortIndel();
                    break;
                case 'D':
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementNormalBpsShortIndel();
                    for (signed int k = 0; k != alignment.getReadBreakpointsSizes()[j]; ++k) {
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos + k)].decrementNormalSpans();
                    }
                    break;
                default:
                    break;
                }
            }
            break;
        case 2:
            if (!(alignment.isLowMapq() || alignment.isNullMapq())) {
                for (auto j = 0u; j < alignment.getReadBreakpoints().size(); ++j) {
                    auto bpPos = alignment.getReadBreakpoints()[j];
                    if (bpPos > maxPos) {
                        coverageProfiles.emplace_back();
                        ++maxPos;
                    }
                    if (alignment.getReadBreakpointTypes()[j] == 'H') {
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementNormalBpsHard();
                    }
                }
            } else {
                for (auto j = 0u; j < alignment.getReadBreakpoints().size(); ++j) {
                    auto bpPos = alignment.getReadBreakpoints()[j];
                    if (bpPos > maxPos) {
                        coverageProfiles.emplace_back();
                        ++maxPos;
                    }
                    switch (alignment.getReadBreakpointTypes()[j]) {
                    case 'S':
                        if (bpPos == alignment.getStartPos()) {
                            coverageProfiles[static_cast<unsigned long>(bpPos - minPos)]
                                .decrementLowQualSpansHard();
                        }
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementLowQualBpsHard();
                        break;
                    case 'H':
                        if (bpPos == alignment.getStartPos()) {
                            coverageProfiles[static_cast<unsigned long>(bpPos - minPos)]
                                .decrementLowQualSpansHard();
                        }
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementLowQualBpsHard();
                        break;
                    case 'I':
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementLowQualBpsHard();
                        break;
                    case 'D':
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementLowQualBpsHard();
                        for (signed int k = 0; k != alignment.getReadBreakpointsSizes()[j]; ++k) {
                            coverageProfiles[static_cast<unsigned long>(bpPos - minPos + k)]
                                .decrementLowQualSpansHard();
                        }
                        break;
                    default:
                        break;
                    }
                }
            }
            break;
        case 5:
            for (auto j = 0u; j < alignment.getReadBreakpoints().size(); ++j) {
                auto bpPos = alignment.getReadBreakpoints()[j];
                if (bpPos > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                switch (alignment.getReadBreakpointTypes()[j]) {
                case 'S':
                case 'H':
                    if (bpPos == alignment.getStartPos()) {
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos)]
                            .decrementLowQualSpansSoft();
                    }
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementLowQualBpsSoft();
                    break;
                case 'I':
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementLowQualBpsSoft();
                    break;
                case 'D':
                    coverageProfiles[static_cast<unsigned long>(bpPos - minPos)].incrementLowQualBpsSoft();
                    for (signed int k = 0; k != alignment.getReadBreakpointsSizes()[j]; ++k) {
                        coverageProfiles[static_cast<unsigned long>(bpPos - minPos + k)]
                            .decrementLowQualSpansSoft();
                    }
                    break;
                default:
                    break;
                }
            }
            break;
        default:
            break;
        }
    }

    void
    SamSegmentMapper::assignBps(shared_ptr<Alignment> &alignment) {
        switch (alignment->getReadType()) {
        case 1:
            for (auto i = 0u; i < alignment->getReadBreakpoints().size(); ++i) {
                if (alignment->getReadBreakpointTypes()[i] == 'S') {
                    ChrSize bpLoc = alignment->getReadBreakpoints()[i];
                    auto it = breakpointsCurrent.find(bpLoc);
                    if (it == breakpointsCurrent.end()) {
                        auto newIt = breakpointsCurrent.emplace(
                            piecewise_construct, forward_as_tuple(bpLoc),
                            forward_as_tuple(chrIndexCurrent, bpLoc));
                        newIt.first->second.addSoftAlignment(alignment);
                    } else {
                        it->second.addSoftAlignment(alignment);
                    }
                }
            }
            break;
        case 2:
            for (auto i = 0u; i < alignment->getReadBreakpoints().size(); ++i) {
                if (alignment->getReadBreakpointTypes()[i] == 'H') {
                    auto bpLoc = alignment->getReadBreakpoints()[i];
                    auto it = breakpointsCurrent.find(bpLoc);
                    if (it == breakpointsCurrent.end()) {
                        auto newIt = breakpointsCurrent.emplace(
                            piecewise_construct, forward_as_tuple(bpLoc),
                            forward_as_tuple(chrIndexCurrent, bpLoc));
                        newIt.first->second.addHardAlignment(alignment);
                    } else {
                        it->second.addHardAlignment(alignment);
                    }
                }
            }
            break;
        default:
            break;
        }
    }

    // void SamSegmentMapper::printMetadata(int ISIZESIGMALEVEL) {
    //	auto elapsedTime = div(difftime(time(nullptr), STARTTIME), 60);
    //	*metaOutputHandle << "#Using soft/hard clip length threshold " <<
    //Alignment::CLIPPED_NUCLEOTIDE_COUNT_THRESHOLD << endl; 	*metaOutputHandle <<
    //"#Using low quality clipped overhang length threshold " <<
    //Alignment::LOW_QUAL_CLIP_THRESHOLD << endl; 	*metaOutputHandle << "#Using Base
    //Quality Threshold " << Alignment::BASE_QUALITY_THRESHOLD << endl;
    //	*metaOutputHandle << "#Using Base Quality Threshold Low " <<
    //Alignment::BASE_QUALITY_THRESHOLD_LOW << endl; 	*metaOutputHandle << "#Using
    //sigmas = " << ISIZESIGMALEVEL << " away from the median insert size for
    //'distant' classification" << endl; 	*metaOutputHandle << "#Using minimum isize
    //for 'distant' classification = " << Alignment::ISIZEMAX << " bps" << endl;
    //	*metaOutputHandle << "#Using minimum reads supporting a breakpoint " <<
    //Breakpoint::BP_SUPPORT_THRESHOLD << endl; 	*metaOutputHandle << "#Using minimum
    //reads supporting a discordant mate contig " << Breakpoint::BP_SUPPORT_THRESHOLD
    //<< endl; 	*metaOutputHandle << "#Using (-F 0x600 -f 0x001)" << endl;
    //	*metaOutputHandle << "#done\t" << printedBps << " lines printed in " <<
    //elapsedTime.quot << " minutes, " << elapsedTime.rem << " seconds" << endl;
    // }

} /* namespace sophia */
