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

    SamSegmentMapper::SamSegmentMapper(int defaultReadLengthIn)
        : STARTTIME{time(nullptr)},
          PROPERPARIRCOMPENSATIONMODE{Breakpoint::PROPERPAIRCOMPENSATIONMODE},
          DISCORDANTLEFTRANGE{static_cast<int>(round(defaultReadLengthIn * 3))},
          DISCORDANTRIGHTRANGE{static_cast<int>(round(defaultReadLengthIn * 2.51))},
          printedBps{0u},
          chrIndexCurrent{0},
          minPos{-1},
          maxPos{-1},
          breakpointsCurrent{},
          discordantAlignmentsPool{},
          discordantAlignmentCandidatesPool{},
          discordantLowQualAlignmentsPool{} {}

    void
    SamSegmentMapper::parseSamStream() {
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        while (true) {
            auto alignment = make_shared<Alignment>();

            if (!chrConverter.isCompressedMrefIndex(alignment->getChrIndex())) {
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
        if (PROPERPARIRCOMPENSATIONMODE) {
            discordantAlignmentCandidatesPool.clear();
        }
        discordantLowQualAlignmentsPool.clear();
        minPos = -1;
        maxPos = -1;
    }

    void
    SamSegmentMapper::printBps(int alignmentStart) {
        for (auto &bp : breakpointsCurrent) {
            if (!bp.second.isCovFinalized() && ((bp.first) + 1 < alignmentStart)) {
                auto posDiff = bp.first - minPos;
                if (bp.first != minPos) {
                    bp.second.setLeftCoverage(
                        coverageProfiles[bp.first - minPos - 1].getCoverage());
                } else {
                    bp.second.setLeftCoverage(0);
                }
                bp.second.setRightCoverage(coverageProfiles[posDiff].getCoverage());
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
        if (minPos != -1) {
            while (minPos + 2 + DISCORDANTLEFTRANGE < alignmentStart) {
                if (minPos != maxPos) {
                    coverageProfiles.pop_front();
                    ++minPos;
                } else {
                    coverageProfiles.clear();
                    minPos = -1;
                    maxPos = -1;
                    break;
                }
            }
        }
        for (auto bpIt = breakpointsCurrent.begin();
             bpIt != breakpointsCurrent.end();) {
            if ((bpIt->first) + DISCORDANTRIGHTRANGE < alignmentStart) {
                if (bpIt->second.finalizeBreakpoint(
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
                    DISCORDANTLEFTRANGE + DISCORDANTRIGHTRANGE <
                alignmentStart)) {
            discordantAlignmentsPool.pop_front();
        }
        if (PROPERPARIRCOMPENSATIONMODE) {
            while (!discordantAlignmentCandidatesPool.empty() &&
                   (discordantAlignmentCandidatesPool.front().readStartPos +
                        DISCORDANTLEFTRANGE + DISCORDANTRIGHTRANGE <
                    alignmentStart)) {
                discordantAlignmentCandidatesPool.pop_front();
            }
        }
        while (!discordantLowQualAlignmentsPool.empty() &&
               (discordantLowQualAlignmentsPool.front().readStartPos +
                    DISCORDANTLEFTRANGE + DISCORDANTRIGHTRANGE <
                alignmentStart)) {
            discordantLowQualAlignmentsPool.pop_front();
        }
    }

    void
    SamSegmentMapper::incrementCoverages(const Alignment &alignment) {
        if (minPos == -1) {
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
                coverageProfiles[i - minPos].incrementCoverage();
                coverageProfiles[i - minPos].incrementNormalSpans();
            }
            if (PROPERPARIRCOMPENSATIONMODE) {
                discordantAlignmentCandidatesPool.emplace_back(
                    alignment.getStartPos(), alignment.getEndPos(), -1, -1, -1,
                    false);
            }
            break;
        case 1:
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                if (i > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                coverageProfiles[i - minPos].incrementCoverage();
                coverageProfiles[i - minPos].incrementNormalSpans();
            }
            break;
        case 4:
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                if (i > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                coverageProfiles[i - minPos].incrementCoverage();
                coverageProfiles[i - minPos].incrementNormalSpans();
            }
            if (!chrConverter.isIgnoredChromosome(alignment.getMateChrIndex()) &&
                !(alignment.getMateChrIndex() == 2 &&
                  (alignment.getMatePos() / 10000 == 3314))) {
                if (PROPERPARIRCOMPENSATIONMODE) {
                    discordantAlignmentCandidatesPool.emplace_back(
                        alignment.getStartPos(), alignment.getEndPos(),
                        alignment.getMateChrIndex(), alignment.getMatePos(), 2,
                        alignment.isInvertedMate());
                }
                if (!alignment.isNullMapq()) {
                    discordantAlignmentsPool.emplace_back(
                        alignment.getStartPos(), alignment.getEndPos(),
                        alignment.getMateChrIndex(), alignment.getMatePos(), 2,
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
                for (auto i = alignment.getStartPos(); i < alignment.getEndPos();
                     ++i) {
                    if (i > maxPos) {
                        coverageProfiles.emplace_back();
                        ++maxPos;
                    }
                }
            } else {
                for (auto i = alignment.getStartPos(); i < alignment.getEndPos();
                     ++i) {
                    if (i > maxPos) {
                        coverageProfiles.emplace_back();
                        ++maxPos;
                    }
                    coverageProfiles[i - minPos].incrementLowQualSpansHard();
                }
            }
            break;
        case 5:
            for (auto i = alignment.getStartPos(); i < alignment.getEndPos(); ++i) {
                if (i > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                coverageProfiles[i - minPos].incrementLowQualSpansSoft();
            }
            if (!alignment.isSupplementary() &&
                !chrConverter.isIgnoredChromosome(alignment.getMateChrIndex()) &&
                alignment.isDistantMate()) {
                if (!(alignment.getMateChrIndex() == 2 &&
                      (alignment.getMatePos() / 10000 == 3314))) {
                    discordantLowQualAlignmentsPool.emplace_back(
                        alignment.getStartPos(), alignment.getEndPos(),
                        alignment.getMateChrIndex(), alignment.getMatePos(), 2,
                        alignment.isInvertedMate(), alignment.getReadBreakpoints());
                }
            }
            break;
        default:
            break;
        }
        switch (alignment.getReadType()) {
        case 1:
            for (auto j = 0u; j < alignment.getReadBreakpoints().size(); ++j) {
                auto bpPos = alignment.getReadBreakpoints()[j];
                if (bpPos > maxPos) {
                    coverageProfiles.emplace_back();
                    ++maxPos;
                }
                switch (alignment.getReadBreakpointTypes()[j]) {
                case 'S':
                    if (bpPos == alignment.getStartPos()) {
                        coverageProfiles[bpPos - minPos].decrementNormalSpans();
                    }
                    coverageProfiles[bpPos - minPos].incrementNormalBpsSoft();
                    break;
                case 'I':
                    coverageProfiles[bpPos - minPos].incrementNormalBpsShortIndel();
                    break;
                case 'D':
                    coverageProfiles[bpPos - minPos].incrementNormalBpsShortIndel();
                    for (auto k = 0; k != alignment.getReadBreakpointsSizes()[j];
                         ++k) {
                        coverageProfiles[bpPos - minPos + k].decrementNormalSpans();
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
                    coverageProfiles[bpPos - minPos].incrementNormalBpsShortIndel();
                    break;
                case 'D':
                    coverageProfiles[bpPos - minPos].incrementNormalBpsShortIndel();
                    for (auto k = 0; k != alignment.getReadBreakpointsSizes()[j];
                         ++k) {
                        coverageProfiles[bpPos - minPos + k].decrementNormalSpans();
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
                        coverageProfiles[bpPos - minPos].incrementNormalBpsHard();
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
                            coverageProfiles[bpPos - minPos]
                                .decrementLowQualSpansHard();
                        }
                        coverageProfiles[bpPos - minPos].incrementLowQualBpsHard();
                        break;
                    case 'H':
                        if (bpPos == alignment.getStartPos()) {
                            coverageProfiles[bpPos - minPos]
                                .decrementLowQualSpansHard();
                        }
                        coverageProfiles[bpPos - minPos].incrementLowQualBpsHard();
                        break;
                    case 'I':
                        coverageProfiles[bpPos - minPos].incrementLowQualBpsHard();
                        break;
                    case 'D':
                        coverageProfiles[bpPos - minPos].incrementLowQualBpsHard();
                        for (auto k = 0;
                             k != alignment.getReadBreakpointsSizes()[j]; ++k) {
                            coverageProfiles[bpPos - minPos + k]
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
                        coverageProfiles[bpPos - minPos]
                            .decrementLowQualSpansSoft();
                    }
                    coverageProfiles[bpPos - minPos].incrementLowQualBpsSoft();
                    break;
                case 'I':
                    coverageProfiles[bpPos - minPos].incrementLowQualBpsSoft();
                    break;
                case 'D':
                    coverageProfiles[bpPos - minPos].incrementLowQualBpsSoft();
                    for (auto k = 0; k != alignment.getReadBreakpointsSizes()[j];
                         ++k) {
                        coverageProfiles[bpPos - minPos + k]
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
                    auto bpLoc = alignment->getReadBreakpoints()[i];
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
    //Alignment::CLIPPEDNUCLEOTIDECOUNTTHRESHOLD << endl; 	*metaOutputHandle <<
    //"#Using low quality clipped overhang length threshold " <<
    //Alignment::LOWQUALCLIPTHRESHOLD << endl; 	*metaOutputHandle << "#Using Base
    //Quality Threshold " << Alignment::BASEQUALITYTHRESHOLD << endl;
    //	*metaOutputHandle << "#Using Base Quality Threshold Low " <<
    //Alignment::BASEQUALITYTHRESHOLDLOW << endl; 	*metaOutputHandle << "#Using
    //sigmas = " << ISIZESIGMALEVEL << " away from the median insert size for
    //'distant' classification" << endl; 	*metaOutputHandle << "#Using minimum isize
    //for 'distant' classification = " << Alignment::ISIZEMAX << " bps" << endl;
    //	*metaOutputHandle << "#Using minimum reads supporting a breakpoint " <<
    //Breakpoint::BPSUPPORTTHRESHOLD << endl; 	*metaOutputHandle << "#Using minimum
    //reads supporting a discordant mate contig " << Breakpoint::BPSUPPORTTHRESHOLD
    //<< endl; 	*metaOutputHandle << "#Using (-F 0x600 -f 0x001)" << endl;
    //	*metaOutputHandle << "#done\t" << printedBps << " lines printed in " <<
    //elapsedTime.quot << " minutes, " << elapsedTime.rem << " seconds" << endl;
    // }

} /* namespace sophia */
