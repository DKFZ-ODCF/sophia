/*
 * Alignment.cpp
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

#include "global.h"
#include "Alignment.h"
#include "GlobalAppConfig.h"
#include "HelperFunctions.h"
#include "MateInfo.h"
#include "Sdust.h"
#include "strtk-wrap.h"
#include <bitset>
#include <iostream>
#include <boost/exception/all.hpp>
#include <boost/lexical_cast.hpp>

namespace sophia {

    using namespace std;

    ChrSize Alignment::LOW_QUAL_CLIP_THRESHOLD{};

    int Alignment::BASE_QUALITY_THRESHOLD{},
        Alignment::BASE_QUALITY_THRESHOLD_LOW{};

    ChrSize Alignment::CLIPPED_NUCLEOTIDE_COUNT_THRESHOLD{},
            Alignment::INDEL_NUCLEOTIDE_COUNT_THRESHOLD{};

    double Alignment::ISIZEMAX{};
    Alignment::Alignment()
        : lowMapq(false),
          nullMapq(true),
          distantMate(0),
          chosenBp(nullptr),
          chrIndex(0),
          readType(0),
          startPos(0),
          endPos(0),
          mateChrIndex(0),
          matePos(0),
          samLine(),
          validLine(error_terminating_getline(cin, samLine)),
          samTabPositions(),
          saCbegin(),
          saCend(),
          hasSa(false),
          supplementary(false),
          fwdStrand(true),
          invertedMate(false),
          qualChecked(false) {

        if (validLine) {
            unsigned int index = 0;
            for (auto it = samLine.cbegin(); it != samLine.cend(); ++it) {
                if (*it == '\t') {
                    samTabPositions.push_back(index);
                }
                ++index;
            }
            try {
                chrIndex = GlobalAppConfig::getInstance().getChrConverter().parseChrAndReturnIndex(
                    next(samLine.cbegin(), static_cast<long>(samTabPositions[1]) + 1),
                    samLine.cend(),
                    '\t');
            } catch (DomainError &e) {
                e << error_info_string("line = " +
                                       std::string(next(samLine.cbegin(),
                                                   static_cast<long>(samTabPositions[1]) + 1),
                                                   samLine.cend()));
                throw e;
            }
        }
    }

    void
    Alignment::continueConstruction() {
        mappingQualityCheck();  // May set the readType to 7!
        for (auto startPos_cit = samLine.cbegin() + 1 + static_cast<int>(samTabPositions[2]);
             startPos_cit != samLine.cbegin() + static_cast<int>(samTabPositions[3]);
             ++startPos_cit) {
            startPos = startPos * 10 + ChrSize(*startPos_cit - '0');
        }
        ChrSize readLength = static_cast<ChrSize>(samTabPositions[9] - samTabPositions[8] - 1);
        if (readLength < 0) {
            throw_with_trace(std::logic_error("Invalid calculated readLength < 0: " +
                                              std::to_string(readLength)));
        }
        endPos = startPos + ChrSize(readLength);

        unsigned short flag = 0;
        for (auto flag_cit = samLine.cbegin() + 1 + static_cast<long>(samTabPositions[0]);
             flag_cit != samLine.cbegin() + static_cast<long>(samTabPositions[1]); ++flag_cit) {
             if (*flag_cit >= '0') {
                flag = flag * 10 + (unsigned short) ((signed short) *flag_cit - '0');
            } else {
                throw_with_trace(std::logic_error("Invalid flag in SAM file: " + samLine));
            }
        }
        auto flags = bitset<12>(flag);
        supplementary = (flags[11] == true);
        fwdStrand = (flags[4] == false);
        auto mateFwdStrand = (flags[5] == false);
        invertedMate = (fwdStrand == mateFwdStrand);

        // Alignment ends in match, or read contains soft-clip, hard-clip, insertion, or deletion.
        bool eventCandidate = isEventCandidate();
        if (eventCandidate) {
            createCigarChunks();
            assignBreakpointsAndOverhangs();
            if (supplementary) {
                auto startCit = next(samLine.cbegin(), 1 + static_cast<long>(samTabPositions[9]));
                auto endCit = next(samLine.cbegin(), static_cast<long>(samTabPositions[10]));
                vector<int> overhangPerBaseQuality{};
                fullMedianQuality(startCit, endCit, overhangPerBaseQuality);
                if (overhangPerBaseQuality.empty() ||
                    getMedian(overhangPerBaseQuality.begin(),
                              overhangPerBaseQuality.end()) <
                        BASE_QUALITY_THRESHOLD) {
                    // eventCandidate, supplementary, medianOverhangPerBaseQualities < BASQUALITYTHRESHOLD
                    readType = 5;
                } else {
                    // eventCandidate, supplementary, medianOverhangPerBaseQualities >= BASQUALITYTHRESHOLD
                    readType = 2;
                }
            }
            if (readType == 7) { // mapq != 0 && mapq < 13
                if (supplementary) {
                    if (uniqueSuppCheck() && hasSa) {
                        // eventCandidate, lowMapqCheckFailed, supplementary, uniqueSuppCheck, hasSa
                        readType = 2;
                    } else {
                        // eventCandidate, lowMapqCheckFailed, supplementary, (!uniqueSuppCheck || !hasSa)
                        readType = 5;
                    }
                } else {
                    readType = 5; // eventCandidate, lowMapqCheckFailed, !supplementary
                    auto rescueCandidate = false;
                    for (const auto &cigarChunk : cigarChunks) {
                        if (cigarChunk.chunkType == 'S') {
                            auto medianQual = overhangMedianQuality(cigarChunk);
                            if (cigarChunk.length > LOW_QUAL_CLIP_THRESHOLD &&
                                medianQual < BASE_QUALITY_THRESHOLD) {
                                rescueCandidate = false;
                                break;
                            }
                            if (cigarChunk.length / (readLength + 0.0) > 0.5) {
                                if (medianQual >= BASE_QUALITY_THRESHOLD) {
                                    rescueCandidate = true;
                                }
                            }
                        }
                    }
                    if (rescueCandidate) {
                        // eventCandidate, lowMapqCheckFailed, !supplementary, rescueCandidate
                        readType = 1;
                    }
                    qualChecked = true;
                }
            }
            if (readType < 5) {
                // Note that here readType is used as ordinal. It might be a score, ...?
                qualityCheckCascade();
            }
        } else if (readType == 7) {
            // !eventCandidate, mapq != 0 && mapq < 13
            readType = 5;
        }

        switch (readType) {
        case 0:
        case 3:
        case 5:
            assessOutlierMateDistance();
            if (distantMate == 1 && readType != 5) {
                readType = 4;
            }
            break;
        default:
            break;
        }

        for (auto mpos_cit = samLine.cbegin() + 1 + static_cast<long>(samTabPositions[6]);
             mpos_cit != samLine.cbegin() + static_cast<long>(samTabPositions[7]); ++mpos_cit) {
            matePos = matePos * 10 + ChrSize(*mpos_cit - '0');
        }
        if (samLine[1 + samTabPositions[5]] == '=') {
            mateChrIndex = chrIndex;
        } else {
            try {
                mateChrIndex = GlobalAppConfig::getInstance().getChrConverter().
                    parseChrAndReturnIndex(
                        next(samLine.cbegin(), 1 + static_cast<long>(samTabPositions[5])),
                        samLine.cend(),
                        '\t');
            } catch (const DomainError &e) {
                throw e << error_info_string(
                    "from = " + std::string(next(samLine.cbegin(),
                                                 1 + static_cast<long>(samTabPositions[5])),
                                            samLine.cend()));
            }
        }
    }

    void
    Alignment::mappingQualityCheck() {
        int mapq = boost::lexical_cast<int>(
            samLine.substr(samTabPositions[3] + 1,
                           samTabPositions[4] - samTabPositions[3] - 1));
        if (mapq == 0) {
            nullMapq = true;
            // readType = 0; lowMapq = true; see constructor
        } else {
            nullMapq = false;

            if (mapq < 13) {
                readType = 7;
                lowMapq = true;
            }
        }
    }

    /** The `Alignment` isEventCandidate` is true, if the last CIGAR code indicates a match,
     *  or if the CIGAR indicates a soft-clip, hard-clip, insertion, or deletion.
     */
    bool
    Alignment::isEventCandidate() const {
        // samTabPositions[0] is the position of the first tabulator. Zero-based index of the
        // CIGAR column in SAM is 5. Therefore, this means: The CIGAR string ends with a match.
        if (samLine[samTabPositions[5] - 1] != 'M') {
            return true;
        } else {
            // If the CIGAR does not end with a match, then (continue parsing the CIGAR string).
            // Return true, if there is a soft-clip, hard-clip, insertion, or deletion.
            for (auto cigarString_it = samLine.cbegin() + static_cast<long>(samTabPositions[4]) + 1;
                 cigarString_it != samLine.cbegin() + static_cast<long>(samTabPositions[5]) - 1;
                 ++cigarString_it) {
                switch (*cigarString_it) {
                case 'S': // soft-clipped
                case 'H': // hard-clipped
                case 'I': // insertion
                case 'D': // deletion
                    // The CIGAR string starts with a soft-clip, hard-clip, insertion or deletion.
                    return true;
                default:
                    // Continue with the next CIGAR code.
                    break;
                }
            }
            return false;
        }
    }

    void
    Alignment::createCigarChunks() {
        auto encounteredM = false;
        auto cumulativeNucleotideCount = 0,
             currentNucleotideCount = 0,
             indelAdjustment = 0,
             leftClipAdjustment = 0,
             rightClipAdjustment = 0;
        for (auto cigarString_cit = samLine.cbegin() + 1 + static_cast<long>(samTabPositions[4]);
             cigarString_cit != samLine.cbegin() + static_cast<long>(samTabPositions[5]);
             ++cigarString_cit) {
            if (isdigit(*cigarString_cit)) {
                currentNucleotideCount =
                    currentNucleotideCount * 10 + (*cigarString_cit - '0');
            } else {
                switch (*cigarString_cit) {
                case 'M':
                    encounteredM = true;
                    cumulativeNucleotideCount += currentNucleotideCount;
                    break;
                case 'S':
                    cigarChunks.emplace_back(
                        *cigarString_cit,
                         encounteredM,
                         cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
                         currentNucleotideCount,
                         indelAdjustment - leftClipAdjustment);
                    cumulativeNucleotideCount += currentNucleotideCount;
                    if (!encounteredM) {
                        leftClipAdjustment = currentNucleotideCount;
                    } else {
                        rightClipAdjustment = currentNucleotideCount;
                    }
                    break;
                case 'H':
                    cigarChunks.emplace_back(
                        *cigarString_cit, encounteredM,
                        cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
                        currentNucleotideCount);
                    break;
                case 'I':
                    cigarChunks.emplace_back(
                        *cigarString_cit,
                        encounteredM,
                        cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
                        currentNucleotideCount);
                    cumulativeNucleotideCount += currentNucleotideCount;
                    indelAdjustment -= currentNucleotideCount;
                    break;
                case 'D':
                    cigarChunks.emplace_back(
                        *cigarString_cit,
                        encounteredM,
                        cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
                        currentNucleotideCount);
                    indelAdjustment += currentNucleotideCount;
                    break;
                default:
                    break;
                }
                currentNucleotideCount = 0;
            }
        }
        endPos += ChrSize(indelAdjustment - leftClipAdjustment - rightClipAdjustment);
    }

    void
    Alignment::assignBreakpointsAndOverhangs() {
        for (const auto &chunk : cigarChunks) {
            switch (chunk.chunkType) {
            case 'S':
                readBreakpointTypes.push_back(chunk.chunkType);
                readBreakpointSizes.push_back(chunk.length);
                readBreakpointsEncounteredM.push_back(chunk.encounteredM);
                if (chunk.encounteredM) {
                    readBreakpoints.push_back(endPos);
                    readOverhangCoords.emplace_back(
                        chunk.encounteredM,
                        endPos,
                        static_cast<long>(chunk.startPosOnRead) - static_cast<long>(chunk.indelAdjustment),
                        chunk.length);
                } else {
                    readBreakpoints.push_back(startPos);
                    readOverhangCoords.emplace_back(
                        chunk.encounteredM,
                        startPos,
                        static_cast<long>(chunk.startPosOnRead) - static_cast<long>(chunk.indelAdjustment),
                        chunk.length);
                }
                break;
            case 'H':
                if (chunk.encounteredM) {
                    readBreakpoints.push_back(endPos);
                } else {
                    readBreakpoints.push_back(startPos);
                }
                readBreakpointSizes.push_back(chunk.length);
                readBreakpointTypes.push_back(chunk.chunkType);
                readBreakpointsEncounteredM.push_back(chunk.encounteredM);
                break;
            case 'I':
                readBreakpoints.push_back(startPos + ChrSize(chunk.startPosOnRead));
                readBreakpointTypes.push_back(chunk.chunkType);
                readBreakpointSizes.push_back(chunk.length);
                readBreakpointsEncounteredM.push_back(chunk.encounteredM);
                break;
            case 'D':
                readBreakpoints.push_back(startPos + ChrSize(chunk.startPosOnRead));
                readBreakpointTypes.push_back(chunk.chunkType);
                readBreakpointSizes.push_back(chunk.length);
                readBreakpointsEncounteredM.push_back(chunk.encounteredM);
                readBreakpoints.push_back(startPos +
                                          chunk.startPosOnRead +
                                          chunk.length);
                readBreakpointSizes.push_back(-1);
                readBreakpointTypes.push_back('#');
                readBreakpointsEncounteredM.push_back(chunk.encounteredM);
                break;
            default:
                break;
            }
        }
    }

    void
    Alignment::qualityCheckCascade() {
        // cerr << "a\n";
        if (!clipCountCheck()) {
            readType = 5;
            return;
        }
        // cerr << "b\n";
        if (!uniqueSuppCheck()) {
            readType = 5;
            return;
        }
        // cerr << "c\n";
        if (!qualChecked) {
            for (const auto &cigarChunk : cigarChunks) {
                if (cigarChunk.chunkType == 'S' &&
                    cigarChunk.length > LOW_QUAL_CLIP_THRESHOLD &&
                    overhangMedianQuality(cigarChunk) < BASE_QUALITY_THRESHOLD) {
                    readType = 5;
                    return;
                }
            }
        }
        // cerr << "d\n";
        assessReadType();
        // cerr << "e\n";
    }

    bool
    Alignment::clipCountCheck() {
        auto hCounts = 0;
        auto sCounts = 0;
        for (const auto &cigarChunk : cigarChunks) {
            switch (cigarChunk.chunkType) {
            case 'H':
                ++hCounts;
                break;
            case 'S':
                ++sCounts;
                break;
            default:
                break;
            }
        }
        if (hCounts + sCounts > 1 && nullMapq) {
            // mapq 0 is treated as a special case, where number of SAs and
            // base qualities will be the sole determinants of read quality
            lowMapq = true;
            return false;
        }
        return (hCounts < 2 && !(hCounts > 0 && sCounts > 0));
    }

    bool
    Alignment::uniqueSuppCheck() {
        auto hCounts = 0, sCounts = 0;
        for (const auto &cigarChunk : cigarChunks) {
            switch (cigarChunk.chunkType) {
            case 'H':
                ++hCounts;
                break;
            case 'S':
                ++sCounts;
                break;
            default:
                break;
            }
        }
        saCbegin = samLine.cend();
        saCend = samLine.cend();
        if (samLine.back() == ';' && samLine[samTabPositions.back() + 1] == 'S' &&
            samLine[samTabPositions.back() + 2] == 'A') {
            saCbegin = samLine.cbegin() + static_cast<long>(samTabPositions.back()) + 6;
            saCend = samLine.cend() - 1;
            hasSa = true;
        } else {
            for (auto i = 10u; i < samTabPositions.size() - 1; ++i) {
                if (samLine[samTabPositions[i + 1] - 1] == ';' &&
                    samLine[samTabPositions[i] + 1] == 'S' &&
                    samLine[samTabPositions[i] + 2] == 'A') {
                    saCbegin = samLine.cbegin() + static_cast<long>(samTabPositions[i]) + 6;
                    saCend = samLine.cbegin() + static_cast<long>(samTabPositions[i + 1]) - 1;
                    hasSa = true;
                    break;
                }
            }
        }
        if (hasSa) {
            auto lowQualSacounts = 0;
            auto block = 0;
            auto mapq = 0;
            auto highQualSa = false;
            for (auto saCit = saCbegin; saCit != saCend; ++saCit) {
                //"SA:Z:10,24753146,+,68S33M,48,1;X,135742083,-,47S22M32S,0,0;8,72637925,-,29S19M53S,0,0;"
                // 0-chr,1-pos,2-orientation,3-cigar,4-mapq,5-whatever
                switch (*saCit) {
                case ',':
                    ++block;
                    if (block == 4) {
                        ++saCit;
                        while (*saCit != ',') {
                            mapq = mapq * 10 + (*saCit - '0');
                            ++saCit;
                        }
                    }
                    break;
                case ';':
                    if (mapq < 13) {
                        ++lowQualSacounts;
                    } else if (mapq > 20) {
                        highQualSa = true;
                    }
                    if (!highQualSa && ((sCounts == 1 && lowQualSacounts == 2) ||
                                        (hCounts == 1 && lowQualSacounts == 2) ||
                                        (sCounts == 2 && lowQualSacounts == 4))) {
                        return false;
                    }
                    block = 0;
                    mapq = 0;
                    break;
                default:
                    break;
                }
            }
        }
        return true;
    }

    double
    Alignment::overhangMedianQuality(const CigarChunk &cigarChunk) const {
        vector<int> overhangPerBaseQuality{};
        if (!cigarChunk.encounteredM) {
            auto startCit = next(
                samLine.cbegin(),
                1 + static_cast<long>(samTabPositions[9]) +
                (static_cast<long>(cigarChunk.startPosOnRead) - static_cast<long>(cigarChunk.indelAdjustment)));
            auto endCit =
                next(samLine.cbegin(),
                     1 + static_cast<long>(samTabPositions[9]) +
                     (static_cast<long>(cigarChunk.startPosOnRead) - static_cast<long>(cigarChunk.indelAdjustment)) +
                     static_cast<long>(cigarChunk.length));
            fullMedianQuality(startCit, endCit, overhangPerBaseQuality);
        } else {
            string::const_reverse_iterator startCrit{
                next(samLine.cbegin(),
                     1 + static_cast<long>(samTabPositions[9]) +
                     (static_cast<long>(cigarChunk.startPosOnRead) - static_cast<long>(cigarChunk.indelAdjustment)) +
                     static_cast<long>(cigarChunk.length))
             };   // dito
            string::const_reverse_iterator endCrit{
                next(samLine.cbegin(),
                     1 + static_cast<long>(samTabPositions[9]) +
                     (static_cast<long>(cigarChunk.startPosOnRead) - static_cast<long>(cigarChunk.indelAdjustment)))};
            fullMedianQuality(startCrit, endCrit, overhangPerBaseQuality);
        }
        if (overhangPerBaseQuality.empty()) {
            return -1.0;
        } else {
            return getMedian(overhangPerBaseQuality.begin(),
                             overhangPerBaseQuality.end());
        }
    }

    void
    Alignment::assessReadType() {
        /* 0 for non-split,
         * 1 for softSplit,
         * 2 for hardSplit,
         * 3 for indel,
         * 4 for distant mate,
         * 5 for low quality overhang
         * 6 for low quality hardClipped
         * (precedence:
         *  decoy mate>
         *  low qual >
         *  soft clips >
         *  hard clips >
         *  distant mate >
         *  indels >
         *  normal)
         */
        auto hardLongClip = false, indelStatus = false;
        for (const auto &chunk : cigarChunks) {
            switch (chunk.chunkType) {
            case 'S':
                if (static_cast<ChrSize>(chunk.length) >= CLIPPED_NUCLEOTIDE_COUNT_THRESHOLD) {
                    readType = 1;
                    return;
                }
                break;
            case 'H':
                if (static_cast<ChrSize>(chunk.length) >= CLIPPED_NUCLEOTIDE_COUNT_THRESHOLD) {
                    hardLongClip = true;
                }
                break;
            case 'I':
            case 'D':
                indelStatus = true;
                break;
            default:
                break;
            }
        }
        if (hardLongClip) {
            readType = 2;
        } else if (indelStatus) {
            readType = 3;
        }
    }

    bool
    Alignment::assessOutlierMateDistance() {
        switch (distantMate) {
        case -1:
            return false;
        case 1:
            return true;
        default:
            if (*(samLine.cbegin() + 1 + static_cast<long>(samTabPositions[5])) != '=') {
                distantMate = 1;
                return true;
            } else {
                auto isize_cit = samLine.cbegin() + 1 + static_cast<long>(samTabPositions[7]);
                if (*isize_cit == '-') {
                    ++isize_cit;
                }
                auto isize = 0;
                for (; isize_cit != samLine.cbegin() + static_cast<long>(samTabPositions[8]);
                     ++isize_cit) {
                    isize = isize * 10 + (*isize_cit - '0');
                }
                if (isize > ISIZEMAX) {
                    distantMate = 1;
                    return true;
                }
            }
            distantMate = -1;
            return false;
        }
    }

    void
    Alignment::setChosenBp(ChrSize chosenBpLoc, int alignmentIndex) {
        auto overhangStartIndex = 0;
        ChrSize overhangLength = 0;
        char bpType{};
        auto bpEncounteredM = false;
        auto bpSize = 0;
        for (auto i = 0u; i < readBreakpoints.size(); ++i) {
            if (readBreakpoints[i] == chosenBpLoc) {
                bpEncounteredM = readBreakpointsEncounteredM[i];
                bpType = readBreakpointTypes[i];
                if (bpType == 'S') {
                    for (const auto &overhang : readOverhangCoords) {
                        if (overhang.bpPos == chosenBpLoc) {
                            overhangStartIndex = 1 + static_cast<int>(samTabPositions[8]) + overhang.startPosOnRead;
                            overhangLength = overhang.length;
                            break;
                        }
                    }
                }
                bpSize = readBreakpointSizes[i];
                break;
            }
        }
        chosenBp.reset();
        chosenBp = make_unique<ChosenBp>(ChosenBp(bpType,
                                                  bpSize,
                                                  bpEncounteredM,
                                                  overhangStartIndex,
                                                  overhangLength,
                                                  alignmentIndex /* origin index */));
    }

    vector<SuppAlignment>
    Alignment::generateSuppAlignments(ChrIndex bpChrIndex, int bpPos) {
        vector<SuppAlignment> suppAlignmentsTmp;
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();

        if (hasSa) {
            vector<string::const_iterator> saBegins = {saCbegin};
            vector<string::const_iterator> saEnds;
            for (auto it = saCbegin; it != saCend; ++it) {
                if (*it == ';') {
                    saEnds.push_back(it);
                    saBegins.push_back(it + 1);
                }
            }
            saEnds.push_back(saCend);
            for (auto i = 0u; i < saBegins.size(); ++i) {
                SuppAlignment saTmp = SuppAlignment::parseSamSaTag(
                                    saBegins[i],
                                    saEnds[i],
                                    !supplementary,
                                    lowMapq,
                                    nullMapq,
                                    fwdStrand,
                                    chosenBp->bpEncounteredM,
                                    chosenBp->selfNodeIndex,
                                    bpChrIndex,
                                    bpPos);
                if (!chrConverter.isTechnical(saTmp.getChrIndex())) {
                    suppAlignmentsTmp.push_back(saTmp);
                }
            }
        }
        if (assessOutlierMateDistance()) {
            if (!chrConverter.isTechnical(getMateChrIndex())) {
                auto foundMatch = false;
                MateInfo tmpPairDummy{
                    0, 0, getMateChrIndex(), getMatePos(), true, invertedMate};
                for (const auto &sa : suppAlignmentsTmp) {
                    if (tmpPairDummy.suppAlignmentFuzzyMatch(sa)) {
                        foundMatch = true;
                        break;
                    }
                }
                if (!foundMatch) {
                    suppAlignmentsTmp.emplace_back(SuppAlignment::create(
                        getMateChrIndex(),
                        getMatePos(),
                        0,
                        0,
                        chosenBp->bpEncounteredM,
                        invertedMate,
                        getMatePos() + 1,
                        !supplementary,
                        lowMapq,
                        nullMapq,
                        chosenBp->selfNodeIndex  /* origin index */));
                }
            }
        }
        return suppAlignmentsTmp;
    }

    string
    Alignment::printOverhang() const {
        string res{};
        res.reserve(static_cast<unsigned long>(chosenBp->overhangLength) + 9);
        if (chosenBp->bpEncounteredM) {
            res.append("|").append(samLine.substr(static_cast<unsigned long>(chosenBp->overhangStartIndex),
                                                  static_cast<unsigned long>(chosenBp->overhangLength)));
        } else {
            res.append(samLine.substr(static_cast<unsigned long>(chosenBp->overhangStartIndex),
                                      static_cast<unsigned long>(chosenBp->overhangLength)))
                .append("|");
        }
        res.append("(")
            .append(strtk::type_to_string<int>(chosenBp->childrenNodes.size()))
            .append(")");
        return res;
    }

    double
    Alignment::overhangComplexityMaskRatio() const {
        auto fullSizesTotal = 0.0;
        auto maskedIntervalsTotal = 0.0;
        vector<int> overhang;
        for (unsigned long i = 0; i < static_cast<unsigned long>(chosenBp->overhangLength); ++i) {
            switch (samLine[static_cast<unsigned long>(chosenBp->overhangStartIndex) + i]) {
            case 'A':
                overhang.push_back(0);
                break;
            case 'T':
                overhang.push_back(1);
                break;
            case 'G':
                overhang.push_back(2);
                break;
            case 'C':
                overhang.push_back(3);
                break;
            case 'N':
                if (!overhang.empty()) {
                    fullSizesTotal += overhang.size();
                    auto res = Sdust{overhang}.getRes();
                    if (!res.empty()) {
                        for (const auto &resInterval : res) {
                            maskedIntervalsTotal +=
                                resInterval.endIndex - resInterval.startIndex + 1;
                        }
                    }
                    overhang.clear();
                }
                break;
            default:
                break;
            }
        }
        if (!overhang.empty()) {
            fullSizesTotal += overhang.size();
            auto res = Sdust{overhang}.getRes();
            if (!res.empty()) {
                for (const auto &resInterval : res) {
                    maskedIntervalsTotal +=
                        resInterval.endIndex - resInterval.startIndex + 1;
                }
            }
        }
        return maskedIntervalsTotal / fullSizesTotal;
    }


    template <typename Iterator>
    void
    Alignment::fullMedianQuality(Iterator qualBegin, Iterator qualEnd,
                                 vector<int> &overhangPerBaseQuality) const {
        overhangPerBaseQuality.reserve((size_t) distance(qualBegin, qualEnd));
        auto consecutiveLowQuals = 0;
        for (auto cit = qualBegin; cit != qualEnd; ++cit) {
            if (*cit < BASE_QUALITY_THRESHOLD_LOW) {   // 33 + phred 11
                if (consecutiveLowQuals == 5) {
                    overhangPerBaseQuality.clear();
                    return;
                }
                ++consecutiveLowQuals;
            } else {
                consecutiveLowQuals = 0;
            }
            overhangPerBaseQuality.push_back(*cit);
        }
    }

    // Median Code taken from http://rosettacode.org/wiki/Averages/Median#C.2B.2B
    template <typename Iterator>
    double
    Alignment::getMedian(Iterator begin, Iterator end) const {
        // this is middle for odd-length, and "upper-middle" for even length
        Iterator middle = begin + (end - begin) / 2;
        // This function runs in O(n) on average, according to the standard
        nth_element(begin, middle, end);
        if ((end - begin) % 2 != 0) {   // odd length
            return *middle;
        } else {   // even length
            // the "lower middle" is the max of the lower half
            Iterator lower_middle = max_element(begin, middle);
            return (*middle + *lower_middle) / 2.0;
        }
    }

} /* namespace sophia */
