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

#include "Alignment.h"
#include "GlobalAppConfig.h"
#include "HelperFunctions.h"
#include "MateInfo.h"
#include "Sdust.h"
#include "strtk-wrap.h"
#include <bitset>
#include <iostream>

namespace sophia {

using namespace std;

int Alignment::LOWQUALCLIPTHRESHOLD{},
    Alignment::BASEQUALITYTHRESHOLD{},
    Alignment::BASEQUALITYTHRESHOLDLOW{},
    Alignment::CLIPPEDNUCLEOTIDECOUNTTHRESHOLD{},
    Alignment::INDELNUCLEOTIDECOUNTTHRESHOLD{};

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
      samChunkPositions(),
      saCbegin(),
      saCend(),
      hasSa(false),
      supplementary(false),
      fwdStrand(true),
      invertedMate(false),
      qualChecked(false) {
    if (validLine) {
        auto index = 0;
        for (auto it = samLine.cbegin(); it != samLine.cend(); ++it) {
            if (*it == '\t') {
                samChunkPositions.push_back(index);
            }
            ++index;
        }
        chrIndex = GlobalAppConfig::getInstance().getChrConverter().parseChrAndReturnIndex(
            next(samLine.cbegin(), samChunkPositions[1] + 1),
            samLine.cend(),
            '\t');
    }
}

void
Alignment::continueConstruction() {
    mappingQualityCheck();
    for (auto startPos_cit = samLine.cbegin() + 1 + samChunkPositions[2];
         startPos_cit != samLine.cbegin() + samChunkPositions[3];
         ++startPos_cit) {
        startPos = startPos * 10 + (*startPos_cit - '0');
    }
    auto readLength = (samChunkPositions[9] - samChunkPositions[8] - 1);
    endPos = startPos + readLength;

    auto flag = 0;
    for (auto flag_cit = samLine.cbegin() + 1 + samChunkPositions[0];
         flag_cit != samLine.cbegin() + samChunkPositions[1]; ++flag_cit) {
        flag = flag * 10 + (*flag_cit - '0');
    }
    auto flags = bitset<12>(flag);
    supplementary = (flags[11] == true);
    fwdStrand = (flags[4] == false);
    auto mateFwdStrand = (flags[5] == false);
    invertedMate = (fwdStrand == mateFwdStrand);

    bool eventCandidate = isEventCandidate();
    if (eventCandidate) {
        createCigarChunks();
        assignBreakpointsAndOverhangs();
        if (supplementary) {
            auto startCit = next(samLine.cbegin(), 1 + samChunkPositions[9]);
            auto endCit = next(samLine.cbegin(), samChunkPositions[10]);
            vector<int> overhangPerBaseQuality{};
            fullMedianQuality(startCit, endCit, overhangPerBaseQuality);
            if (overhangPerBaseQuality.empty() ||
                getMedian(overhangPerBaseQuality.begin(),
                          overhangPerBaseQuality.end()) <
                    BASEQUALITYTHRESHOLD) {
                readType = 5;
            } else {
                readType = 2;
            }
        }
        if (readType == 7) {
            if (supplementary) {
                if (uniqueSuppCheck() && hasSa) {
                    readType = 2;
                } else {
                    readType = 5;
                }
            } else {
                readType = 5;
                auto rescueCandidate = false;
                for (const auto &cigarChunk : cigarChunks) {
                    if (cigarChunk.chunkType == 'S') {
                        auto medianQual = overhangMedianQuality(cigarChunk);
                        if (cigarChunk.length > LOWQUALCLIPTHRESHOLD &&
                            medianQual < BASEQUALITYTHRESHOLD) {
                            rescueCandidate = false;
                            break;
                        }
                        if (cigarChunk.length / (readLength + 0.0) > 0.5) {
                            if (medianQual >= BASEQUALITYTHRESHOLD) {
                                rescueCandidate = true;
                            }
                        }
                    }
                }
                if (rescueCandidate) {
                    readType = 1;
                }
                qualChecked = true;
            }
        }
        if (readType < 5) {
            qualityCheckCascade();
        }
    } else if (readType == 7) {
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
    for (auto mpos_cit = samLine.cbegin() + 1 + samChunkPositions[6];
         mpos_cit != samLine.cbegin() + samChunkPositions[7]; ++mpos_cit) {
        matePos = matePos * 10 + (*mpos_cit - '0');
    }
    if (samLine[1 + samChunkPositions[5]] == '=') {
        mateChrIndex = chrIndex;
    } else {
        mateChrIndex = GlobalAppConfig::getInstance().getChrConverter().
            parseChrAndReturnIndex(
                next(samLine.cbegin(), 1 + samChunkPositions[5]),
                samLine.cend(),
                '\t');
    }
}

void
Alignment::mappingQualityCheck() {
    if (samLine[1 + samChunkPositions[3]] !=
        '0') {   // mapq 0 is treated as a special case, where number of SAs and
                 // base qualities will be the sole determinants of read quality
        nullMapq = false;
        switch (samChunkPositions[4] - samChunkPositions[3]) {
        case 2:   // this checks if the mapq is a single-digit number
            readType = 7;
            lowMapq = true;
            break;
        case 3:
            if (samLine[1 + samChunkPositions[3]] == '1' &&
                ((samLine[2 + samChunkPositions[3]] - '0') < 3)) {
                // this checks if the mapq is a two-digit number, and if the
                // first digit is a "1", that the second digit is less than 3
                readType = 7;
                lowMapq = true;
            }
            break;
        default:
            break;
        }
    }
}

bool
Alignment::isEventCandidate() const {
    if (samLine[samChunkPositions[5] - 1] != 'M') {
        return true;
    } else {
        for (auto cigarString_cit = samLine.cbegin() + 1 + samChunkPositions[4];
             cigarString_cit != samLine.cbegin() + samChunkPositions[5] - 1;
             ++cigarString_cit) {
            switch (*cigarString_cit) {
            case 'S':
            case 'H':
            case 'I':
            case 'D':
                return true;
            default:
                break;
            }
        }
        return false;
    }
}

void
Alignment::createCigarChunks() {
    auto encounteredM = false;
    auto cumulativeNucleotideCount = 0, currentNucleotideCount = 0,
         indelAdjustment = 0, leftClipAdjustment = 0, rightClipAdjustment = 0;
    for (auto cigarString_cit = samLine.cbegin() + 1 + samChunkPositions[4];
         cigarString_cit != samLine.cbegin() + samChunkPositions[5];
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
                cigarChunks.emplace_back(*cigarString_cit, encounteredM,
                                         cumulativeNucleotideCount +
                                             indelAdjustment -
                                             leftClipAdjustment,
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
                cigarChunks.emplace_back(*cigarString_cit, encounteredM,
                                         cumulativeNucleotideCount +
                                             indelAdjustment -
                                             leftClipAdjustment,
                                         currentNucleotideCount);
                break;
            case 'I':
                cigarChunks.emplace_back(*cigarString_cit, encounteredM,
                                         cumulativeNucleotideCount +
                                             indelAdjustment -
                                             leftClipAdjustment,
                                         currentNucleotideCount);
                cumulativeNucleotideCount += currentNucleotideCount;
                indelAdjustment -= currentNucleotideCount;
                break;
            case 'D':
                cigarChunks.emplace_back(*cigarString_cit, encounteredM,
                                         cumulativeNucleotideCount +
                                             indelAdjustment -
                                             leftClipAdjustment,
                                         currentNucleotideCount);
                indelAdjustment += currentNucleotideCount;
                break;
            default:
                break;
            }
            currentNucleotideCount = 0;
        }
    }
    endPos += indelAdjustment - leftClipAdjustment - rightClipAdjustment;
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
                    chunk.encounteredM, endPos,
                    chunk.startPosOnRead - chunk.indelAdjustment, chunk.length);
            } else {
                readBreakpoints.push_back(startPos);
                readOverhangCoords.emplace_back(
                    chunk.encounteredM, startPos,
                    chunk.startPosOnRead - chunk.indelAdjustment, chunk.length);
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
            readBreakpoints.push_back(startPos + chunk.startPosOnRead);
            readBreakpointTypes.push_back(chunk.chunkType);
            readBreakpointSizes.push_back(chunk.length);
            readBreakpointsEncounteredM.push_back(chunk.encounteredM);
            break;
        case 'D':
            readBreakpoints.push_back(startPos + chunk.startPosOnRead);
            readBreakpointTypes.push_back(chunk.chunkType);
            readBreakpointSizes.push_back(chunk.length);
            readBreakpointsEncounteredM.push_back(chunk.encounteredM);
            readBreakpoints.push_back(startPos + chunk.startPosOnRead +
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
                cigarChunk.length > LOWQUALCLIPTHRESHOLD &&
                overhangMedianQuality(cigarChunk) < BASEQUALITYTHRESHOLD) {
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
    if (samLine.back() == ';' && samLine[samChunkPositions.back() + 1] == 'S' &&
        samLine[samChunkPositions.back() + 2] == 'A') {
        saCbegin = samLine.cbegin() + samChunkPositions.back() + 6;
        saCend = samLine.cend() - 1;
        hasSa = true;
    } else {
        for (auto i = 10u; i < samChunkPositions.size() - 1; ++i) {
            if (samLine[samChunkPositions[i + 1] - 1] == ';' &&
                samLine[samChunkPositions[i] + 1] == 'S' &&
                samLine[samChunkPositions[i] + 2] == 'A') {
                saCbegin = samLine.cbegin() + samChunkPositions[i] + 6;
                saCend = samLine.cbegin() + samChunkPositions[i + 1] - 1;
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
        auto startCit = next(samLine.cbegin(), 1 + samChunkPositions[9] +
                                                   cigarChunk.startPosOnRead -
                                                   cigarChunk.indelAdjustment);
        auto endCit =
            next(samLine.cbegin(),
                 1 + samChunkPositions[9] + cigarChunk.startPosOnRead -
                     cigarChunk.indelAdjustment + cigarChunk.length);
        fullMedianQuality(startCit, endCit, overhangPerBaseQuality);
    } else {
        string::const_reverse_iterator startCrit{
            next(samLine.cbegin(),
                 1 + samChunkPositions[9] + cigarChunk.startPosOnRead -
                     cigarChunk.indelAdjustment + cigarChunk.length)};
        string::const_reverse_iterator endCrit{
            next(samLine.cbegin(), 1 + samChunkPositions[9] +
                                       cigarChunk.startPosOnRead -
                                       cigarChunk.indelAdjustment)};
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
            if (chunk.length >= CLIPPEDNUCLEOTIDECOUNTTHRESHOLD) {
                readType = 1;
                return;
            }
            break;
        case 'H':
            if (chunk.length >= CLIPPEDNUCLEOTIDECOUNTTHRESHOLD) {
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
        if (*(samLine.cbegin() + 1 + samChunkPositions[5]) != '=') {
            distantMate = 1;
            return true;
        } else {
            auto isize_cit = samLine.cbegin() + 1 + samChunkPositions[7];
            if (*isize_cit == '-') {
                ++isize_cit;
            }
            auto isize = 0;
            for (; isize_cit != samLine.cbegin() + samChunkPositions[8];
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
Alignment::setChosenBp(int chosenBpLoc, int alignmentIndex) {
    auto overhangStartIndex = 0;
    auto overhangLength = 0;
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
                        overhangStartIndex =
                            1 + samChunkPositions[8] + overhang.startPosOnRead;
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
    chosenBp = make_unique<ChosenBp>(bpType, bpSize, bpEncounteredM,
                                     overhangStartIndex, overhangLength,
                                     alignmentIndex);
}
vector<SuppAlignment>
Alignment::generateSuppAlignments(int bpChrIndex, int bpPos) {
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
            SuppAlignment saTmp{saBegins[i],
                                saEnds[i],
                                !supplementary,
                                lowMapq,
                                nullMapq,
                                fwdStrand,
                                chosenBp->bpEncounteredM,
                                chosenBp->selfNodeIndex,
                                bpChrIndex,
                                bpPos};
            if (!chrConverter.isIgnoredChromosome(saTmp.getChrIndex())) {
                suppAlignmentsTmp.push_back(saTmp);
            }
        }
    }
    if (assessOutlierMateDistance()) {
        if (!chrConverter.isIgnoredChromosome(getMateChrIndex())) {
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
                suppAlignmentsTmp.emplace_back(
                    getMateChrIndex(), getMatePos(), 0, 0,
                    chosenBp->bpEncounteredM, invertedMate, getMatePos() + 1,
                    !supplementary, lowMapq, nullMapq, chosenBp->selfNodeIndex);
            }
        }
    }
    return suppAlignmentsTmp;
}
string
Alignment::printOverhang() const {
    string res{};
    res.reserve(chosenBp->overhangLength + 9);
    if (chosenBp->bpEncounteredM) {
        res.append("|").append(samLine.substr(chosenBp->overhangStartIndex,
                                              chosenBp->overhangLength));
    } else {
        res.append(samLine.substr(chosenBp->overhangStartIndex,
                                  chosenBp->overhangLength))
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
    for (auto i = 0; i < chosenBp->overhangLength; ++i) {
        switch (samLine[chosenBp->overhangStartIndex + i]) {
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
    overhangPerBaseQuality.reserve(distance(qualBegin, qualEnd));
    auto consecutiveLowQuals = 0;
    for (auto cit = qualBegin; cit != qualEnd; ++cit) {
        if (*cit < BASEQUALITYTHRESHOLDLOW) {   // 33 + phred 11
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
