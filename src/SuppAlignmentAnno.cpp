/*
 * SuppAlignment.cpp
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

#include "SuppAlignmentAnno.h"
#include "GlobalAppConfig.h"
#include "ChrConverter.h"
#include "strtk.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
// #include <iostream>

namespace sophia {

using namespace std;

double SuppAlignmentAnno::ISIZEMAX{};
int SuppAlignmentAnno::DEFAULTREADLENGTH{};

SuppAlignmentAnno::SuppAlignmentAnno(const string &saStrIn)
    : chrIndex{0}, pos{0}, extendedPos{0}, support{0}, secondarySupport{0},
      mateSupport{0}, expectedDiscordants{0}, encounteredM{saStrIn[0] == '|'},
      toRemove{false}, inverted{false}, fuzzy{false}, strictFuzzy{false},
      strictFuzzyCandidate{false}, distant{false}, suspicious{false},
      semiSuspicious{false}, properPairErrorProne{saStrIn.back() == '#'},
      supportingIndices{} {
    auto index = 0;
    if (encounteredM) {
        ++index;
    }
    chrIndex =
        GlobalAppConfig::getInstance().
            getChrConverter().parseChrAndReturnIndex(next(saStrIn.cbegin(), index), ':');
    if (chrIndex > 1001) {
        return;
    }
    while (saStrIn[index] != ':') {
        ++index;
    }
    ++index;
    for (; saStrIn[index] != '('; ++index) {
        if (saStrIn[index] == '-') {
            fuzzy = true;
        } else if (saStrIn[index] == '_') {
            inverted = true;
            while (saStrIn[index] != '(') {
                ++index;
            }
            break;
        } else if (saStrIn[index] != '|') {
            if (!fuzzy) {
                pos = 10 * pos + (saStrIn[index] - '0');
            } else {
                extendedPos = 10 * extendedPos + (saStrIn[index] - '0');
            }
        }
    }
    if (!fuzzy) {
        extendedPos = pos;
    }
    ++index;
    for (; saStrIn[index] != ','; ++index) {
        support = 10 * support + (saStrIn[index] - '0');
    }
    ++index;
    for (; saStrIn[index] != ','; ++index) {
        secondarySupport = 10 * secondarySupport + (saStrIn[index] - '0');
    }
    ++index;
    if (saStrIn[index] == '!') {
        suspicious = true;
        index += 2;
    } else {
        for (; saStrIn[index] != '/'; ++index) {
            if (saStrIn[index] == '?') {
                semiSuspicious = true;
            } else {
                mateSupport = 10 * mateSupport + (saStrIn[index] - '0');
            }
        }
        ++index;
    }
    for (; saStrIn[index] != ')'; ++index) {
        expectedDiscordants = 10 * expectedDiscordants + (saStrIn[index] - '0');
    }
    if (support + secondarySupport == 0) {
        fuzzy = true;
    }
    distant = expectedDiscordants > 0 || suspicious;
    strictFuzzyCandidate = (support + secondarySupport) < 3;
    strictFuzzy = fuzzy || strictFuzzyCandidate;
}

SuppAlignmentAnno::SuppAlignmentAnno(const SuppAlignment &saIn)
    : chrIndex{saIn.getChrIndex()}, pos{saIn.getPos()},
      extendedPos{saIn.getExtendedPos()}, support{saIn.getSupport()},
      secondarySupport{saIn.getSecondarySupport()},
      mateSupport{saIn.getMateSupport()},
      expectedDiscordants{saIn.getExpectedDiscordants()},
      encounteredM{saIn.isEncounteredM()}, toRemove{false},
      inverted{saIn.isInverted()}, fuzzy{saIn.isFuzzy()}, strictFuzzy{false},
      strictFuzzyCandidate{false}, distant{false},
      suspicious{saIn.isSuspicious()}, semiSuspicious{saIn.isSemiSuspicious()},
      properPairErrorProne{saIn.isProperPairErrorProne()},
      supportingIndices{saIn.getSupportingIndices()} {
    distant = expectedDiscordants > 0 || suspicious;
    if (support + secondarySupport == 0) {
        fuzzy = true;
    }
    strictFuzzyCandidate = (support + secondarySupport) < 3;
    strictFuzzy = fuzzy || strictFuzzyCandidate;
}

SuppAlignmentAnno::SuppAlignmentAnno(const SuppAlignmentAnno &saAnnoIn)
    : chrIndex{saAnnoIn.getChrIndex()}, pos{saAnnoIn.getPos()},
      extendedPos{saAnnoIn.getExtendedPos()}, support{saAnnoIn.getSupport()},
      secondarySupport{saAnnoIn.getSecondarySupport()},
      mateSupport{saAnnoIn.getMateSupport()},
      expectedDiscordants{saAnnoIn.getExpectedDiscordants()},
      encounteredM{saAnnoIn.isEncounteredM()}, toRemove{false},
      inverted{saAnnoIn.isInverted()}, fuzzy{saAnnoIn.isFuzzy()},
      strictFuzzy{saAnnoIn.isStrictFuzzy()},
      strictFuzzyCandidate{saAnnoIn.isStrictFuzzyCandidate()},
      distant{saAnnoIn.isDistant()}, suspicious{saAnnoIn.isSuspicious()},
      semiSuspicious{saAnnoIn.isSemiSuspicious()},
      properPairErrorProne{saAnnoIn.isProperPairErrorProne()},
      supportingIndices{saAnnoIn.getSupportingIndices()} {}

SuppAlignmentAnno::SuppAlignmentAnno(int emittingBpChrIndex, int emittingBpPos,
                                     const SuppAlignmentAnno &saAnnoIn)
    : chrIndex{emittingBpChrIndex},
      pos{saAnnoIn.isDistant()
              ? max(1, static_cast<int>(
                           round(emittingBpPos - 1.5 * DEFAULTREADLENGTH)))
              : emittingBpPos},
      extendedPos{
          saAnnoIn.isDistant()
              ? static_cast<int>(emittingBpPos + 1.5 * DEFAULTREADLENGTH)
              : emittingBpPos},
      support{0}, secondarySupport{0}, mateSupport{0}, expectedDiscordants{0},
      encounteredM{true}, toRemove{true}, inverted{saAnnoIn.isInverted()},
      fuzzy{saAnnoIn.isDistant()}, strictFuzzy{saAnnoIn.isDistant()},
      strictFuzzyCandidate{true}, distant{saAnnoIn.isDistant()},
      suspicious{saAnnoIn.isSuspicious()},
      semiSuspicious{saAnnoIn.isSemiSuspicious()},
      properPairErrorProne{saAnnoIn.isProperPairErrorProne()},
      supportingIndices{} {}

string
SuppAlignmentAnno::print() const {
    string outStr;
    outStr.reserve(36);
    string invStr{};
    if (inverted) {
        invStr.append("_INV");
    }
    if (encounteredM) {
        outStr.append("|");
    } else {
        invStr.append("|");
    }
    const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
    if (!fuzzy || pos == extendedPos) {
        outStr.append(chrConverter.indexToChrName(chrIndex))
            .append(":")
            .append(strtk::type_to_string<int>(pos));
    } else {
        outStr.append(chrConverter.indexToChrName(chrIndex))
            .append(":")
            .append(strtk::type_to_string<int>(pos))
            .append("-")
            .append(strtk::type_to_string<int>(extendedPos));
    }
    outStr.append(invStr)
        .append("(")
        .append(strtk::type_to_string<int>(support))
        .append(",")
        .append(strtk::type_to_string<int>(secondarySupport))
        .append(",");
    if (!suspicious) {
        outStr.append(strtk::type_to_string<int>(mateSupport));
        if (semiSuspicious) {
            outStr.append("?");
        }
    } else {
        outStr.append("!");
    }
    outStr.append("/")
        .append(strtk::type_to_string<int>(expectedDiscordants))
        .append(")");
    if (properPairErrorProne) {
        outStr.append("#");
    }
    return outStr;
}

bool
SuppAlignmentAnno::saCloseness(const SuppAlignmentAnno &rhs,
                               int fuzziness) const {
    // inverted == rhs.isInverted()  &&   //test
    if (chrIndex == rhs.getChrIndex()) {
        if (strictFuzzy || rhs.isStrictFuzzy()) {
            if (pos <= rhs.getExtendedPos() && rhs.getPos() <= extendedPos) {
                return true;
            }
            fuzziness = 2.5 * DEFAULTREADLENGTH;
            if (rhs.getPos() >= extendedPos) {
                return (rhs.getPos() - extendedPos) <= fuzziness;
            }
            return (pos - rhs.getExtendedPos()) <= fuzziness;
            //			return (rhs.getPos() - fuzziness) <= (extendedPos +
            //fuzziness) && (pos - fuzziness) <= (rhs.getExtendedPos() +
            //fuzziness);
        } else {
            return abs(pos - rhs.getPos()) <= fuzziness;
        }
    } else {
        return false;
    }
}

bool
SuppAlignmentAnno::saClosenessDirectional(const SuppAlignmentAnno &rhs,
                                          int fuzziness) const {
    // inverted == rhs.isInverted()  &&   //test
    if (chrIndex == rhs.getChrIndex() && encounteredM == rhs.isEncounteredM()) {
        if (strictFuzzy || rhs.isStrictFuzzy()) {
            if (pos <= rhs.getExtendedPos() && rhs.getPos() <= extendedPos) {
                return true;
            }
            fuzziness = 2.5 * DEFAULTREADLENGTH;
            if (rhs.getPos() >= extendedPos) {
                return (rhs.getPos() - extendedPos) <= fuzziness;
            }
            return (pos - rhs.getExtendedPos()) <= fuzziness;
        } else {
            return abs(pos - rhs.getPos()) <= fuzziness;
        }
    } else {
        return false;
    }
}

void
SuppAlignmentAnno::mergeMrefSa(const SuppAlignmentAnno &mrefSa) {
    support = max(support, mrefSa.getSupport());
    secondarySupport = max(secondarySupport, mrefSa.getSecondarySupport());
    for (auto index : mrefSa.getSupportingIndices()) {
        supportingIndices.push_back(index);
    }
    sort(supportingIndices.begin(), supportingIndices.end());
    supportingIndices.erase(
        unique(supportingIndices.begin(), supportingIndices.end()),
        supportingIndices.end());
    if (mrefSa.getExpectedDiscordants() > 0 && expectedDiscordants > 0) {
        if ((0.0 + mrefSa.getMateSupport()) / mrefSa.getExpectedDiscordants() >
            (0.0 + mateSupport) / expectedDiscordants) {
            mateSupport = mrefSa.getMateSupport();
            expectedDiscordants = mrefSa.getExpectedDiscordants();
        }
    } else if (mrefSa.getExpectedDiscordants() > 0) {
        mateSupport = mrefSa.getMateSupport();
        expectedDiscordants = mrefSa.getExpectedDiscordants();
    }
    if (!mrefSa.isSemiSuspicious() && semiSuspicious) {
        semiSuspicious = false;
    }
}

void
SuppAlignmentAnno::finalizeSupportingIndices() {
    sort(supportingIndices.begin(), supportingIndices.end());
    supportingIndices.erase(
        unique(supportingIndices.begin(), supportingIndices.end()),
        supportingIndices.end());
    support = static_cast<int>(supportingIndices.size());
    secondarySupport = 0;
}

} /* namespace sophia */
