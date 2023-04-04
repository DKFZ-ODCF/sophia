/*
 * SuppAlignment.h
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

#ifndef SUPPALIGNMENT_H_
#define SUPPALIGNMENT_H_
#include "CigarChunk.h"
#include <algorithm>
#include <array>
#include <string>
#include <unordered_set>
#include <vector>

namespace sophia {

using namespace std;

class SuppAlignment {
  public:
    SuppAlignment(string::const_iterator saCbegin,
                  string::const_iterator saCend, bool primaryIn,
                  bool lowMapqSourceIn, bool nullMapqSourceIn,
                  bool alignmentOnForwardStrand, bool bpEncounteredM,
                  int originIndexIn, int bpChrIndex, int bpPos);
    SuppAlignment(int chrIndexIn, int posIn, int mateSupportIn,
                  int expectedDiscordantsIn, bool encounteredMIn,
                  bool invertedIn, int extendedPosIn, bool primaryIn,
                  bool lowMapqSourceIn, bool nullMapqSourceIn,
                  int originIndexIn);
    SuppAlignment(const string &saIn);
    ~SuppAlignment() = default;
    static double ISIZEMAX;
    static int DEFAULTREADLENGTH;
    string print() const;
    void extendSuppAlignment(int minPos, int maxPos) {
        pos = min(pos, minPos);
        extendedPos = max(extendedPos, maxPos);
    }
    bool saCloseness(const SuppAlignment &rhs, int fuzziness) const;
    bool saDistHomologyRescueCloseness(const SuppAlignment &rhs,
                                       int fuzziness) const;
    void padMateSupportHomologyRescue() { expectedDiscordants = mateSupport; }
    void removeFuzziness(const SuppAlignment &sa) {
        pos = sa.getPos();
        extendedPos = pos;
        fuzzy = false;
        if (!distant && sa.isDistant()) {
            distant = true;
        }
    }
    int getChrIndex() const { return chrIndex; }
    bool isEncounteredM() const { return encounteredM; }
    bool isInverted() const { return inverted; }
    int getMateSupport() const { return mateSupport; }
    void incrementDistinctReads() { ++distinctReads; }
    void incrementMateSupport(int incrementIn) { mateSupport += incrementIn; }
    void setMateSupport(int mateSupportIn) { mateSupport = mateSupportIn; }
    int getPos() const { return pos; }
    bool isPrimary() const { return primary; }
    int getSupport() const { return support; }
    void addSupportingIndices(const vector<int> &supportingIndicesIn) {
        supportingIndices.insert(supportingIndices.end(),
                                 supportingIndicesIn.cbegin(),
                                 supportingIndicesIn.cend());
    }
    void addSecondarySupportIndices(int supportingIndicesSecondaryIn) {
        supportingIndicesSecondary.push_back(supportingIndicesSecondaryIn);
    }
    void addSecondarySupportIndices(
        const vector<int> &supportingIndicesSecondaryIn) {
        supportingIndicesSecondary.insert(supportingIndicesSecondary.end(),
                                          supportingIndicesSecondaryIn.cbegin(),
                                          supportingIndicesSecondaryIn.cend());
    }
    void finalizeSupportingIndices();
    int getSecondarySupport() const { return secondarySupport; }
    bool isToRemove() const { return toRemove; }
    void setToRemove(bool toRemove) { this->toRemove = toRemove; }
    int getMapq() const { return mapq; }
    void setMapq(int mapq) { this->mapq = mapq; }
    bool isSuspicious() const { return suspicious; }
    void setSuspicious(bool suspicious) { this->suspicious = suspicious; }
    bool isDistant() const { return distant; }
    const vector<int> &getSupportingIndices() const {
        return supportingIndices;
    }
    const vector<int> &getSupportingIndicesSecondary() const {
        return supportingIndicesSecondary;
    }

    void setExpectedDiscordants(int expectedDiscordants) {
        this->expectedDiscordants = expectedDiscordants;
    }

    int getExpectedDiscordants() const { return expectedDiscordants; }

    int getDistinctReads() const { return distinctReads; }

    int getMatchFuzziness() const { return matchFuzziness; }

    bool isFuzzy() const { return fuzzy; }
    bool isStrictFuzzy() const { return strictFuzzy; }
    int getExtendedPos() const { return extendedPos; }

    bool isLowMapqSource() const { return lowMapqSource; }
    void mrefSaTransform(int fileIndex) {
        support = 0;
        secondarySupport = 0;
        supportingIndices.clear();
        supportingIndices.push_back(fileIndex);
    }
    void mrefSaConsensus(const unordered_set<short> &fileIndices) {
        supportingIndices.clear();
        for (const auto &index : fileIndices) {
            supportingIndices.push_back(index);
        }
    }
    void mergeSa(const SuppAlignment &rhs) {
        support = max(support, rhs.getSupport());
        secondarySupport = max(secondarySupport, rhs.getSecondarySupport());
        if (rhs.getExpectedDiscordants() > 0 && expectedDiscordants > 0) {
            if ((0.0 + rhs.getMateSupport()) / rhs.getExpectedDiscordants() >
                (0.0 + mateSupport) / expectedDiscordants) {
                mateSupport = rhs.getMateSupport();
                expectedDiscordants = rhs.getExpectedDiscordants();
            }
        } else if (rhs.getExpectedDiscordants() > 0) {
            mateSupport = rhs.getMateSupport();
            expectedDiscordants = rhs.getExpectedDiscordants();
        }
    }
    void mergeMrefSa(const SuppAlignment &mrefSa) {
        for (auto index : mrefSa.getSupportingIndices()) {
            supportingIndices.push_back(index);
        }
        sort(supportingIndices.begin(), supportingIndices.end());
        sort(supportingIndicesSecondary.begin(),
             supportingIndicesSecondary.end());
        supportingIndices.erase(
            unique(supportingIndices.begin(), supportingIndices.end()),
            supportingIndices.end());
        if (mrefSa.getExpectedDiscordants() > 0 && expectedDiscordants > 0) {
            if ((0.0 + mrefSa.getMateSupport()) /
                    mrefSa.getExpectedDiscordants() >
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

    bool isSemiSuspicious() const { return semiSuspicious; }

    void setSemiSuspicious(bool semiSuspicious) {
        this->semiSuspicious = semiSuspicious;
    }

    bool isNullMapqSource() const { return nullMapqSource; }

    void setNullMapqSource(bool nullMapqSource) {
        this->nullMapqSource = nullMapqSource;
    }

    bool isProperPairErrorProne() const { return properPairErrorProne; }

    void setProperPairErrorProne(bool properPairErrorProne) {
        this->properPairErrorProne = properPairErrorProne;
    }

  private:
    int matchFuzziness;
    int chrIndex;
    int pos;
    int extendedPos;
    int mapq;
    vector<int> supportingIndices;
    vector<int> supportingIndicesSecondary;
    int distinctReads;
    int support;
    int secondarySupport;
    int mateSupport;
    int expectedDiscordants;
    bool encounteredM;
    bool toRemove;
    bool inverted;
    bool fuzzy;
    bool strictFuzzy;
    bool distant;
    bool lowMapqSource;
    bool nullMapqSource;
    bool suspicious;
    bool semiSuspicious;
    bool properPairErrorProne;
    bool primary;
};
} /* namespace sophia */

#endif /* SUPPALIGNMENT_H_ */
