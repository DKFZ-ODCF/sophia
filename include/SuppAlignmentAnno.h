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

#ifndef SUPPALIGNMENTANNO_H_
#define SUPPALIGNMENTANNO_H_
#include "CigarChunk.h"
#include "SuppAlignment.h"
#include "global.h"
#include <algorithm>
#include <array>
#include <string>
#include <unordered_set>
#include <vector>

namespace sophia {

    using namespace std;

    /**
     * @brief The SuppAlignmentAnno class
     * This is similar to SuppAlignment.
     *
     * Note that this class is under size constraints, as it will be instantiated via MrefEntry
     * once for each genome position in MasterRefEntry.
     **/
    class SuppAlignmentAnno {

      public:

        SuppAlignmentAnno(const string &saStrIn);

        SuppAlignmentAnno(const SuppAlignment &saIn);

        SuppAlignmentAnno(const SuppAlignmentAnno &saAnnoIn);

        SuppAlignmentAnno(ChrIndex emittingBpChrIndex,
                          ChrSize emittingBpPos,
                          const SuppAlignmentAnno &saAnnoIn);

        ~SuppAlignmentAnno() = default;

        static double ISIZEMAX;

        static ChrSize DEFAULT_READ_LENGTH;

        string print() const;

        void extendSuppAlignment(ChrSize minPos, ChrSize maxPos) {
            pos = min(pos, minPos);
            extendedPos = max(extendedPos, maxPos);
        }

        bool saCloseness(const SuppAlignmentAnno &rhs, int fuzziness) const;

        bool saClosenessDirectional(const SuppAlignmentAnno &rhs,
                                    int fuzziness) const;

        void removeFuzziness(const SuppAlignmentAnno &sa) {
            pos = sa.getPos();
            extendedPos = pos;
            fuzzy = false;
            if (!distant && sa.isDistant()) {
                distant = true;
            }
        }

        ChrIndex getChrIndex() const { return chrIndex; }

        bool isEncounteredM() const { return encounteredM; }

        bool isInverted() const { return inverted; }

        int getMateSupport() const { return mateSupport; }

        void incrementMateSupport() { ++mateSupport; }

        void setMateSupport(int mateSupportIn) { mateSupport = mateSupportIn; }

        ChrSize getPos() const { return pos; }

        int getSupport() const { return support; }

        int getSecondarySupport() const { return secondarySupport; }

        bool isToRemove() const { return toRemove; }

        void setToRemove(bool toRemove) { this->toRemove = toRemove; }

        bool isSuspicious() const { return suspicious; }

        void setSuspicious(bool suspicious) { this->suspicious = suspicious; }

        bool isDistant() const { return distant; }

        void setExpectedDiscordants(int expectedDiscordants) {
            this->expectedDiscordants = expectedDiscordants;
        }

        int getExpectedDiscordants() const { return expectedDiscordants; }

        bool isFuzzy() const { return fuzzy; }

        bool isStrictFuzzy() const { return strictFuzzy; }

        ChrSize getExtendedPos() const { return extendedPos; }

        bool isSemiSuspicious() const { return semiSuspicious; }

        void setSemiSuspicious(bool semiSuspicious) {
            this->semiSuspicious = semiSuspicious;
        }

        void setFuzzy(bool fuzzy) { this->fuzzy = fuzzy; }

        bool isProperPairErrorProne() const { return properPairErrorProne; }

        bool isStrictFuzzyCandidate() const { return strictFuzzyCandidate; }

        void addSupportingIndices(const vector<int> &supportingIndicesIn) {
            supportingIndices.insert(supportingIndices.end(),
                                     supportingIndicesIn.cbegin(),
                                     supportingIndicesIn.cend());
        }

        const vector<int> &getSupportingIndices() const {
            return supportingIndices;
        }

        void mergeMrefSa(const SuppAlignmentAnno &mrefSa);

        void finalizeSupportingIndices();

        void mrefSaTransform(int fileIndex) {
            supportingIndices.clear();
            supportingIndices.push_back(fileIndex);
        }

        void mrefSaConsensus(const unordered_set<unsigned short> &fileIndices) {
            supportingIndices.clear();
            for (const auto &index : fileIndices) {
                supportingIndices.push_back(index);
            }
        }

        void addFileIndex(int fileIndex) { supportingIndices.push_back(fileIndex); }

        void setSecondarySupport(int secondarySupport) {
            this->secondarySupport = secondarySupport;
        }

        void setSupport(int support) { this->support = support; }

      private:
        ChrIndex chrIndex;
        ChrSize pos;
        ChrSize extendedPos;
        int support;
        int secondarySupport;
        int mateSupport;
        int expectedDiscordants;
        bool encounteredM;
        bool toRemove;
        bool inverted;
        bool fuzzy;
        bool strictFuzzy;
        bool strictFuzzyCandidate;
        bool distant;
        bool suspicious;
        bool semiSuspicious;
        bool properPairErrorProne;
        vector<int> supportingIndices;

        static const string STOP_CHARS;
        inline bool isStopChar(char c);

    };
} /* namespace sophia */

#endif /* SUPPALIGNMENTANNO_H_ */
