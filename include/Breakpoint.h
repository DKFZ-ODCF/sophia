/*
 * Breakpoint.h
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

#ifndef BREAKPOINT_H_
#define BREAKPOINT_H_
#include "Alignment.h"
#include "MateInfo.h"
#include "SuppAlignment.h"
#include "SuppAlignmentAnno.h"
#include "global.h"
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace sophia {

    using namespace std;

    class Breakpoint {
      public:

        Breakpoint(ChrIndex chrIndexIn,
                   ChrSize posIn);

        static Breakpoint parse(const string &bpIn,
                                bool ignoreOverhang);

        ~Breakpoint() = default;

        static const int PERMISSIBLE_MISMATCHES = 2;

        static const int MAX_PERMISSIBLE_SOFTCLIPS = 2000;

        static const int MAX_PERMISSIBLE_HARDCLIPS = 2000;

        static const int MAX_PERMISSIBLE_LOW_MAPQ_HARDCLIPS = 50;

        static int BP_SUPPORT_THRESHOLD;

        static ChrSize DEFAULT_READ_LENGTH;

        static ChrSize DISCORDANT_LOW_QUAL_LEFT_RANGE;

        static ChrSize DISCORDANT_LOW_QUAL_RIGHT_RANGE;

        static double IMPROPER_PAIR_RATIO;

        static bool PROPER_PAIR_COMPENSATION_MODE;

        static int bpindex;

        static const string COLUMN_STR;

        void addSoftAlignment(shared_ptr<Alignment> alignmentIn);

        void addHardAlignment(shared_ptr<Alignment> alignmentIn);

        bool finalizeBreakpoint(
            const deque<MateInfo> &discordantAlignmentsPool,
            const deque<MateInfo> &discordantLowQualAlignmentsPool,
            const deque<MateInfo> &discordantAlignmentCandidatesPool);

        void setLeftCoverage(int leftCoverageIn) { leftCoverage = leftCoverageIn; }

        void setRightCoverage(int rightCoverageIn) {
            rightCoverage = rightCoverageIn;
        }

        void setLowQualBreaksSoft(int lowQualBreaksSoftIn) {
            lowQualBreaksSoft = lowQualBreaksSoftIn;
        }

        void setLowQualBreaksHard(int lowQualBreaksHardIn) {
            lowQualBreaksHard = lowQualBreaksHardIn;
        }

        void setLowQualSpansSoft(int lowQualSpansSoftIn) {
            lowQualSpansSoft = lowQualSpansSoftIn;
        }

        void setLowQualSpansHard(int lowQualSpansHardIn) {
            lowQualSpansHard = lowQualSpansHardIn;
        }

        void setNormalSpans(int normalSpansIn) { normalSpans = normalSpansIn; }

        void setUnpairedBreaksSoft(int unpairedBreaksSoftIn) {
            unpairedBreaksSoft = unpairedBreaksSoftIn;
        }

        void setUnpairedBreaksHard(int unpairedBreaksHardIn) {
            unpairedBreaksHard = unpairedBreaksHardIn;
        }

        void setBreaksShortIndel(int breaksShortIndelIn) {
            breaksShortIndel = breaksShortIndelIn;
        }

        bool isCovFinalized() const { return covFinalized; }

        void setCovFinalized(bool covFinalizedIn) { covFinalized = covFinalizedIn; }

        template <typename T> bool operator<(const T &rhs) const {
            if (chrIndex < rhs.getChrIndex())
                return true;
            if (chrIndex > rhs.getChrIndex())
                return false;
            return (pos < rhs.getPos());
        }

        bool closeToSupp(const SuppAlignment &compIn, ChrSize fuzziness) const {
            if (chrIndex == compIn.getChrIndex()) {
                if (compIn.isFuzzy()) {
                    fuzziness = ChrSize(2.5 * DEFAULT_READ_LENGTH);  // Truncates to lower integer.
                    return (pos - fuzziness) <=
                               (compIn.getExtendedPos() + fuzziness) &&
                           (compIn.getPos() - fuzziness) <= (pos + fuzziness);
                } else {
                    return (unsigned long) abs((long) pos - (long) compIn.getPos()) <= fuzziness;
                }
            } else {
                return false;
            }
        }

        ChrSize distanceToSupp(const SuppAlignmentAnno &compIn) const {
            ChrSize result;
            if (chrIndex == compIn.getChrIndex()) {
                if (compIn.isFuzzy()) {
                    if (compIn.getPos() <= pos && pos <= compIn.getExtendedPos()) {
                        result = 0;
                    } else {
                        if (pos < compIn.getPos()) {
                            result = ChrSize(compIn.getPos() - pos);
                        } else {
                            result = ChrSize(pos - compIn.getExtendedPos());
                        }
                    }
                } else {
                    result = ChrSize(abs((long) pos - (long) compIn.getPos()));
                }
            } else {
                result = 1000000;
            }
            return result;
        }

        template <typename T> int distanceToBp(const T &compIn) const {
            if (chrIndex == compIn.getChrIndex()) {
                return abs(pos - compIn.getPos());
            } else {
                return -1;
            }
        }

        bool operator==(const Breakpoint &rhs) const {
            return chrIndex == rhs.getChrIndex() && pos == rhs.getPos();
        }

        ChrIndex getChrIndex() const { return chrIndex; }

        ChrSize getPos() const { return pos; }

        bool isMissingInfoBp() const { return missingInfoBp; }

        const vector<SuppAlignment> &getDoubleSidedMatches() const {
            return doubleSidedMatches;
        }

        vector<SuppAlignment *> getDoubleSidedMatchesPtr() {
            vector<SuppAlignment *> res{};
            for (auto &sa : doubleSidedMatches) {
                res.push_back(&sa);
            }
            return res;
        }

        const vector<SuppAlignment> &getSupplementsPrimary() const {
            return supplementsPrimary;
        }

        vector<SuppAlignment *> getSupplementsPrimaryPtr() {
            vector<SuppAlignment *> res{};
            for (auto &sa : supplementsPrimary) {
                res.push_back(&sa);
            }
            return res;
        }

        bool isGermline() const { return germline; }

        int getHitsInMref() const { return hitsInMref; }

        int getLeftCoverage() const { return leftCoverage; }

        int getRightCoverage() const { return rightCoverage; }

        int getBreaksShortIndel() const { return breaksShortIndel; }

        const vector<string> &getConsensusOverhangs() const {
            return consensusOverhangs;
        }

        int getLowQualBreaksSoft() const { return lowQualBreaksSoft; }

        int getLowQualBreaksHard() const { return lowQualBreaksHard; }

        int getRepetitiveOverhangBreaks() const { return repetitiveOverhangBreaks; }

        int getLowQualSpansSoft() const { return lowQualSpansSoft; }

        int getLowQualSpansHard() const { return lowQualSpansHard; }

        int getMateSupport() const { return mateSupport; }

        int getNormalSpans() const { return normalSpans; }

        int getPairedBreaksSoft() const { return pairedBreaksSoft; }

        int getPairedBreaksHard() const { return pairedBreaksHard; }

        int getUnpairedBreaksHard() const { return unpairedBreaksHard; }

        int getUnpairedBreaksSoft() const { return unpairedBreaksSoft; }

        void removeMarkedFuzzies() {
            cleanUpVector(doubleSidedMatches);
            cleanUpVector(supplementsPrimary);
        }

        SuppAlignment *searchFuzzySa(const SuppAlignment &fuzzySa);

        void setGermline(bool germlineIn) { this->germline = germlineIn; }

        void setHitsInMref(int hitsInMref) { this->hitsInMref = hitsInMref; }

      private:

        // Compose the string that will be printed as column 8 into the breakpoint BED.
        string finalizeOverhangs();

        // Actually prints to stdout.
        void printBreakpointReport(const string &overhangStr);

        bool matchDetector(const shared_ptr<Alignment> &longAlignment,
                           const shared_ptr<Alignment> &shortAlignment) const;

        void detectDoubleSupportSupps();

        void collapseSuppRange(string &res, const vector<SuppAlignment> &vec) const;

        template <typename T> void cleanUpVector(vector<T> &objectPool);

        void fillMatePool(const deque<MateInfo> &discordantAlignmentsPool,
                          const deque<MateInfo> &discordantLowQualAlignmentsPool,
                          const deque<MateInfo> &discordantAlignmentCandidatesPool);

        void collectMateSupport();

        void compressMatePool(vector<MateInfo> &discordantAlignmentsPool);

        void
        collectMateSupportHelper(SuppAlignment &sa,
                                 vector<MateInfo> &discordantAlignmentsPool,
                                 vector<MateInfo> &discordantLowQualAlignmentsPool);

        void saHomologyClashSolver();

        bool covFinalized;

        bool missingInfoBp;

        ChrIndex chrIndex;

        ChrSize pos;

        int normalSpans,
            lowQualSpansSoft,
            lowQualSpansHard,
            unpairedBreaksSoft,
            unpairedBreaksHard,
            breaksShortIndel,
            lowQualBreaksSoft,
            lowQualBreaksHard,
            repetitiveOverhangBreaks;

        int pairedBreaksSoft,
            pairedBreaksHard;

        int leftSideDiscordantCandidates,
            rightSideDiscordantCandidates;

        int mateSupport;

        int leftCoverage,
            rightCoverage;

        int totalLowMapqHardClips;

        int hitsInMref;

        bool germline;

        vector<shared_ptr<Alignment>> supportingSoftAlignments;

        vector<shared_ptr<Alignment>> supportingHardAlignments;

        vector<shared_ptr<Alignment>> supportingHardLowMapqAlignments;

        vector<SuppAlignment> supplementsPrimary;

        vector<SuppAlignment> doubleSidedMatches;

        vector<string> consensusOverhangs;

        vector<MateInfo> poolLeft, poolRight, poolLowQualLeft, poolLowQualRight;

        vector<SuppAlignment> supplementsSecondary;
    };

} /* namespace sophia */

#endif /* BREAKPOINT_H_ */
