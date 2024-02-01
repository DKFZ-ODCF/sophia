/*
 * MrefEntry.h
 *
 *  Created on: 27 Nov 2016
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

#ifndef BREAKPOINTREDUCED_H_
#define BREAKPOINTREDUCED_H_
#include "global.h"
#include "Breakpoint.h"
#include "GermlineMatch.h"
#include "MrefMatch.h"
#include "SuppAlignmentAnno.h"
#include <boost/format.hpp>
#include <iostream>
#include <string>

namespace sophia {

    using namespace std;

    class BreakpointReduced {

      public:

        static ChrSize DEFAULT_READ_LENGTH;

        static double CLONALITY_STRICT_LOW_THRESHOLD;

        static double ARTIFACT_FREQ_HIGH_THRESHOLD;

        static string PIDS_IN_MREF_STR;

        static boost::format doubleFormatter;

        BreakpointReduced(const Breakpoint &tmpBp, int lineIndexIn,
                          bool hasOverhangIn);

        BreakpointReduced(const SuppAlignmentAnno &sa,
                          const BreakpointReduced &emittingBp, bool fuzzySecondary);

        template <typename T> bool operator<(const T &rhs) const {
            return pos < rhs.getPos();
        }

        /** This seems to be a function used for sorting positions, that assumes that the ordering
          * of chromosomes produces a total order, i.e., e.g., that all positions on chr1 are
          * "smaller" than all positions on chr2 (if chr1 has a smaller index than chr2).
          *
          * I doubt there is a biological meaning in the ordering of the chromosomes in the
          * classic Hg37ChrConverter. */
        bool fullSmaller(const BreakpointReduced &rhs) const {
            if (chrIndex < rhs.getChrIndex()) {
                return true;
            }
            if (chrIndex > rhs.getChrIndex()) {
                return false;
            }
            return pos < rhs.getPos();
        }

        template <typename T> int distanceTo(const T &rhs) const {
            if (chrIndex != rhs.getChrIndex()) {
                // Just any big value?
                return 1000000;
            } else {
                // Effectively always >= 0.
                return abs(pos - rhs.getPos());
            }
        }

        template <typename T> int distanceToBp(const T &compIn) const {
            if (chrIndex == compIn.getChrIndex()) {
                return abs((int) pos - (int) compIn.getPos());
            } else {
                // Ups. -1 is used in < comparisons. Check usages, before refactoring this.
                return -1;
            }
        }

        ChrIndex getChrIndex() const { return chrIndex; }

        ChrSize getPos() const { return pos; }


        bool isToRemove() const { return toRemove; }

        void setToRemove(bool toRemove) { this->toRemove = toRemove; }

        void removeMarkedFuzzies();

        SuppAlignmentAnno *searchFuzzySa(const SuppAlignmentAnno &fuzzySa);

        int getBreaksShortIndel() const { return breaksShortIndel; }

        int getLeftCoverage() const { return leftCoverage; }

        int getLowQualBreaksHard() const { return lowQualBreaksHard; }

        int getLowQualBreaksSoft() const { return lowQualBreaksSoft; }

        int getLowQualSpansHard() const { return lowQualSpansHard; }

        int getLowQualSpansSoft() const { return lowQualSpansSoft; }

        int getMateSupport() const { return mateSupport; }

        int getNormalSpans() const { return normalSpans; }

        int getPairedBreaksHard() const { return pairedBreaksHard; }

        int getPairedBreaksSoft() const { return pairedBreaksSoft; }

        int getRepetitiveOverhangBreaks() const { return repetitiveOverhangBreaks; }

        int getRightCoverage() const { return rightCoverage; }

        const vector<SuppAlignmentAnno> &getSuppAlignments() const {
            return suppAlignments;
        }


        int getUnpairedBreaksHard() const { return unpairedBreaksHard; }

        int getUnpairedBreaksSoft() const { return unpairedBreaksSoft; }

        int getLineIndex() const { return lineIndex; }

        vector<SuppAlignmentAnno *> getSupplementsPtr() {
            vector<SuppAlignmentAnno *> res{};
            for (auto &sa : suppAlignments) {
                res.push_back(&sa);
            }
            return res;
        }

        bool closeToSupp(const SuppAlignmentAnno &compIn, ChrSize fuzziness) const {
            if (chrIndex == compIn.getChrIndex()) {
                if (compIn.isFuzzy()) {
                    fuzziness = ChrSize(2.5 * DEFAULT_READ_LENGTH);  // truncate
                    return ((long) pos - (long) fuzziness) <= (long) (compIn.getExtendedPos() + fuzziness) &&
                           ((long) compIn.getPos() - (long) fuzziness) <= (long) (pos + fuzziness);
                } else {
                    return abs((long) pos - (long) compIn.getPos()) <= (long) fuzziness;
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
                            // TODO Why here getExtendenPos(), but getPos() above?
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

        const MrefMatch &getMrefHits() const { return mrefHits; }

        void setMrefHits(MrefMatch mrefHits) { this->mrefHits = mrefHits; }

        void setGermlineInfo(GermlineMatch germlineInfo) {
            this->germlineInfo = germlineInfo;
        }

        bool testOverhangBasedCandidacy() const;

        string printOverhang(double germlineClonality, int numHits,
                             const string &overhang) const;

        const GermlineMatch &getGermlineInfo() const { return germlineInfo; }

        void addFileIndex(int fileIndex) {
            for (auto &sa : suppAlignments) {
                sa.addFileIndex(fileIndex);
            }
        }

        void complexRearrangementMateRatioRescue(bool encounteredM);

        bool hasOverhang;

        void addDummySa(const SuppAlignmentAnno &sa,
                        const BreakpointReduced &emittingBp);

        const SuppAlignmentAnno &getDummySa();

      private:
        bool toRemove;
        int lineIndex;
        ChrIndex chrIndex;
        ChrSize pos;
        int normalSpans, lowQualSpansSoft, lowQualSpansHard, unpairedBreaksSoft,
            unpairedBreaksHard, breaksShortIndel, lowQualBreaksSoft,
            lowQualBreaksHard, repetitiveOverhangBreaks;
        int pairedBreaksSoft, pairedBreaksHard;
        int mateSupport;
        int leftCoverage, rightCoverage;
        MrefMatch mrefHits;
        GermlineMatch germlineInfo;
        vector<SuppAlignmentAnno> suppAlignments;
    };

} /* namespace sophia */

#endif /* BREAKPOINTREDUCED_H_ */
