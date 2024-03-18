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

        /** This is used for sorting breakpoints. No biological meaning. */
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
                return abs(static_cast<int>(pos) - static_cast<int>(compIn.getPos()));
            } else {
                // This seems to be a special value. It is not explicitly used in comparisons.
                // Check usages, before refactoring this.
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

        bool closeToSupp(const SuppAlignmentAnno &compIn, ChrDistance fuzziness) const {
            if (chrIndex == compIn.getChrIndex()) {
                if (compIn.isFuzzy()) {
                    fuzziness = ChrDistance(2.5 * DEFAULT_READ_LENGTH);  // truncate
                    return (ChrDistance(pos) - fuzziness) <= (ChrDistance(compIn.getExtendedPos()) + fuzziness) &&
                           (ChrDistance(compIn.getPos()) - fuzziness) <= (ChrDistance(pos) + fuzziness);
                } else {
                    return abs(ChrDistance(pos) - ChrDistance(compIn.getPos())) <= fuzziness;
                }
            } else {
                return false;
            }
        }

        ChrDistance distanceToSupp(const SuppAlignmentAnno &compIn) const {
            ChrDistance result;
            if (chrIndex == compIn.getChrIndex()) {
                if (compIn.isFuzzy()) {
                    if (compIn.getPos() <= pos && pos <= compIn.getExtendedPos()) {
                        result = 0;
                    } else {
                        if (pos < compIn.getPos()) {
                            result = ChrDistance(compIn.getPos() - pos);
                        } else {
                            // TODO Why here getExtendenPos(), but getPos() above?
                            result = ChrDistance(pos - compIn.getExtendedPos());
                        }
                    }
                } else {
                    result = ChrDistance(abs(static_cast<long>(pos) - static_cast<long>(compIn.getPos())));
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
        int mateSupport;
        int leftCoverage,
            rightCoverage;
        MrefMatch mrefHits;
        GermlineMatch germlineInfo;
        vector<SuppAlignmentAnno> suppAlignments;
    };

} /* namespace sophia */

#endif /* BREAKPOINTREDUCED_H_ */
