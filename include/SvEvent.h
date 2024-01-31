/*
 * PairInfo.h
 *
 *  Created on: 24 Oct 2016
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

#ifndef SVEVENT_H_
#define SVEVENT_H_

#include "global.h"
#include "Breakpoint.h"
#include "GermlineMatch.h"
#include "MrefMatch.h"
#include "SuppAlignmentAnno.h"
#include <BreakpointReduced.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace sophia {

    using namespace std;

    enum ArtifactStatus { ARTIFACT, BORDERLINE, CLEAN, AS_UNKNOWN };

    enum ClonalityStatus { HOMO, HETERO, SUBCLONAL, EXTREME_SUBCLONAL, CS_UNKNOWN };

    class SvEvent {
      public:
        static boost::format doubleFormatter;

        static int GERMLINE_OFFSET_THRESHOLD;

        static double RELAXED_BP_FREQ_THRESHOLD;

        static double BP_FREQ_THRESHOLD;

        static double ARTIFACT_FREQ_LOW_THRESHOLD;

        static double ARTIFACT_FREQ_HIGH_THRESHOLD;

        static double CLONALITY_LOW_THRESHOLD;

        static double CLONALITY_STRICT_LOW_THRESHOLD;

        static double CLONALITY_HIGH_THRESHOLD;

        static string PIDS_IN_MREF_STR;

        static int HALFDEFAULT_READ_LENGTH;

        static int GERMLINE_DB_LIMIT;

        static bool ABRIDGED_OUTPUT;

        static bool NO_CONTROL_MODE;

        static bool DEBUG_MODE;

        const static vector<string> EVENTTYPES;

        SvEvent(const BreakpointReduced &bp1In,
                const BreakpointReduced &bp2In,
                const SuppAlignmentAnno &sa1In,
                const SuppAlignmentAnno &sa2In,
                const vector<pair<int, string>> &overhangDb);

        SvEvent(const BreakpointReduced &bp1In,
                const BreakpointReduced &bp2In,
                const SuppAlignmentAnno &sa1In,
                const vector<pair<int, string>> &overhangDb,
                const SuppAlignmentAnno &dummySaIn);

        SvEvent(const BreakpointReduced &bp1In, const SuppAlignmentAnno &sa1In,
                GermlineMatch germlineInfo2, MrefMatch hitsInMref2In,
                const vector<pair<int, string>> &overhangDb,
                const SuppAlignmentAnno &dummySaIn);

        //	vector<int> getKey() const;
        string getKey() const;

        bool isGermline() const { return germline; }

        int getEventSize() const { return eventSize; }

        bool isInverted() const { return inverted; }

        int getTotalEvidence1() const { return totalEvidence1; }

        int getTotalEvidence2() const { return totalEvidence2; }

        int getEventScore() const { return eventScore; }

        int getSuspicious() const { return suspicious; }

        double getMateRatio1() const { return mateRatio1; }

        double getMateRatio2() const { return mateRatio2; }

        short getEvidenceLevel1() const { return evidenceLevel1; }

        short getEvidenceLevel2() const { return evidenceLevel2; }

        bool isSemiSuspicious() const { return semiSuspicious; }

        bool isDistant() const { return distant; }

        const SuppAlignmentAnno &getSelectedSa1() const { return selectedSa1; }

        const SuppAlignmentAnno &getSelectedSa2() const { return selectedSa2; }

        string printMatch(const vector<pair<int, string>> &overhangDb) const;

        bool isToRemove() const { return toRemove; }

        void setToRemove(bool toRemove) { this->toRemove = toRemove; }

        int getContaminationCandidate() const { return contaminationCandidate; }

        void setEventScore(int eventScore) { this->eventScore = eventScore; }

        void setEventType(unsigned int eventType) { this->eventType = eventType; }

        bool isOverhang1Compensation() const { return overhang1Compensation; }

        double getOverhang1lengthRatio() const { return overhang1lengthRatio; }
        double getOverhang2lengthRatio() const { return overhang2lengthRatio; }

      private:

        pair<int, double> mateQualityConditions(const SuppAlignmentAnno &sa);

        pair<bool, int> assessOverhangQualityCompensation(
            int lineIndex,
            const vector<pair<int, string>> &overhangDb) const;

        pair<bool, short> processMrefHits(const MrefMatch &hitsInMref,
                                          const SuppAlignmentAnno &sa,
                                          int evidenceLevelIn) const;

        double determineGermlineClonalityBp(const BreakpointReduced &bp1,
                                            const SuppAlignmentAnno &sa,
                                            double clonalityInit) const;

        void determineEventTypeAndSize(ChrPosition pos1,
                                       ChrPosition pos2,
                                       bool matchEncounteredM);

        int filterMatch(const BreakpointReduced &bp1,
                        const BreakpointReduced &bp2);

        int filterMatchSingle(const BreakpointReduced &bp1,
                              const BreakpointReduced &bp2);

        int filterMatchUnknown(const BreakpointReduced &bp1);

        pair<double, double> assessSvClonality(const BreakpointReduced &bp,
                                               int eventSupportTotal) const;

        ClonalityStatus
        assessBreakpointClonalityStatus(double clonalityRatioIn,
                                        const BreakpointReduced &bp1,
                                        const BreakpointReduced &bp2) const;

        ClonalityStatus
        assessBreakpointClonalityStatusSingle(double clonalityRatioIn,
                                              const BreakpointReduced &bp1,
                                              const BreakpointReduced &bp2) const;

        ClonalityStatus
        assessBreakpointClonalityStatusUnknown(double clonalityRatioIn,
                                               const BreakpointReduced &bp1) const;

        void assessSvArtifactStatus(const BreakpointReduced &bp1,
                                    const BreakpointReduced &bp2);

        void assessSvArtifactStatusUnknown();

        int assessEventScore(bool hardClipSuspiciousCall, int inputScoreCategory);

        void assessContamination(const vector<pair<int, string>> &overhangDb);

        pair<int, double>
        assessContaminationSingleBp(int overhangIndex,
                                    const vector<pair<int, string>> &overhangDb,
                                    const SuppAlignmentAnno &selectedSa);

        string collapseRange(const vector<string> &vec,
                             const string &delimiter) const {
            if (vec.empty()) {
                return "_";
            } else {
                return boost::join(vec, delimiter);
            }
        }

        bool toRemove;


        int contaminationCandidate;

        ChrIndex chrIndex1;

        ChrPosition pos1;

        ChrIndex chrIndex2;

        ChrPosition pos2;

        int lineIndex1;

        int lineIndex2;

        // Could probably be an enum.
        unsigned int eventType;

        int eventSize;

        bool inverted;

        bool doubleSupport;

        bool distant;

        bool overhang1Compensation;

        bool overhang2Compensation;

        // Can have value -1
        int overhang1Index;
        int overhang2Index;

        double overhang1lengthRatio;

        double overhang2lengthRatio;

        int inputScore;

        int eventScore;

        int totalEvidence1;

        int span1;

        int totalEvidence2;

        int span2;

        short evidenceLevel1;
        short evidenceLevel2;

        short mrefHits1;

        bool mrefHits1Conservative;
        short mrefHits2;

        bool mrefHits2Conservative;

        bool germline;

        double germlineClonality1;

        bool germlineStatus1;
        double germlineClonality2;

        bool germlineStatus2;

        SuppAlignmentAnno selectedSa1;

        SuppAlignmentAnno selectedSa2;

        double mateRatio1;

        double mateRatio2;

        int suspicious;

        bool semiSuspicious;

        double artifactRatio1;

        double clonalityRatio1;

        ClonalityStatus clonalityStatus1;

        double artifactRatio2;

        double clonalityRatio2;

        ClonalityStatus clonalityStatus2;

        ArtifactStatus artifactStatus;
    };

} /* namespace sophia */

#endif /* MATCHINFO_H_ */
