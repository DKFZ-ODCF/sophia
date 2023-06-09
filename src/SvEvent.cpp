/*
 * PairInfo.cpp
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

#include "ChrConverter.h"
#include "strtk.hpp"
#include <SvEvent.h>

namespace sophia {

using namespace std;

int SvEvent::GERMLINEOFFSETTHRESHOLD{};
double SvEvent::RELAXEDBPFREQTHRESHOLD{};
double SvEvent::BPFREQTHRESHOLD{};
double SvEvent::ARTIFACTFREQLOWTHRESHOLD{};
double SvEvent::ARTIFACTFREQHIGHTHRESHOLD{};
double SvEvent::CLONALITYLOWTHRESHOLD{};
double SvEvent::CLONALITYSTRICTLOWTHRESHOLD{};
double SvEvent::CLONALITYHIGHTHRESHOLD{};
int SvEvent::HALFDEFAULTREADLENGTH{};
int SvEvent::GERMLINEDBLIMIT{};
string SvEvent::PIDSINMREFSTR{};
boost::format SvEvent::doubleFormatter{"%.3f"};
bool SvEvent::ABRIDGEDOUTPUT{true};
bool SvEvent::NOCONTROLMODE{false};
bool SvEvent::DEBUGMODE{false};
const vector<string> SvEvent::EVENTTYPES{"UNKNOWN", "DEL", "DUP",
                                         "TRA",     "INV", "CONTAMINATION"};

SvEvent::SvEvent(const BreakpointReduced &bp1In, const BreakpointReduced &bp2In,
                 const SuppAlignmentAnno &sa1In, const SuppAlignmentAnno &sa2In,
                 const vector<pair<int, string>> &overhangDb)
    : toRemove{false}, contaminationCandidate{0},
      chrIndex1{bp1In.getChrIndex()}, pos1{bp1In.getPos()},
      chrIndex2{bp2In.getChrIndex()}, pos2{bp2In.getPos()},
      lineIndex1{bp1In.getLineIndex()},
      lineIndex2{bp2In.getLineIndex()}, eventType{0}, eventSize{0},
      inverted{sa1In.isInverted() || sa2In.isInverted()}, distant{false},
      overhang1Compensation{false}, overhang2Compensation{false},
      overhang1Index{-1}, overhang2Index{-1}, overhang1lengthRatio{0},
      overhang2lengthRatio{0}, inputScore{2}, eventScore{0},
      totalEvidence1{sa1In.getSupport() + sa1In.getSecondarySupport() +
                     sa1In.getMateSupport()},
      span1{bp1In.getNormalSpans()},
      totalEvidence2{sa2In.getSupport() + sa2In.getSecondarySupport() +
                     sa2In.getMateSupport()},
      span2{bp2In.getNormalSpans()}, evidenceLevel1{0},
      evidenceLevel2{0}, mrefHits1{bp1In.getMrefHits().getNumConsevativeHits()},
      mrefHits1Conservative{true},
      mrefHits2{bp2In.getMrefHits().getNumConsevativeHits()},
      mrefHits2Conservative{true}, germline{false},
      germlineClonality1{bp1In.getGermlineInfo().getConservativeClonality()},
      germlineStatus1{bp1In.getGermlineInfo().getConservativeClonality() >
                      0.15},
      germlineClonality2{bp2In.getGermlineInfo().getConservativeClonality()},
      germlineStatus2{bp2In.getGermlineInfo().getConservativeClonality() >
                      0.15},
      selectedSa1{sa1In}, selectedSa2{sa2In},
      mateRatio1{sa1In.getExpectedDiscordants() > 0
                     ? sa1In.getMateSupport() /
                           (0.0 + sa1In.getExpectedDiscordants())
                     : 1.0},
      mateRatio2{sa2In.getExpectedDiscordants() > 0
                     ? sa2In.getMateSupport() /
                           (0.0 + sa2In.getExpectedDiscordants())
                     : 1.0},
      suspicious{0}, semiSuspicious{sa1In.isSemiSuspicious() ||
                                    sa2In.isSemiSuspicious()} {
    //	auto messageMode = selectedSa1.getChrIndex() == 11 &&
    //(selectedSa1.getPos() == 2261373 && selectedSa2.getPos() == 2148480);
    auto posDifferential = pos1 - pos2;
    determineEventTypeAndSize(posDifferential, selectedSa2.isEncounteredM());
    if (chrIndex1 != chrIndex2) {
        distant = true;
    } else if (selectedSa1.getExpectedDiscordants() > 0 ||
               selectedSa2.getExpectedDiscordants() > 0) {
        distant = true;
    } else if (eventSize > 1500) {
        distant = true;
    }
    if (distant && chrIndex1 == chrIndex2 &&
        (selectedSa1.isFuzzy() || selectedSa1.isStrictFuzzyCandidate() ||
         selectedSa2.isFuzzy() || selectedSa2.isStrictFuzzyCandidate())) {
        if (eventSize < 5000) {
            distant = false;
        }
    }
    auto clonalityRes1 = assessSvClonality(
        bp1In, selectedSa1.getSupport() + selectedSa1.getSecondarySupport() +
                   selectedSa1.getMateSupport());
    artifactRatio1 = clonalityRes1.first;
    clonalityRatio1 = clonalityRes1.second;
    clonalityStatus1 =
        assessBreakpointClonalityStatus(clonalityRatio1, bp1In, bp2In);
    auto clonalityRes2 = assessSvClonality(
        bp2In, selectedSa2.getSupport() + selectedSa2.getSecondarySupport() +
                   selectedSa2.getMateSupport());
    artifactRatio2 = clonalityRes2.first;
    clonalityRatio2 = clonalityRes2.second;
    clonalityStatus2 =
        assessBreakpointClonalityStatus(clonalityRatio2, bp1In, bp2In);

    auto res1 = assessOverhangQualityCompensation(lineIndex1, overhangDb);
    overhang1Index = res1.second;
    overhang1Compensation = (clonalityStatus1 != EXTREME_SUBCLONAL) &&
                            selectedSa1.isDistant() && res1.first;
    auto res2 = assessOverhangQualityCompensation(lineIndex2, overhangDb);
    overhang2Index = res2.second;
    overhang2Compensation = (clonalityStatus2 != EXTREME_SUBCLONAL) &&
                            selectedSa2.isDistant() && res2.first;

    auto doubleSemiSuspicious =
        (selectedSa1.isSemiSuspicious() && selectedSa2.isSemiSuspicious());
    germlineClonality1 =
        determineGermlineClonalityBp(bp1In, selectedSa1, germlineClonality1);
    germlineStatus1 = germlineClonality1 > 0.15;
    germlineClonality2 =
        determineGermlineClonalityBp(bp2In, selectedSa2, germlineClonality2);
    germlineStatus2 = germlineClonality2 > 0.15;

    auto strictNonDecoy = !selectedSa1.isProperPairErrorProne() &&
                          !selectedSa2.isProperPairErrorProne() &&
                          ChrConverter::indexConverter[chrIndex1] < 23 &&
                          ChrConverter::indexConverter[chrIndex2] < 23;
    auto splitSupportThreshold1 =
        (strictNonDecoy && !selectedSa1.isSemiSuspicious() &&
         (mateRatio1 >= 0.6))
            ? 0
            : 2;
    auto splitSupportThreshold2 =
        (strictNonDecoy && !selectedSa2.isSemiSuspicious() &&
         (mateRatio2 >= 0.6))
            ? 0
            : 2;

    if (selectedSa1.getSupport() > splitSupportThreshold1) {
        ++evidenceLevel1;
    } else {
        if (strictNonDecoy && selectedSa1.getSupport() > 0 &&
            selectedSa1.getSecondarySupport() > splitSupportThreshold1) {
            ++evidenceLevel1;
        }
    }
    if (selectedSa1.getSecondarySupport() > splitSupportThreshold1) {
        ++evidenceLevel1;
    } else {
        if (strictNonDecoy && selectedSa1.getSecondarySupport() > 0 &&
            selectedSa1.getSupport() > splitSupportThreshold1) {
            ++evidenceLevel1;
        }
    }
    if (selectedSa1.isDistant()) {
        auto mateQuality = mateQualityConditions(selectedSa1);
        if (selectedSa1.getMateSupport() > mateQuality.first &&
            mateRatio1 >= mateQuality.second) {
            ++evidenceLevel1;
            if (evidenceLevel1 < 3 &&
                (overhang1Compensation ||
                 (!semiSuspicious && overhang2Compensation)) &&
                !doubleSemiSuspicious && strictNonDecoy) {
                if ((selectedSa1.getMateSupport() > 2) ||
                    (selectedSa1.getMateSupport() < 3 &&
                     selectedSa1.getExpectedDiscordants() ==
                         selectedSa1.getMateSupport())) {
                    ++evidenceLevel1;
                }
            }
        }
    }

    if (selectedSa2.getSupport() > splitSupportThreshold2) {
        ++evidenceLevel2;
    } else {
        if (strictNonDecoy && selectedSa2.getSupport() > 0 &&
            selectedSa2.getSecondarySupport() > splitSupportThreshold2) {
            ++evidenceLevel2;
        }
    }
    if (selectedSa2.getSecondarySupport() > splitSupportThreshold2) {
        ++evidenceLevel2;
    } else {
        if (strictNonDecoy && selectedSa2.getSecondarySupport() > 0 &&
            selectedSa2.getSupport() > splitSupportThreshold2) {
            ++evidenceLevel2;
        }
    }
    if (selectedSa2.isDistant()) {
        auto mateQuality = mateQualityConditions(selectedSa2);
        if (selectedSa2.getMateSupport() > mateQuality.first &&
            mateRatio2 >= mateQuality.second) {
            ++evidenceLevel2;
            if (evidenceLevel2 < 3 &&
                ((!semiSuspicious && overhang1Compensation) ||
                 overhang2Compensation) &&
                !doubleSemiSuspicious && strictNonDecoy) {
                if ((selectedSa2.getMateSupport() > 2) ||
                    (selectedSa2.getMateSupport() < 3 &&
                     selectedSa2.getExpectedDiscordants() ==
                         selectedSa2.getMateSupport())) {
                    ++evidenceLevel2;
                }
            }
        }
    }
    auto mrefHits1Tmp =
        processMrefHits(bp1In.getMrefHits(), selectedSa1, evidenceLevel1);
    mrefHits1 = mrefHits1Tmp.second;
    mrefHits1Conservative = mrefHits1Tmp.first;
    auto mrefHits2Tmp =
        processMrefHits(bp2In.getMrefHits(), selectedSa2, evidenceLevel2);
    mrefHits2 = mrefHits2Tmp.second;
    mrefHits2Conservative = mrefHits2Tmp.first;
    if (!germlineStatus1 && germlineClonality1 > 0 &&
        mrefHits1 > GERMLINEDBLIMIT) {
        germlineStatus1 = true;
    }
    if (!germlineStatus2 && germlineClonality2 > 0 &&
        mrefHits2 > GERMLINEDBLIMIT) {
        germlineStatus2 = true;
    }
    germline =
        (germlineStatus1 || germlineStatus2) &&
        !((selectedSa1.getSupport() + selectedSa1.getSecondarySupport()) >
              200 &&
          (selectedSa1.getSupport() + selectedSa2.getSecondarySupport()) > 200);

    assessSvArtifactStatus(bp1In, bp2In);
    if (!selectedSa1.isSemiSuspicious() && selectedSa2.isSemiSuspicious()) {
        //		auto messageMode = selectedSa1.getChrIndex() == 21 &&
        //selectedSa1.getPos() == 48119076;
        auto score1 = filterMatch(bp1In, bp2In);
        auto score2 = filterMatchUnknown(bp1In);
        //		if (messageMode) cerr << score1 << " " << score2 <<
        //"\n";
        if (score2 == 0 && score1 != 0) {
            suspicious = score2;
            semiSuspicious = false;
        } else {
            suspicious = score1;
        }
    } else {
        suspicious = filterMatch(bp1In, bp2In);
    }
    eventScore = assessEventScore(false, inputScore);
    if (suspicious == 0 && eventScore > 2) {
        assessContamination(overhangDb);
    }
}

SvEvent::SvEvent(const BreakpointReduced &bp1In, const BreakpointReduced &bp2In,
                 const SuppAlignmentAnno &sa1In,
                 const vector<pair<int, string>> &overhangDb,
                 const SuppAlignmentAnno &dummySaIn)
    : toRemove{false}, contaminationCandidate{0},
      chrIndex1{bp1In.getChrIndex()}, pos1{bp1In.getPos()},
      chrIndex2{bp2In.getChrIndex()}, pos2{bp2In.getPos()},
      lineIndex1{bp1In.getLineIndex()}, lineIndex2{bp2In.getLineIndex()},
      eventType{0}, eventSize{0}, inverted{sa1In.isInverted()}, distant{false},
      overhang1Compensation{false}, overhang2Compensation{false},
      overhang1Index{-1}, overhang2Index{-1}, overhang1lengthRatio{0},
      overhang2lengthRatio{0}, inputScore{1}, eventScore{0},
      totalEvidence1{sa1In.getSupport() + sa1In.getSecondarySupport() +
                     sa1In.getMateSupport()},
      span1{bp1In.getNormalSpans()},
      totalEvidence2{bp2In.getPairedBreaksSoft() + bp2In.getPairedBreaksHard() +
                     bp2In.getUnpairedBreaksSoft() +
                     bp2In.getUnpairedBreaksHard() +
                     bp2In.getBreaksShortIndel() + bp2In.getMateSupport()},
      span2{bp2In.getNormalSpans()}, evidenceLevel1{0},
      evidenceLevel2{0}, mrefHits1{bp1In.getMrefHits().getNumConsevativeHits()},
      mrefHits1Conservative{true},
      mrefHits2{bp2In.getMrefHits().getNumConsevativeHits()},
      mrefHits2Conservative{true}, germline{false},
      germlineClonality1{bp1In.getGermlineInfo().getConservativeClonality()},
      germlineStatus1{bp1In.getGermlineInfo().getConservativeClonality() >
                      0.15},
      germlineClonality2{bp2In.getGermlineInfo().getConservativeClonality()},
      germlineStatus2{bp2In.getGermlineInfo().getConservativeClonality() >
                      0.15},
      selectedSa1{sa1In}, selectedSa2{dummySaIn},
      mateRatio1{sa1In.getExpectedDiscordants() > 0
                     ? sa1In.getMateSupport() /
                           (0.0 + sa1In.getExpectedDiscordants())
                     : 1.0},
      mateRatio2{1.0}, suspicious{0}, semiSuspicious{sa1In.isSemiSuspicious()} {

    auto posDifferential = pos1 - pos2;
    determineEventTypeAndSize(
        posDifferential, (bp2In.getLeftCoverage() > bp2In.getRightCoverage()));
    if (chrIndex1 != chrIndex2) {
        distant = true;
    } else if (selectedSa1.getExpectedDiscordants() > 0) {
        distant = true;
    } else if (eventSize > 1500) {
        distant = true;
    }
    if (distant && chrIndex1 == chrIndex2 &&
        (selectedSa1.isFuzzy() || selectedSa1.isStrictFuzzyCandidate())) {
        if (eventSize < 5000) {
            distant = false;
        }
    }
    if (distant && (eventSize > 0) && (eventSize < 1000)) {
        distant = false;
    }

    auto res1 = assessOverhangQualityCompensation(lineIndex1, overhangDb);
    overhang1Index = res1.second;
    auto res2 = assessOverhangQualityCompensation(lineIndex2, overhangDb);
    overhang2Index = res2.second;

    auto mateEvidence1 = false;
    auto mateEvidence2 = false;
    germlineClonality1 =
        determineGermlineClonalityBp(bp1In, selectedSa1, germlineClonality1);
    germlineStatus1 = germlineClonality1 > 0.15;

    auto strictNonDecoy = !selectedSa1.isProperPairErrorProne() &&
                          ChrConverter::indexConverter[chrIndex1] < 23 &&
                          ChrConverter::indexConverter[chrIndex2] < 23;
    auto splitSupportThreshold =
        (strictNonDecoy && !selectedSa2.isSemiSuspicious() &&
         (mateRatio1 >= 0.66))
            ? 0
            : 2;

    if (selectedSa1.getSupport() > splitSupportThreshold) {
        ++evidenceLevel1;
    } else {
        if (strictNonDecoy && selectedSa1.getSupport() > 0 &&
            selectedSa1.getSecondarySupport() > splitSupportThreshold) {
            ++evidenceLevel1;
        }
    }
    if (selectedSa1.getSecondarySupport() > splitSupportThreshold) {
        ++evidenceLevel1;
    } else {
        if (strictNonDecoy && selectedSa1.getSecondarySupport() > 0 &&
            selectedSa1.getSupport() > splitSupportThreshold) {
            ++evidenceLevel1;
        }
    }
    if (selectedSa1.isDistant()) {
        if (!selectedSa1.isStrictFuzzyCandidate() ||
            (selectedSa1.isStrictFuzzyCandidate() &&
             selectedSa1.getMateSupport() > 4)) {
            auto mateQualityCriteria = selectedSa1.isProperPairErrorProne() ||
                                       selectedSa1.isSemiSuspicious() ||
                                       selectedSa1.isStrictFuzzyCandidate();
            if (mateRatio1 >= 0.4 &&
                (!mateQualityCriteria ||
                 (mateQualityCriteria && selectedSa1.getMateSupport() > 4))) {
                ++evidenceLevel1;
                mateEvidence1 = true;
            }
        }
    }

    if (bp2In.getPairedBreaksSoft() + bp2In.getUnpairedBreaksSoft() > 0) {
        ++evidenceLevel2;
    }
    if (bp2In.getPairedBreaksHard() + bp2In.getUnpairedBreaksHard() > 0) {
        ++evidenceLevel2;
    }
    if (selectedSa1.isDistant()) {
        if (bp2In.getMateSupport() > 4) {
            ++evidenceLevel2;
            mateEvidence2 = true;
        }
    }
    auto mrefHits1Tmp =
        processMrefHits(bp1In.getMrefHits(), selectedSa1, evidenceLevel1);
    mrefHits1 = mrefHits1Tmp.second;
    mrefHits1Conservative = mrefHits1Tmp.first;
    auto mrefHits2Tmp =
        processMrefHits(bp2In.getMrefHits(), selectedSa2, evidenceLevel2);
    mrefHits2 = mrefHits2Tmp.second;
    mrefHits2Conservative = mrefHits2Tmp.first;
    if (!germlineStatus1 && germlineClonality1 > 0 &&
        mrefHits1 > GERMLINEDBLIMIT) {
        germlineStatus1 = true;
    }
    if (!germlineStatus2 && germlineClonality2 > 0 &&
        mrefHits2 > GERMLINEDBLIMIT) {
        germlineStatus2 = true;
    }
    germline = (germlineStatus1 || germlineStatus2) &&   //
               !((selectedSa1.getSupport() +
                  selectedSa1.getSecondarySupport()) > 200 &&
                 (bp2In.getPairedBreaksSoft() + bp2In.getPairedBreaksHard() +
                  bp2In.getUnpairedBreaksSoft() +
                  bp2In.getUnpairedBreaksHard()) > 200);

    if (!distant) {
        if (mateEvidence1) {
            --evidenceLevel1;
            totalEvidence1 -= selectedSa1.getMateSupport();
        }
        if (mateEvidence2) {
            --evidenceLevel2;
            totalEvidence2 -= bp2In.getMateSupport();
        }
    }
    auto clonalityRes1 = assessSvClonality(
        bp1In, selectedSa1.getSupport() + selectedSa1.getSecondarySupport() +
                   selectedSa1.getMateSupport());
    artifactRatio1 = clonalityRes1.first;
    clonalityRatio1 = clonalityRes1.second;
    clonalityStatus1 =
        assessBreakpointClonalityStatusSingle(clonalityRatio1, bp1In, bp2In);
    auto clonalityRes2 = assessSvClonality(
        bp2In, bp2In.getPairedBreaksSoft() + bp2In.getPairedBreaksHard() +
                   bp2In.getUnpairedBreaksSoft());
    artifactRatio2 = clonalityRes2.first;
    clonalityRatio2 = clonalityRes2.second;
    clonalityStatus2 =
        assessBreakpointClonalityStatusSingle(clonalityRatio2, bp1In, bp2In);

    assessSvArtifactStatus(bp1In, bp2In);
    suspicious = filterMatchSingle(bp1In, bp2In);
    auto hardClipSuspiciousCall =
        selectedSa1.getSupport() == 0 &&
        selectedSa1.getSecondarySupport() > 0 &&
        (selectedSa1.getMateSupport() <= selectedSa1.getSecondarySupport() + 4);
    eventScore = assessEventScore(hardClipSuspiciousCall, inputScore);
    if (suspicious == 0 && eventScore > 2) {
        assessContamination(overhangDb);
    }
}

SvEvent::SvEvent(const BreakpointReduced &bp1In, const SuppAlignmentAnno &sa1In,
                 GermlineMatch germlineInfo2, MrefMatch hitsInMref2In,
                 const vector<pair<int, string>> &overhangDb,
                 const SuppAlignmentAnno &dummySaIn)
    : toRemove{false}, contaminationCandidate{0},
      chrIndex1{bp1In.getChrIndex()}, pos1{bp1In.getPos()},
      chrIndex2{sa1In.getChrIndex()}, pos2{sa1In.getPos()},
      lineIndex1{bp1In.getLineIndex()}, lineIndex2{-1}, eventType{0},
      eventSize{0}, inverted{sa1In.isInverted()}, distant{false},
      overhang1Compensation{false}, overhang2Compensation{false},
      overhang1Index{-1}, overhang2Index{-1}, overhang1lengthRatio{0},
      overhang2lengthRatio{0}, inputScore{0}, eventScore{0},
      totalEvidence1{sa1In.getSupport() + sa1In.getSecondarySupport() +
                     sa1In.getMateSupport()},
      span1{bp1In.getNormalSpans()}, totalEvidence2{0}, evidenceLevel1{0},
      evidenceLevel2{0}, mrefHits1{bp1In.getMrefHits().getNumConsevativeHits()},
      mrefHits1Conservative{true},
      mrefHits2{hitsInMref2In.getNumConsevativeHits()},
      mrefHits2Conservative{true}, germline{false},
      germlineClonality1{bp1In.getGermlineInfo().getConservativeClonality()},
      germlineStatus1{bp1In.getGermlineInfo().getConservativeClonality() >
                      0.15},
      germlineClonality2{germlineInfo2.getClonality()},
      germlineStatus2{germlineInfo2.getClonality() > 0.15}, selectedSa1{sa1In},
      selectedSa2{dummySaIn},
      mateRatio1{sa1In.getExpectedDiscordants() > 0
                     ? sa1In.getMateSupport() /
                           (0.0 + sa1In.getExpectedDiscordants())
                     : 1.0},
      mateRatio2{1.0}, suspicious{0}, semiSuspicious{sa1In.isSemiSuspicious()} {
    auto truePos2 = pos2;
    if (chrIndex1 == chrIndex2) {
        if (abs(pos1 - pos2) > abs(pos1 - sa1In.getExtendedPos())) {
            truePos2 = sa1In.getExtendedPos();
        }
    }
    auto posDifferential = pos1 - truePos2;
    determineEventTypeAndSize(posDifferential, !selectedSa1.isEncounteredM());
    if (chrIndex1 != chrIndex2) {
        distant = true;
    } else if (selectedSa1.getExpectedDiscordants() > 0) {
        distant = true;
    } else if (eventSize > 2500) {
        distant = true;
    }
    if (distant && selectedSa1.isFuzzy() &&
        bp1In.getChrIndex() == selectedSa1.getChrIndex() &&
        selectedSa1.isStrictFuzzyCandidate()) {
        auto fuzDiff = selectedSa1.getExtendedPos() - selectedSa1.getPos();
        if (max(0, pos1 - fuzDiff) <= selectedSa1.getExtendedPos() &&
            selectedSa1.getPos() <= (pos1 + fuzDiff)) {
            distant = false;
        } else if (eventSize < 5000) {
            distant = false;
        }
    }
    if (distant && (eventSize > 0) && (eventSize < 2000)) {
        distant = false;
    }
    if (chrIndex1 == chrIndex2 && !distant && eventType == 3) {
        suspicious = 3;
        return;
    }
    auto clonalityRes1 = assessSvClonality(
        bp1In, selectedSa1.getSupport() + selectedSa1.getSecondarySupport() +
                   selectedSa1.getMateSupport());
    artifactRatio1 = clonalityRes1.first;
    clonalityRatio1 = clonalityRes1.second;
    clonalityStatus1 =
        assessBreakpointClonalityStatusUnknown(clonalityRatio1, bp1In);

    auto res1 = assessOverhangQualityCompensation(lineIndex1, overhangDb);
    overhang1Index = res1.second;
    overhang1Compensation = (clonalityStatus1 != EXTREME_SUBCLONAL) &&
                            selectedSa1.isDistant() && res1.first;

    auto mateEvidence1 = false;
    auto additionalEvidence1 = false;
    germlineClonality1 =
        determineGermlineClonalityBp(bp1In, selectedSa1, germlineClonality1);
    germlineStatus1 = germlineClonality1 > 0.15;

    auto strictNonDecoy = !selectedSa1.isProperPairErrorProne() &&
                          ChrConverter::indexConverter[chrIndex1] < 23 &&
                          ChrConverter::indexConverter[chrIndex2] < 23;
    auto splitSupportThreshold =
        (strictNonDecoy && (mateRatio1 >= 0.66) ? 0 : 2);

    if (selectedSa1.getSupport() > splitSupportThreshold) {
        ++evidenceLevel1;
    } else {
        if (strictNonDecoy && selectedSa1.getSupport() > 0 &&
            selectedSa1.getSecondarySupport() > splitSupportThreshold) {
            ++evidenceLevel1;
        }
    }
    if (selectedSa1.getSecondarySupport() > splitSupportThreshold) {
        ++evidenceLevel1;
    } else {
        if (strictNonDecoy && selectedSa1.getSecondarySupport() > 0 &&
            selectedSa1.getSupport() > splitSupportThreshold) {
            ++evidenceLevel1;
        }
    }
    if (selectedSa1.isDistant()) {
        if (mateRatio1 >= 0.4) {
            if (!(selectedSa1.isStrictFuzzyCandidate() ||
                  selectedSa1.isProperPairErrorProne()) ||
                ((selectedSa1.isStrictFuzzyCandidate() ||
                  selectedSa1.isProperPairErrorProne()) &&
                 selectedSa1.getMateSupport() > 4)) {
                ++evidenceLevel1;
                mateEvidence1 = true;
                if (evidenceLevel1 < 3 && overhang1Compensation &&
                    strictNonDecoy) {
                    if ((selectedSa1.getMateSupport() > 2) ||
                        (selectedSa1.getMateSupport() < 3 &&
                         selectedSa1.getExpectedDiscordants() ==
                             selectedSa1.getMateSupport())) {
                        ++evidenceLevel1;
                        additionalEvidence1 = true;
                    }
                }
            }
        }
    }
    auto mrefHits1Tmp =
        processMrefHits(bp1In.getMrefHits(), selectedSa1, evidenceLevel1);
    mrefHits1 = mrefHits1Tmp.second;
    mrefHits1Conservative = mrefHits1Tmp.first;
    auto mrefHits2Tmp = processMrefHits(hitsInMref2In, selectedSa2, 0);
    mrefHits2 = mrefHits2Tmp.second;
    mrefHits2Conservative = mrefHits2Tmp.first;
    if (!germlineStatus1 && germlineClonality1 > 0 &&
        mrefHits1 > GERMLINEDBLIMIT) {
        germlineStatus1 = true;
    }
    if (!germlineStatus2 && germlineClonality2 > 0 &&
        mrefHits2 > GERMLINEDBLIMIT) {
        germlineStatus2 = true;
    }

    germline =
        (germlineStatus1 || germlineStatus2) &&
        !((selectedSa1.getSupport() + selectedSa1.getSecondarySupport()) > 200);
    if (!distant) {
        if (mateEvidence1) {
            --evidenceLevel1;
            totalEvidence1 -= selectedSa1.getMateSupport();
        }
        if (additionalEvidence1) {
            --evidenceLevel1;
        }
    }
    assessSvArtifactStatusUnknown();
    suspicious = filterMatchUnknown(bp1In);
    auto hardClipSuspiciousCall =
        selectedSa1.getSupport() == 0 &&
        selectedSa1.getSecondarySupport() > 0 &&
        (selectedSa1.getMateSupport() <= selectedSa1.getSecondarySupport() + 4);
    eventScore = assessEventScore(hardClipSuspiciousCall, inputScore);
    if (suspicious == 0 && eventScore > 2) {
        assessContamination(overhangDb);
    }
}

void
SvEvent::determineEventTypeAndSize(int posDifferential,
                                   bool matchEncounteredM) {
    //	static vector<string> EVENTTYPES { "UNKNOWN", "DEL", "DUP", "TRA",
    //"INV", "CONTAMINATION" };
    if (chrIndex1 != chrIndex2) {
        eventType = 3;
        eventSize = -1;
    } else {
        eventSize = abs(posDifferential);
        if (posDifferential < 0) {
            if (inverted) {
                eventType = 4;
            } else {
                if (selectedSa1.isEncounteredM() && !matchEncounteredM) {
                    eventType = 1;
                } else if (!selectedSa1.isEncounteredM() && matchEncounteredM) {
                    eventType = 2;
                } else {
                    eventType = 3;
                }
            }
        } else if (posDifferential > 0) {
            if (inverted) {
                eventType = 4;
            } else {
                if (selectedSa1.isEncounteredM() && !matchEncounteredM) {
                    eventType = 2;
                } else if (!selectedSa1.isEncounteredM() && matchEncounteredM) {
                    eventType = 1;
                } else {
                    eventType = 3;
                }
            }
        } else {
            eventType = 3;
            suspicious = 3;
        }
    }
}

pair<int, double>
SvEvent::mateQualityConditions(const SuppAlignmentAnno &sa) {
    //	auto messageMode = sa.getChrIndex() == 11 && (sa.getPos() == 2261373 ||
    //sa.getPos() == 2148480);
    auto doubleSemiSuspicious =
        (selectedSa1.isSemiSuspicious() && selectedSa2.isSemiSuspicious()) ||
        (selectedSa1.isSemiSuspicious() &&
         selectedSa2.isProperPairErrorProne()) ||
        (selectedSa1.isProperPairErrorProne() &&
         selectedSa2.isSemiSuspicious());
    auto mateLowQualityCriteriaTier1 = sa.isProperPairErrorProne() ||
                                       doubleSemiSuspicious ||
                                       (mateRatio1 + mateRatio2) < 1.1;
    auto mateLowQualityCriteriaTier2 =
        sa.isSemiSuspicious() && sa.isStrictFuzzyCandidate();
    auto mateLowQualityCriteriaTier3 = sa.isSemiSuspicious() || sa.isFuzzy();
    auto mateLowQualityCriteriaTier4 = sa.isStrictFuzzyCandidate();
    if (mateLowQualityCriteriaTier1 && mateLowQualityCriteriaTier2) {
        return {9, 0.8};
    } else if (mateLowQualityCriteriaTier1 || inputScore == 0) {
        return {4, 0.8};
    } else if (mateLowQualityCriteriaTier2) {
        return {4, 0.6};
    } else if (mateLowQualityCriteriaTier3) {
        return {2, 0.6};
    } else if (mateLowQualityCriteriaTier4) {
        return {0, 0.6};
    } else if (sa.getMateSupport() < 10) {
        return {0, 0.4};
    } else {
        return {0, 0.33};
    }
}

pair<bool, int>
SvEvent::assessOverhangQualityCompensation(
    int lineIndex, const vector<pair<int, string>> &overhangDb) const {
    auto overhangIndex = -1;
    auto compensation = false;
    pair<int, string> dummy{lineIndex, ""};
    auto lower = lower_bound(overhangDb.cbegin(), overhangDb.cend(), dummy);
    auto it = overhangDb.cend();
    if (lower != overhangDb.cend()) {
        if (lower->first == lineIndex) {
            it = lower;
        } else if (next(lower) != overhangDb.cend() &&
                   next(lower)->first == lineIndex) {
            it = next(lower);
        } else if (lower != overhangDb.cbegin() &&
                   prev(lower)->first == lineIndex) {
            it = prev(lower);
        }
    }
    if (it != overhangDb.cend()) {
        auto overhangLength = 0;
        auto maxOverhangLength = 0;
        auto counts = 0;
        overhangIndex = it - overhangDb.cbegin();
        for (const auto c : it->second) {
            switch (c) {
            case '(':
                maxOverhangLength = max(maxOverhangLength, overhangLength);
                overhangLength = 0;
                break;
            case ':':
                overhangLength = 0;
                ++counts;
                break;
            default:
                ++overhangLength;
                break;
            }
        }
        if (counts < 3) {
            auto lengthRatio = (0.0 + maxOverhangLength) /
                               SuppAlignmentAnno::DEFAULTREADLENGTH;
            compensation = 0.25 <= lengthRatio && lengthRatio <= 0.8;
        }
    }
    return {compensation, overhangIndex};
}

pair<bool, short>
SvEvent::processMrefHits(const MrefMatch &hitsInMref,
                         const SuppAlignmentAnno &sa,
                         int evidenceLevelIn) const {
    auto initScore = hitsInMref.getNumConsevativeHits();
    auto distantHit = false;
    auto saMatch = false;
    short maxScore{0};
    short hits{initScore};
    for (const auto &saRef : hitsInMref.getSuppMatches()) {
        if (saRef.saCloseness(sa, SuppAlignmentAnno::DEFAULTREADLENGTH * 6)) {
            saMatch = true;
            auto score = saRef.getSecondarySupport();
            if (score > maxScore) {
                maxScore = score;
                if (maxScore > initScore) {
                    hits = maxScore;
                    distantHit = true;
                }
            }
        }
    }
    if (!distantHit &&
        (selectedSa1.isStrictFuzzy() || selectedSa2.isStrictFuzzy()) &&
        evidenceLevelIn < 3) {
        hits = hitsInMref.getNumHits();
    }
    return {!saMatch, hits};
}

double
SvEvent::determineGermlineClonalityBp(const BreakpointReduced &bp1,
                                      const SuppAlignmentAnno &sa,
                                      double clonalityInit) const {
    auto maxClonality = clonalityInit;
    if (sa.isDistant()) {
        auto i = 0;
        for (const auto &saRef : bp1.getGermlineInfo().getSuppMatches()) {
            if (saRef.saCloseness(sa,
                                  SuppAlignmentAnno::DEFAULTREADLENGTH * 6)) {
                auto clonality = bp1.getGermlineInfo().getClonalities()[i];
                if (clonality > maxClonality) {
                    maxClonality = clonality;
                }
                break;
            }
            ++i;
        }
    }
    return maxClonality;
}

int
SvEvent::filterMatch(const BreakpointReduced &bp1,
                     const BreakpointReduced &bp2) {
    if (suspicious != 0) {
        return suspicious;
    }
    if (selectedSa1.isSuspicious() || selectedSa2.isSuspicious() ||
        (selectedSa1.isSemiSuspicious() && selectedSa2.isSemiSuspicious() &&
         !(evidenceLevel1 == 3 && evidenceLevel2 == 3))) {
        return 1;
    }
    if ((selectedSa1.getExpectedDiscordants() > 0 && mateRatio1 < 0.1) ||
        (selectedSa2.getExpectedDiscordants() > 0 && mateRatio2 < 0.1)) {
        return 2;
    }
    if (mrefHits1 > GERMLINEDBLIMIT) {
        if (selectedSa1.isSemiSuspicious() && evidenceLevel1 < 3) {
            if (!(selectedSa1.getSupport() > 9 ||
                  selectedSa1.getSecondarySupport() > 9 ||
                  selectedSa1.getMateSupport() > 9 || overhang1Compensation)) {
                return 4;
            }
        }
    }
    if (mrefHits2 > GERMLINEDBLIMIT) {
        if (selectedSa2.isSemiSuspicious() && evidenceLevel2 < 3) {
            if (!(selectedSa2.getSupport() > 9 ||
                  selectedSa2.getSecondarySupport() > 9 ||
                  selectedSa1.getMateSupport() > 9 || overhang2Compensation)) {
                return 5;
            }
        }
    }
    if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
        if (clonalityRatio1 < CLONALITYSTRICTLOWTHRESHOLD &&
            clonalityRatio2 < CLONALITYSTRICTLOWTHRESHOLD) {
            return 4;
        } else if ((!germlineStatus1 || !germlineStatus2) &&   //
                   (mrefHits1 > GERMLINEDBLIMIT ||
                    mrefHits2 > GERMLINEDBLIMIT) &&
                   (eventSize > 0) && (eventSize < HALFDEFAULTREADLENGTH)) {
            return 5;
        }
    }
    if (!distant) {
        if (semiSuspicious) {
            return 13;
        }
        if (selectedSa1.isStrictFuzzy() || selectedSa2.isStrictFuzzy()) {
            return 13;
        }
        if (mrefHits1 > GERMLINEDBLIMIT || mrefHits2 > GERMLINEDBLIMIT) {
            if (totalEvidence1 < 5 || totalEvidence2 < 5) {
                return 13;
            }
        }
    } else {
        if (bp1.getChrIndex() != bp2.getChrIndex()) {
            if (selectedSa1.getMateSupport() == 0) {
                return 6;
            } else if (selectedSa2.getMateSupport() == 0) {
                return 7;
            }
        }
        auto threshold =
            (semiSuspicious || selectedSa1.isProperPairErrorProne() ||
             selectedSa2.isProperPairErrorProne())
                ? 5
                : 3;
        if (selectedSa1.getSupport() < threshold &&
            selectedSa1.getSecondarySupport() < threshold &&
            selectedSa2.getMateSupport() < threshold &&
            selectedSa2.getSupport() < threshold &&
            selectedSa2.getSecondarySupport() < threshold &&
            selectedSa2.getMateSupport() < threshold) {
            return 8;
        }
        if (mateRatio1 < 0.4 && mateRatio2 < 0.4) {
            return 11;
        }
        if (mateRatio1 < 0.25 || mateRatio2 < 0.25) {
            return 12;
        }
        if (selectedSa1.isFuzzy() || selectedSa2.isFuzzy()) {
            if (bp1.getChrIndex() == bp2.getChrIndex()) {
                if (selectedSa1.isStrictFuzzy() ||
                    selectedSa2.isStrictFuzzy()) {
                    if (abs(bp1.getPos() - bp2.getPos()) < 5000) {
                        return 14;
                    }
                }
            }
            if (mrefHits1 > GERMLINEDBLIMIT || mrefHits2 > GERMLINEDBLIMIT) {
                return 173;
            }
            if (selectedSa1.isFuzzy() && selectedSa2.isFuzzy()) {
                if (mateRatio1 < 0.5 || mateRatio2 < 0.5) {
                    return 174;
                }
            }
            if (selectedSa1.isFuzzy() && selectedSa1.getSupport() == 0 &&
                selectedSa1.getSecondarySupport() == 0 &&
                selectedSa2.isSemiSuspicious()) {
                return 175;
            }
            if (selectedSa2.isFuzzy() && selectedSa2.getSupport() == 0 &&
                selectedSa2.getSecondarySupport() == 0 &&
                selectedSa1.isSemiSuspicious()) {
                return 176;
            }
        }
    }
    if (bp1.getChrIndex() != bp2.getChrIndex() ||
        abs(bp1.getPos() - bp2.getPos()) > 150) {
        auto eventTotal1 = bp1.getPairedBreaksSoft() +
                           bp1.getPairedBreaksHard() +
                           bp1.getUnpairedBreaksSoft();
        auto eventTotal2 = bp2.getPairedBreaksSoft() +
                           bp2.getPairedBreaksHard() +
                           bp2.getUnpairedBreaksSoft();
        if (bp1.getBreaksShortIndel() > eventTotal1 ||
            bp2.getBreaksShortIndel() > eventTotal2) {
            return 16;
        }
    }
    if (clonalityStatus1 == EXTREME_SUBCLONAL) {
        return 39;
    }
    if (clonalityStatus2 == EXTREME_SUBCLONAL) {
        return 40;
    }
    if (artifactStatus == ARTIFACT) {
        return 41;
    }
    if (NOCONTROLMODE || germlineStatus1 || germlineStatus2) {
        if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
            return 42;
        }
    }
    if (mrefHits1 > BPFREQTHRESHOLD) {
        if (NOCONTROLMODE || germlineStatus1 || germlineStatus2 ||
            mrefHits2 > GERMLINEDBLIMIT) {
            return 431;
        }
        if (mrefHits1 > RELAXEDBPFREQTHRESHOLD * 2.5) {
            return 431;
        }
        if (!(mrefHits1Conservative && evidenceLevel1 == 3 &&
              evidenceLevel2 > 1 && mrefHits2 == 0 && !selectedSa1.isFuzzy() &&
              !selectedSa2.isFuzzy())) {
            if (mrefHits1 > RELAXEDBPFREQTHRESHOLD) {
                return 431;
            }
        }
    }
    if (mrefHits2 > BPFREQTHRESHOLD) {
        if (NOCONTROLMODE || germlineStatus1 || germlineStatus2 ||
            mrefHits1 > GERMLINEDBLIMIT) {
            return 432;
        }
        if (mrefHits2 > RELAXEDBPFREQTHRESHOLD * 2.5) {
            return 432;
        }
        if (!(mrefHits2Conservative && evidenceLevel1 > 1 &&
              evidenceLevel2 == 3 && mrefHits1 == 0 && !selectedSa1.isFuzzy() &&
              !selectedSa2.isFuzzy())) {
            if (mrefHits2 > RELAXEDBPFREQTHRESHOLD) {
                return 432;
            }
        }
    }
    return suspicious;
}

int
SvEvent::filterMatchSingle(const BreakpointReduced &bp1,
                           const BreakpointReduced &bp2) {
    if (suspicious != 0) {
        return suspicious;
    }
    if (selectedSa1.isSuspicious() || semiSuspicious) {
        return 1;
    }
    if (selectedSa1.getExpectedDiscordants() > 0 && mateRatio1 < 0.1) {
        return 2;
    }
    if (mrefHits1 > GERMLINEDBLIMIT &&
        !(bp1.getChrIndex() == 999 || bp1.getChrIndex() == 1000)) {
        return 171;
    }
    if (mrefHits2 > GERMLINEDBLIMIT &&
        !(bp2.getChrIndex() == 999 || bp2.getChrIndex() == 1000)) {
        return 171;
    }
    if (!distant) {
        if (inverted) {
            return 22;
        }
        if (eventSize > 0 && eventSize < HALFDEFAULTREADLENGTH) {
            return 23;
        }
        if (eventType == 3) {
            return 23;
        }
        if (eventSize > 150) {
            auto eventTotal1 = bp1.getPairedBreaksSoft() +
                               bp1.getPairedBreaksHard() +
                               bp1.getUnpairedBreaksSoft();
            auto eventTotal2 = bp2.getPairedBreaksSoft() +
                               bp2.getPairedBreaksHard() +
                               bp2.getUnpairedBreaksSoft();
            if ((bp1.getBreaksShortIndel() / (0.0 + eventTotal1) > 0.5)) {
                return 20;
            } else if ((bp2.getBreaksShortIndel() / (0.0 + eventTotal2)) >
                       0.5) {
                return 21;
            }
        }
        if (selectedSa1.isFuzzy() || selectedSa1.isSemiSuspicious() ||
            evidenceLevel1 == 1 ||
            (selectedSa1.getSupport() + selectedSa1.getSecondarySupport()) <
                3) {
            return 25;
        }
        if (mrefHits1 > GERMLINEDBLIMIT || mrefHits2 > GERMLINEDBLIMIT) {
            if (totalEvidence1 < 5) {
                return 25;
            }
        }
    } else {
        if (selectedSa1.getSupport() < 3 &&
            selectedSa1.getSecondarySupport() < 3 &&
            selectedSa1.getMateSupport() < 3) {
            return 18;
        }
        if (bp1.getChrIndex() != bp2.getChrIndex()) {
            if (selectedSa1.getMateSupport() == 0) {
                return 18;
            } else if (evidenceLevel1 == 1 && bp2.getMateSupport() == 0) {
                return 19;
            }
        }
        if (selectedSa1.getMateSupport() < 3 || mateRatio1 < 0.4) {
            return 24;
        }
        if (selectedSa1.isFuzzy() || (selectedSa1.getSupport() +
                                      selectedSa1.getSecondarySupport()) < 3) {
            if (bp1.getChrIndex() == bp2.getChrIndex()) {
                if (abs(bp1.getPos() - bp2.getPos()) < 5000) {
                    return 26;
                }
                if (inverted && abs(bp1.getPos() - bp2.getPos()) < 10000) {
                    return 27;
                }
            }
            if (selectedSa1.isFuzzy() && selectedSa1.getSupport() == 0 &&
                selectedSa1.getSecondarySupport() == 0) {
                return 271;
            }
        }
    }
    auto eventTotal1 = bp1.getPairedBreaksSoft() + bp1.getPairedBreaksHard() +
                       bp1.getUnpairedBreaksSoft();
    auto eventTotal2 = bp2.getPairedBreaksSoft() + bp2.getPairedBreaksHard() +
                       bp2.getUnpairedBreaksSoft();
    if (bp1.getBreaksShortIndel() > eventTotal1 ||
        bp2.getBreaksShortIndel() > eventTotal2) {
        return 28;
    }
    if (clonalityStatus1 == EXTREME_SUBCLONAL) {
        return 39;
    }
    if (clonalityStatus2 == EXTREME_SUBCLONAL) {
        return 40;
    }
    if (artifactStatus == ARTIFACT) {
        return 41;
    }
    if (NOCONTROLMODE || germlineStatus1 || germlineStatus2) {
        if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
            return 42;
        }
    }
    if (mrefHits1 > BPFREQTHRESHOLD) {
        if (chrIndex1 != 999 || chrIndex2 == 999 ||
            mrefHits2 > BPFREQTHRESHOLD) {
            return 431;
        }
        if (mrefHits1 > RELAXEDBPFREQTHRESHOLD || NOCONTROLMODE ||
            germlineStatus1 || germlineStatus2 || mrefHits2 > GERMLINEDBLIMIT) {
            return 431;
        }
    }
    if (mrefHits2 > BPFREQTHRESHOLD) {
        if (chrIndex1 == 999 || chrIndex2 != 999 ||
            mrefHits1 > BPFREQTHRESHOLD) {
            return 432;
        }
        if (mrefHits2 > RELAXEDBPFREQTHRESHOLD || NOCONTROLMODE ||
            germlineStatus1 || germlineStatus2 || mrefHits1 > GERMLINEDBLIMIT) {
            return 432;
        }
    }
    return suspicious;
}

int
SvEvent::filterMatchUnknown(const BreakpointReduced &bp1) {
    if (suspicious != 0) {
        return suspicious;
    }
    if (selectedSa1.isSuspicious() || selectedSa1.isSemiSuspicious()) {
        return 1;
    }
    if (selectedSa1.getExpectedDiscordants() > 0 && mateRatio1 < 0.1) {
        return 2;
    }
    auto eventTotal1 = bp1.getPairedBreaksSoft() + bp1.getPairedBreaksHard() +
                       bp1.getUnpairedBreaksSoft();
    if (eventTotal1 + selectedSa1.getMateSupport() + bp1.getNormalSpans() <
        10) {
        return 29;
    }
    if (mrefHits1 > GERMLINEDBLIMIT &&
        !(bp1.getChrIndex() == 999 || bp1.getChrIndex() == 1000)) {
        return 301;
    }
    if (mrefHits2 > GERMLINEDBLIMIT && !(selectedSa1.getChrIndex() == 999 ||
                                         selectedSa1.getChrIndex() == 1000)) {
        return 302;
    }
    if (!distant) {
        if (inverted || totalEvidence1 < 5 || evidenceLevel1 == 1) {
            return 30;
        }
        if (eventSize > 0 && eventSize < HALFDEFAULTREADLENGTH) {
            return 31;
        }
        if (eventType == 3) {
            return 31;
        }
        if (selectedSa1.isFuzzy() || selectedSa1.isSemiSuspicious() ||
            selectedSa1.isStrictFuzzyCandidate()) {
            return 31;
        }
        if (mrefHits1 > GERMLINEDBLIMIT || mrefHits2 > GERMLINEDBLIMIT) {
            if (totalEvidence1 < 5) {
                return 31;
            }
        }
    } else {
        if (selectedSa1.getSupport() < 3 &&
            selectedSa1.getSecondarySupport() < 3 &&
            selectedSa1.getMateSupport() < 3) {
            return 18;
        }
        if (selectedSa1.getMateSupport() < 3) {
            if ((mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) ||
                (evidenceLevel1 < 3 && selectedSa1.getSupport() < 5)) {
                return 32;
            }
        }
        if (mateRatio1 < 0.33 ||
            (evidenceLevel1 != 3 && totalEvidence1 < 5 && mateRatio1 < 0.5)) {
            return 33;
        }
        if (selectedSa1.isFuzzy() || selectedSa1.isStrictFuzzyCandidate()) {
            if (selectedSa1.getSupport() == 0 &&
                selectedSa1.getSecondarySupport() == 0) {
                return 36;
            }
        }
    }
    if ((bp1.getBreaksShortIndel() > eventTotal1)) {
        return 38;
    }
    if (clonalityStatus1 == EXTREME_SUBCLONAL) {
        return 44;
    }
    if (artifactStatus == ARTIFACT) {
        return 45;
    }
    if (NOCONTROLMODE || (germlineStatus1 || germlineStatus2)) {
        if (NOCONTROLMODE) {
            if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
                return 46;
            }
        } else {
            if (mrefHits1 > GERMLINEDBLIMIT || mrefHits2 > GERMLINEDBLIMIT) {
                return 47;
            }
        }
    }
    if (mrefHits1 > BPFREQTHRESHOLD) {
        if (chrIndex1 != 999 || selectedSa1.getChrIndex() == 999 ||
            mrefHits2 > BPFREQTHRESHOLD) {
            return 471;
        }
        if (mrefHits1 > 3 * BPFREQTHRESHOLD || NOCONTROLMODE ||
            germlineStatus1 || germlineStatus2 || mrefHits2 > GERMLINEDBLIMIT) {
            return 471;
        }
    }
    if (mrefHits2 > BPFREQTHRESHOLD) {
        if (chrIndex2 == 999 || selectedSa1.getChrIndex() != 999 ||
            mrefHits1 > BPFREQTHRESHOLD) {
            return 472;
        }
        if (mrefHits2 > 3 * BPFREQTHRESHOLD || NOCONTROLMODE ||
            germlineStatus1 || germlineStatus2 || mrefHits1 > GERMLINEDBLIMIT) {
            return 472;
        }
    }
    return suspicious;
}

pair<double, double>
SvEvent::assessSvClonality(const BreakpointReduced &bp,
                           int eventSupportTotal) const {
    auto artifactTotal1 = bp.getLowQualSpansSoft() + bp.getLowQualBreaksSoft() +
                          bp.getRepetitiveOverhangBreaks();
    auto eventTotal1 = eventSupportTotal + bp.getUnpairedBreaksSoft();
    auto artifactRatio =
        (artifactTotal1 + 0.0) / (artifactTotal1 + eventTotal1);
    auto clonalityRatio =
        (eventTotal1 + 0.0) / (eventTotal1 + bp.getNormalSpans());
    return {artifactRatio, clonalityRatio};
}

ClonalityStatus
SvEvent::assessBreakpointClonalityStatus(double clonalityRatioIn,
                                         const BreakpointReduced &bp1,
                                         const BreakpointReduced &bp2) const {
    if (clonalityRatioIn < CLONALITYLOWTHRESHOLD) {
        if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
            return EXTREME_SUBCLONAL;
        }
        if (distant) {
            if ((selectedSa1.getMateSupport() > 9 && mateRatio1 >= 0.4) ||
                (selectedSa1.getMateSupport() > 4 && mateRatio1 >= 0.6)) {
                return SUBCLONAL;
            }
            if ((selectedSa2.getMateSupport() > 9 && mateRatio2 >= 0.4) ||
                (selectedSa2.getMateSupport() > 4 && mateRatio2 >= 0.6)) {
                return SUBCLONAL;
            }
        }
        return EXTREME_SUBCLONAL;
    } else if (clonalityRatioIn >= CLONALITYHIGHTHRESHOLD) {
        return HOMO;
    } else {
        return HETERO;
    }
}

ClonalityStatus
SvEvent::assessBreakpointClonalityStatusSingle(
    double clonalityRatioIn, const BreakpointReduced &bp1,
    const BreakpointReduced &bp2) const {
    if (clonalityRatioIn < CLONALITYLOWTHRESHOLD) {
        if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
            return EXTREME_SUBCLONAL;
        }
        if (distant &&
            ((selectedSa1.getMateSupport() > 9 && mateRatio1 >= 0.5) ||
             (selectedSa1.getMateSupport() > 4 && mateRatio1 >= 0.8))) {
            return SUBCLONAL;
        } else {
            return EXTREME_SUBCLONAL;
        }
    } else if (clonalityRatioIn >= CLONALITYHIGHTHRESHOLD) {
        return HOMO;
    } else {
        return HETERO;
    }
}

ClonalityStatus
SvEvent::assessBreakpointClonalityStatusUnknown(
    double clonalityRatioIn, const BreakpointReduced &bp1) const {
    if (clonalityRatioIn < CLONALITYLOWTHRESHOLD) {
        if (mrefHits1 > GERMLINEDBLIMIT) {
            return EXTREME_SUBCLONAL;
        }
        if (distant &&
            ((selectedSa1.getMateSupport() > 9 && mateRatio1 >= 0.5) ||
             (selectedSa1.getMateSupport() > 4 && mateRatio1 >= 0.8))) {
            return SUBCLONAL;
        } else {
            return EXTREME_SUBCLONAL;
        }
    } else if (clonalityRatioIn >= CLONALITYHIGHTHRESHOLD) {
        return HOMO;
    } else {
        return HETERO;
    }
}

void
SvEvent::assessSvArtifactStatus(const BreakpointReduced &bp1,
                                const BreakpointReduced &bp2) {
    if (artifactRatio1 < ARTIFACTFREQLOWTHRESHOLD &&
        artifactRatio2 < ARTIFACTFREQLOWTHRESHOLD) {
        artifactStatus = CLEAN;
        return;
    }
    if (artifactRatio1 > 0.85 && artifactRatio2 > 0.85) {
        artifactStatus = ARTIFACT;
        return;
    }
    if (artifactRatio1 > 0.85 || artifactRatio2 > 0.85) {
        if (!(evidenceLevel1 >= 2 || evidenceLevel2 >= 2)) {
            artifactStatus = ARTIFACT;
            return;
        }
    }
    if (artifactRatio1 > ARTIFACTFREQHIGHTHRESHOLD &&
        artifactRatio2 > ARTIFACTFREQHIGHTHRESHOLD &&
        (mrefHits1 > GERMLINEDBLIMIT || mrefHits1 > GERMLINEDBLIMIT)) {
        artifactStatus = ARTIFACT;
        return;
    }
    if ((artifactRatio1 > ARTIFACTFREQHIGHTHRESHOLD ||
         artifactRatio2 > ARTIFACTFREQHIGHTHRESHOLD) &&
        (clonalityStatus1 == EXTREME_SUBCLONAL ||
         clonalityStatus2 == EXTREME_SUBCLONAL)) {
        artifactStatus = ARTIFACT;
        return;
    }
    artifactStatus = BORDERLINE;
}

void
SvEvent::assessSvArtifactStatusUnknown() {
    if (artifactRatio1 < ARTIFACTFREQLOWTHRESHOLD &&
        artifactRatio2 < ARTIFACTFREQLOWTHRESHOLD) {
        artifactStatus = CLEAN;
        return;
    }
    if (artifactRatio1 > ARTIFACTFREQHIGHTHRESHOLD) {
        artifactStatus = ARTIFACT;
        return;
    }
    artifactStatus = BORDERLINE;
}

int
SvEvent::assessEventScore(bool hardClipSuspiciousCall, int inputScoreCategory) {
    // inputscore 0: UNKNOWN partner
    // inputscore 1: known partner with no SA signal
    // inputscore 2: known partner with matching SA signal
    if (inputScoreCategory == 0) {
        if (distant) {
            if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
                if (selectedSa1.getSupport() > 30 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 10) {
                    eventType = 5;
                    return 1;
                }
                if (selectedSa1.getSupport() > 100 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 20) {
                    eventType = 5;
                    return 1;
                }
            } else if (mrefHits1 > GERMLINEDBLIMIT) {
                if (selectedSa1.getSupport() > 50 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 10) {
                    eventType = 5;
                    return 1;
                }
                if (selectedSa1.getSupport() > 100 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 20) {
                    eventType = 5;
                    return 1;
                }
            }
            if ((selectedSa1.isFuzzy() ||
                 selectedSa1.isStrictFuzzyCandidate())) {
                if (evidenceLevel1 == 3 ||
                    !selectedSa1.isStrictFuzzyCandidate()) {
                    if ((selectedSa1.getMateSupport() > 9 &&
                         mateRatio1 >= 0.5)   //
                        || (selectedSa1.getMateSupport() > 5 &&
                            (selectedSa1.getExpectedDiscordants() -
                             selectedSa1.getMateSupport()) < 3)) {
                        if (mrefHits1 == 0 || mrefHits2 == 0) {
                            return 4;
                        }
                    }
                }
                return 1;
            }
            if (!hardClipSuspiciousCall &&
                !selectedSa1.isStrictFuzzyCandidate()) {
                if (selectedSa1.getSupport() > 4 || evidenceLevel1 == 3 ||
                    selectedSa1.getMateSupport() > 9) {
                    return 4;
                }
            }
        }
        return 1;
    } else if (inputScoreCategory == 2) {
        //		auto messageMode = selectedSa1.getChrIndex() == 11 &&
        //(selectedSa1.getPos() == 2261373 && selectedSa2.getPos() == 2148480);
        if (!distant) {
            if (totalEvidence1 < 5 && totalEvidence2 < 5) {
                if (mrefHits1 > GERMLINEDBLIMIT ||
                    mrefHits2 > GERMLINEDBLIMIT) {
                    return 1;
                }
                return 3;
            }
        } else {
            if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
                if (selectedSa1.getSupport() > 30 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 10) {
                    eventType = 5;
                    return 2;
                }
                if (selectedSa1.getSupport() > 100 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 20) {
                    eventType = 5;
                    return 2;
                }
                if (selectedSa2.getSupport() > 30 &&
                    selectedSa2.getSecondarySupport() < 10 &&
                    selectedSa2.getMateSupport() < 10) {
                    eventType = 5;
                    return 2;
                }
                if (selectedSa2.getSupport() > 100 &&
                    selectedSa2.getSecondarySupport() < 10 &&
                    selectedSa2.getMateSupport() < 20) {
                    eventType = 5;
                    return 2;
                }
            } else if (mrefHits1 > GERMLINEDBLIMIT) {
                if (selectedSa1.getSupport() > 50 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 10) {
                    eventType = 5;
                    return 2;
                }
                if (selectedSa1.getSupport() > 100 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 20) {
                    eventType = 5;
                    return 2;
                }

            } else if (mrefHits2 > GERMLINEDBLIMIT) {
                if (selectedSa2.getSupport() > 50 &&
                    selectedSa2.getSecondarySupport() < 10 &&
                    selectedSa2.getMateSupport() < 10) {
                    eventType = 5;
                    return 2;
                }
                if (selectedSa2.getSupport() > 100 &&
                    selectedSa2.getSecondarySupport() < 10 &&
                    selectedSa2.getMateSupport() < 20) {
                    eventType = 5;
                    return 2;
                }
            }
            if (semiSuspicious) {
                if (selectedSa1.isSemiSuspicious() && evidenceLevel1 < 3) {
                    return 2;
                }
                if (selectedSa1.isSemiSuspicious() &&
                    (selectedSa1.getSecondarySupport() == 0 ||
                     selectedSa1.getSupport() == 0) &&
                    evidenceLevel2 == 0) {
                    return 2;
                }
                if (selectedSa2.isSemiSuspicious()) {
                    if (evidenceLevel1 < 3) {
                        auto oneSidedScore =
                            assessEventScore(hardClipSuspiciousCall, 0);
                        if (oneSidedScore > 2) {
                            semiSuspicious = false;
                            return oneSidedScore;
                        } else {
                            return 2;
                        }
                    }
                }
            }
            if (selectedSa1.isFuzzy() && selectedSa2.isFuzzy()) {
                if (!semiSuspicious && mrefHits1 == 0 && mrefHits2 == 0 &&
                    !(selectedSa1.isProperPairErrorProne() ||
                      selectedSa2.isProperPairErrorProne())) {
                    if (evidenceLevel1 > 1 || evidenceLevel2 > 1) {
                        if (mateRatio1 > 0.5 && mateRatio2 > 0.5) {
                            return 3;
                        }
                    }
                }
                if (evidenceLevel1 < 3 || evidenceLevel2 < 3) {
                    return 2;
                }
            } else {
                if (selectedSa1.isFuzzy()) {
                    if (evidenceLevel1 < 2 && evidenceLevel2 < 2) {
                        return 2;
                    }
                }
                if (selectedSa2.isFuzzy()) {
                    if (evidenceLevel2 < 2 && evidenceLevel1 < 2) {
                        return 2;
                    }
                }
            }
            if (evidenceLevel1 < 3 && evidenceLevel2 < 3) {
                if ((selectedSa1.getMateSupport() < 10 && mateRatio1 < 0.5) ||
                    (selectedSa2.getMateSupport() < 10 && mateRatio2 < 0.5)) {
                    return 2;
                }
                if (selectedSa1.getMateSupport() < 10 ||
                    selectedSa2.getMateSupport() < 10) {
                    if ((mateRatio1 + mateRatio2) < 1.1) {
                        return 2;
                    }
                }
                if (!(evidenceLevel1 > 1 || evidenceLevel2 > 1)) {
                    if (!semiSuspicious && mrefHits1 == 0 && mrefHits2 == 0 &&
                        !(selectedSa1.isProperPairErrorProne() ||
                          selectedSa2.isProperPairErrorProne())) {
                        if (mateRatio1 >= 0.8 && mateRatio2 >= 0.8) {
                            if (selectedSa1.getMateSupport() > 4 &&
                                selectedSa2.getMateSupport() > 4) {
                                if (selectedSa1.getMateSupport() > 9 ||
                                    selectedSa2.getMateSupport() > 9) {
                                    return 3;
                                }
                            }
                        }
                    }
                    return 2;
                }
                if (evidenceLevel2 == 1) {
                    if ((mateRatio2 >= 0.8 &&
                         selectedSa2.getMateSupport() > 9)) {
                        auto oneSidedScore =
                            assessEventScore(hardClipSuspiciousCall, 0);
                        if (oneSidedScore > 2) {
                            return oneSidedScore;
                        }
                    }
                    return 2;
                }
                if (evidenceLevel1 == 1) {
                    if (mateRatio1 >= 0.8 && selectedSa1.getMateSupport() > 9 &&
                        mateRatio2 >= 0.6) {
                        return 5;
                    }
                    return 2;
                }
                if (!semiSuspicious && mrefHits1 < GERMLINEDBLIMIT &&
                    mrefHits2 < GERMLINEDBLIMIT &&
                    !(selectedSa1.isProperPairErrorProne() ||
                      selectedSa2.isProperPairErrorProne())) {
                    if (mateRatio1 >= 0.8 && mateRatio2 >= 0.8) {
                        if (selectedSa1.getMateSupport() > 4 &&
                            selectedSa2.getMateSupport() > 4) {
                            return 5;
                        }
                    }
                }
                return 2;
            }
        }
        return 5;
    } else if (inputScoreCategory == 1) {
        if (distant) {
            if (mrefHits1 > GERMLINEDBLIMIT && mrefHits2 > GERMLINEDBLIMIT) {
                if (selectedSa1.getSupport() > 30 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 10) {
                    eventType = 5;
                    return 1;
                }
                if (selectedSa1.getSupport() > 100 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 20) {
                    eventType = 5;
                    return 1;
                }
            } else if (mrefHits1 > GERMLINEDBLIMIT) {
                if (selectedSa1.getSupport() > 50 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 10) {
                    eventType = 5;
                    return 1;
                }
                if (selectedSa1.getSupport() > 100 &&
                    selectedSa1.getSecondarySupport() < 10 &&
                    selectedSa1.getMateSupport() < 20) {
                    eventType = 5;
                    return 1;
                }
            }
            if (selectedSa1.isFuzzy() || selectedSa1.isStrictFuzzyCandidate()) {
                if (evidenceLevel1 == 3) {
                    return 3;
                }
                return 0;
            }
        } else {
            if (totalEvidence1 < 5) {
                return 1;
            }
            if (selectedSa1.getSupport() < 2 ||
                selectedSa1.getSecondarySupport() < 2) {
                return 1;
            }
        }
        if (!hardClipSuspiciousCall && !selectedSa1.isStrictFuzzyCandidate()) {
            if (evidenceLevel1 == 3) {
                return 3;
            } else {
                return 1;
            }
        }
        return 0;
    }
    return 0;
}

void
SvEvent::assessContamination(const vector<pair<int, string>> &overhangDb) {
    if (eventType == 5) {
        contaminationCandidate = 2;
        return;
    }
    auto score1 = 0;
    auto score2 = 0;
    if (inputScore == 2) {
        if (overhang2Index != -1) {
            auto res = assessContaminationSingleBp(overhang2Index, overhangDb,
                                                   selectedSa2);
            score2 = res.first;
            overhang2lengthRatio = res.second;
        }
    }
    if (contaminationCandidate < 2) {
        if (overhang1Index != -1) {
            auto res = assessContaminationSingleBp(overhang1Index, overhangDb,
                                                   selectedSa1);
            score1 = res.first;
            overhang1lengthRatio = res.second;
        }
    }
    if (score1 > 1 || score2 > 1) {
        contaminationCandidate = 2;
    } else if (score1 == 1 && score2 == 1) {
        contaminationCandidate = 1;
    } else if (score1 == 1 && score2 == 0) {
        contaminationCandidate = 0;
    } else if (score1 == 0 && score2 == 1) {
        contaminationCandidate = 0;
    }
}

pair<int, double>
SvEvent::assessContaminationSingleBp(
    int overhangIndex, const vector<pair<int, string>> &overhangDb,
    const SuppAlignmentAnno &selectedSa) {
    auto overhangLengthMax = 0;
    auto overhangLength = 0;
    for (auto cit = overhangDb[overhangIndex].second.cbegin();
         cit != overhangDb[overhangIndex].second.cend(); ++cit) {
        switch (*cit) {
        case '(':
            overhangLengthMax = max(overhangLengthMax, overhangLength);
            overhangLength = 0;
            break;
        case ':':
            overhangLength = 0;
            break;
        default:
            ++overhangLength;
            break;
        }
    }
    auto maxOverhangLengthRatio =
        (overhangLengthMax + 0.0) / SuppAlignmentAnno::DEFAULTREADLENGTH;
    if (selectedSa.getSecondarySupport() > 4) {
        return {0, maxOverhangLengthRatio};
    }
    if (maxOverhangLengthRatio > 0.7) {
        if (selectedSa.getSecondarySupport() == 0) {
            return {3, maxOverhangLengthRatio};
        }
        if (selectedSa.getSecondarySupport() < 3 ||
            (inputScore < 2 && selectedSa.getSecondarySupport() < 10)) {
            return {2, maxOverhangLengthRatio};
        }
    } else if (maxOverhangLengthRatio > 0.6) {
        if (selectedSa.getSecondarySupport() < 3) {
            return {1, maxOverhangLengthRatio};
        }
    }
    return {0, maxOverhangLengthRatio};
}

string
SvEvent::printMatch(const vector<pair<int, string>> &overhangDb) const {
    vector<string> outputFields;
    outputFields.reserve(20);
    outputFields.emplace_back(ChrConverter::indexToChr[chrIndex1]);
    outputFields.emplace_back(strtk::type_to_string<int>(pos1 - 1));
    outputFields.emplace_back(strtk::type_to_string<int>(pos1));
    outputFields.emplace_back(ChrConverter::indexToChr[chrIndex2]);
    outputFields.emplace_back(strtk::type_to_string<int>(pos2 - 1));
    outputFields.emplace_back(
        inputScore > 0
            ? strtk::type_to_string<int>(pos2)
            : strtk::type_to_string<int>(selectedSa1.getExtendedPos()));

    if (!germlineStatus1) {
        outputFields.emplace_back(
            "SOMATIC(" + strtk::type_to_string<int>(mrefHits1) + "/" +
            PIDSINMREFSTR +
            "):" + boost::str(doubleFormatter % germlineClonality1));
    } else {
        if (germline || (mrefHits1 > GERMLINEDBLIMIT)) {
            outputFields.emplace_back(
                "GERMLINE(" + strtk::type_to_string<int>(mrefHits1) + "/" +
                PIDSINMREFSTR +
                "):" + boost::str(doubleFormatter % germlineClonality1));
        } else {
            outputFields.emplace_back(
                "RESCUED(" + strtk::type_to_string<int>(mrefHits1) + "/" +
                PIDSINMREFSTR +
                "):" + boost::str(doubleFormatter % germlineClonality1));
        }
    }
    if (!germlineStatus2) {
        if (inputScore > 0) {
            outputFields.emplace_back(
                "SOMATIC(" + strtk::type_to_string<int>(mrefHits2) + "/" +
                PIDSINMREFSTR +
                "):" + boost::str(doubleFormatter % germlineClonality2));
        } else {
            outputFields.emplace_back(
                "UNKNOWN(" + strtk::type_to_string<int>(mrefHits2) + "/" +
                PIDSINMREFSTR +
                "):" + boost::str(doubleFormatter % germlineClonality2));
        }

    } else {
        if (germline || (mrefHits2 > GERMLINEDBLIMIT)) {
            outputFields.emplace_back(
                "GERMLINE(" + strtk::type_to_string<int>(mrefHits2) + "/" +
                PIDSINMREFSTR +
                "):" + boost::str(doubleFormatter % germlineClonality2));
        } else {
            outputFields.emplace_back(
                "RESCUED(" + strtk::type_to_string<int>(mrefHits2) + "/" +
                PIDSINMREFSTR +
                "):" + boost::str(doubleFormatter % germlineClonality2));
        }
    }
    outputFields.emplace_back(EVENTTYPES[eventType]);
    outputFields.emplace_back((suspicious == 0)
                                  ? strtk::type_to_string<int>(eventScore)
                                  : strtk::type_to_string<int>(-suspicious));
    outputFields.emplace_back(
        (eventSize > 0) ? strtk::type_to_string<int>(eventSize) : "NA");
    outputFields.emplace_back(inverted ? "INV" : "NORMAL");

    outputFields.emplace_back(strtk::type_to_string<int>(totalEvidence1));
    outputFields.emplace_back(boost::str(
        doubleFormatter % (totalEvidence1 / (totalEvidence1 + span1 + 0.0))));
    if (inputScore > 0) {
        outputFields.emplace_back(strtk::type_to_string<int>(totalEvidence2));
        outputFields.emplace_back(
            (totalEvidence2 == 0 && span2 == 0)
                ? "0.000"
                : boost::str(
                      doubleFormatter %
                      (totalEvidence2 / (totalEvidence2 + span2 + 0.0))));
    } else {
        outputFields.emplace_back("UNKNOWN");
        outputFields.emplace_back("UNKNOWN");
    }

    outputFields.emplace_back(selectedSa1.print());
    outputFields.emplace_back(inputScore == 2 ? selectedSa2.print() : "_");

    outputFields.emplace_back(
        overhang1Index != -1 ? overhangDb[overhang1Index].second : ".");
    outputFields.emplace_back(
        overhang2Index != -1 ? overhangDb[overhang2Index].second : ".");

    return collapseRange(outputFields, "\t").append("\n");
}
// vector<int> SvEvent::getKey() const {
//	if (!DEBUGMODE && (suspicious != 0 || eventScore == 0)) {
//		return {};
//	}
//	auto keyScore = (suspicious == 0) ? eventScore : suspicious;
//	if (chrIndex1 < chrIndex2 || (chrIndex1 == chrIndex2 && pos1 < pos2)) {
//		return {chrIndex1,pos1,chrIndex2,pos2,keyScore};
//	} else {
//		return {chrIndex2,pos2,chrIndex1,pos1,keyScore};
//	}
// }

string
SvEvent::getKey() const {
    if (!DEBUGMODE && (suspicious != 0 || eventScore == 0)) {
        return {};
    }
    auto keyScore = (suspicious == 0) ? eventScore : suspicious;
    if (chrIndex1 < chrIndex2 || (chrIndex1 == chrIndex2 && pos1 < pos2)) {
        return collapseRange({to_string(chrIndex1), to_string(pos1),
                              to_string(chrIndex2), to_string(pos2),
                              to_string(keyScore)},
                             "_");
    } else {
        return collapseRange({to_string(chrIndex2), to_string(pos2),
                              to_string(chrIndex1), to_string(pos1),
                              to_string(keyScore)},
                             "_");
    }
}

} /* namespace sophia */
