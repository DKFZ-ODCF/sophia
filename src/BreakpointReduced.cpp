/*
 * MrefEntry.cpp
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

#include "global.h"
#include "Breakpoint.h"
#include "GlobalAppConfig.h"
#include "strtk-wrap.h"
#include <BreakpointReduced.h>
#include <boost/algorithm/string/join.hpp>
#include <unordered_set>

namespace sophia {

using namespace std;

boost::format BreakpointReduced::doubleFormatter{"%.3f"};
int BreakpointReduced::DEFAULTREADLENGTH{};
double BreakpointReduced::CLONALITYSTRICTLOWTHRESHOLD{};
double BreakpointReduced::ARTIFACTFREQHIGHTHRESHOLD{};
string BreakpointReduced::PIDSINMREFSTR{};

sophia::BreakpointReduced::BreakpointReduced(const Breakpoint &tmpBp,
                                             int lineIndexIn,
                                             bool hasOverhangIn)
    : hasOverhang{hasOverhangIn},
      toRemove{false},
      lineIndex{lineIndexIn},
      chrIndex{tmpBp.getChrIndex()},
      pos{tmpBp.getPos()},
      normalSpans{tmpBp.getNormalSpans()},
      lowQualSpansSoft{tmpBp.getLowQualBreaksSoft()},
      lowQualSpansHard{tmpBp.getLowQualSpansHard()},
      unpairedBreaksSoft{tmpBp.getUnpairedBreaksSoft()},
      unpairedBreaksHard{tmpBp.getUnpairedBreaksHard()},
      breaksShortIndel{tmpBp.getBreaksShortIndel()},
      lowQualBreaksSoft{tmpBp.getLowQualBreaksSoft()},
      lowQualBreaksHard{tmpBp.getLowQualBreaksHard()},
      repetitiveOverhangBreaks{tmpBp.getRepetitiveOverhangBreaks()},
      pairedBreaksSoft{tmpBp.getPairedBreaksSoft()},
      pairedBreaksHard{tmpBp.getPairedBreaksHard()},
      mateSupport{tmpBp.getMateSupport()},
      leftCoverage{tmpBp.getLeftCoverage()},
      rightCoverage{tmpBp.getRightCoverage()},
      mrefHits{MrefMatch{-1, -1, 10000, {}, }},
      germlineInfo{GermlineMatch{0.0, 0.0, {}, }},
      suppAlignments{} {
    const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
    for (const auto &sa : tmpBp.getDoubleSidedMatches()) {
        if (!chrConverter.isIgnoredChromosome(sa.getChrIndex())) {
            suppAlignments.emplace_back(sa);
        }
    }
    for (const auto &sa : tmpBp.getSupplementsPrimary()) {
        if (!chrConverter.isIgnoredChromosome(sa.getChrIndex())) {
            suppAlignments.emplace_back(sa);
        }
    }
    complexRearrangementMateRatioRescue(true);
    complexRearrangementMateRatioRescue(false);
}

void
BreakpointReduced::complexRearrangementMateRatioRescue(bool encounteredM) {
    auto candidateCount = 0;
    auto cumulativeMateSupport = 0.0;
    auto maxExpectedDiscordants = 0;
    for (const auto &sa : suppAlignments) {
        if (sa.isDistant() && sa.isEncounteredM() == encounteredM &&
            !sa.isSuspicious() && sa.getMateSupport() > 4) {
            ++candidateCount;
            if (candidateCount == 3) {
                return;
            }
            cumulativeMateSupport += sa.getMateSupport();
            maxExpectedDiscordants =
                max(maxExpectedDiscordants, sa.getExpectedDiscordants());
        }
    }
    if (candidateCount == 2 &&
        cumulativeMateSupport / maxExpectedDiscordants > 0.7) {
        for (auto &sa : suppAlignments) {
            if (sa.isDistant() && sa.isEncounteredM() == encounteredM &&
                !sa.isSuspicious() && sa.getMateSupport() > 4) {
                sa.setExpectedDiscordants(sa.getMateSupport());
            }
        }
    }
}

sophia::BreakpointReduced::BreakpointReduced(
    const SuppAlignmentAnno &sa,
    const BreakpointReduced &emittingBp,
    bool fuzzySecondary)
    : hasOverhang{false},
      toRemove{false},
      lineIndex{-1},
      chrIndex{sa.getChrIndex()},
      pos{!fuzzySecondary ? sa.getPos(): sa.getExtendedPos()},
      normalSpans{},
      lowQualSpansSoft{},
      lowQualSpansHard{},
      unpairedBreaksSoft{},
      unpairedBreaksHard{},
      breaksShortIndel{},
      lowQualBreaksSoft{},
      lowQualBreaksHard{},
      repetitiveOverhangBreaks{},
      pairedBreaksSoft{},
      pairedBreaksHard{},
      mateSupport{},
      leftCoverage{},
      rightCoverage{},
      mrefHits{MrefMatch{-1, -1, 10000, {}}},
      germlineInfo{GermlineMatch{0.0, 0.0, {}}},
      suppAlignments{} {
    addDummySa(sa, emittingBp);
}

void
sophia::BreakpointReduced::addDummySa(const SuppAlignmentAnno &sa,
                                      const BreakpointReduced &emittingBp) {
    suppAlignments.emplace_back(emittingBp.getChrIndex(), emittingBp.getPos(),
                                sa);
}

const SuppAlignmentAnno &
sophia::BreakpointReduced::getDummySa() {
    return suppAlignments.back();
}

SuppAlignmentAnno *
BreakpointReduced::searchFuzzySa(const SuppAlignmentAnno &fuzzySa) {
    SuppAlignmentAnno *match = nullptr;
    for (auto &sa : suppAlignments) {
        if (sa.saClosenessDirectional(fuzzySa, DEFAULTREADLENGTH * 0.2)) {
            match = &sa;
            return match;
        }
    }
    return nullptr;
}

bool
BreakpointReduced::testOverhangBasedCandidacy() const {
    if (pairedBreaksSoft > 0) {
        return false;
    }
    if (breaksShortIndel > 0) {
        return false;
    }
    if (unpairedBreaksSoft < 5) {
        return false;
    }
    if (((0.0 + unpairedBreaksSoft) / normalSpans) <
        CLONALITYSTRICTLOWTHRESHOLD) {
        return false;
    }
    auto artifactTotal =
        0.0 + lowQualSpansSoft + lowQualBreaksSoft + repetitiveOverhangBreaks;
    if ((artifactTotal / (unpairedBreaksSoft + artifactTotal)) >
        ARTIFACTFREQHIGHTHRESHOLD) {
        return false;
    }
    return true;
}

string
BreakpointReduced::printOverhang(double germlineClonality,
                                 int numHits,
                                 const string &overhang) const {
    string res{"##"};
    res.append(GlobalAppConfig::getInstance().getChrConverter().indexToChrName(chrIndex)).append("\t");
    res.append(strtk::type_to_string<int>(pos - 1)).append("\t");
    res.append(strtk::type_to_string<int>(pos)).append("\t");
    if (germlineClonality > 0.1) {
        res.append("GERMLINE(");
    } else {
        res.append("SOMATIC(");
    }
    res.append(strtk::type_to_string<int>(numHits))
        .append("/")
        .append(PIDSINMREFSTR)
        .append("):");
    res.append(boost::str(doubleFormatter % germlineClonality)).append("\t");
    res.append(overhang).append("\n");
    return res;
}

void
BreakpointReduced::removeMarkedFuzzies() {
    suppAlignments.erase(
        remove_if(suppAlignments.begin(), suppAlignments.end(),
                  [](const SuppAlignmentAnno &sa) { return sa.isToRemove(); }),
        suppAlignments.end());
}

} /* namespace sophia */
