/*
 * Breakpoint.cpp
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

#include "Breakpoint.h"
#include "ChrConverter.h"
#include "strtk.hpp"
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <cmath>
#include <iostream>
#include <memory>

namespace sophia {

using namespace std;

const string Breakpoint::COLUMNSSTR = boost::join(
    vector<string>{
        //
        "#chr", "start", "end",   //
        "covTypes(pairedBreaksSoft,pairedBreaksHard,mateReadSupport,"
        "unpairedBreaksSoft,unpairedBreaksHard,shortIndelReads,normalSpans,"
        "lowQualSpansSoft,lowQualSpansHard,lowQualBreaksSoft,lowQualBreaksHard,"
        "repetitiveOverhangs)",         //
        "leftCoverage,rightCoverage",   //
        "suppAlignmentsDoubleSupport(primary,secondary,mate)",
        "suppAlignments(primary,0,mate)", "significantOverhangs\n"},
    "\t");
int Breakpoint::BPSUPPORTTHRESHOLD{};
int Breakpoint::DEFAULTREADLENGTH{};
int Breakpoint::DISCORDANTLOWQUALLEFTRANGE{};
int Breakpoint::DISCORDANTLOWQUALRIGHTRANGE{};
double Breakpoint::IMPROPERPAIRRATIO{0.0};
bool Breakpoint::PROPERPAIRCOMPENSATIONMODE{false};
int Breakpoint::bpindex{0};

Breakpoint::Breakpoint(int chrIndexIn, int posIn)
    : covFinalized{false}, missingInfoBp{false}, chrIndex{chrIndexIn},
      pos{posIn}, normalSpans{0}, lowQualSpansSoft{0}, lowQualSpansHard{0},
      unpairedBreaksSoft{0}, unpairedBreaksHard{0}, breaksShortIndel{0},
      lowQualBreaksSoft{0}, lowQualBreaksHard{0}, repetitiveOverhangBreaks{0},
      pairedBreaksSoft{0}, pairedBreaksHard{0}, leftSideDiscordantCandidates{0},
      rightSideDiscordantCandidates{0}, mateSupport{0}, leftCoverage{0},
      rightCoverage{0}, totalLowMapqHardClips{0},
      hitsInMref{-1}, germline{false}, poolLeft{}, poolRight{},
      poolLowQualLeft{}, poolLowQualRight{} {}

void
Breakpoint::addSoftAlignment(shared_ptr<Alignment> alignmentIn) {
    if (!alignmentIn->isSupplementary()) {
        if (supportingSoftAlignments.size() <= MAXPERMISSIBLESOFTCLIPS) {
            supportingSoftAlignments.push_back(alignmentIn);
        }
    }
}

void
Breakpoint::addHardAlignment(shared_ptr<Alignment> alignmentIn) {
    if (alignmentIn->isSupplementary()) {
        if (!(alignmentIn->isLowMapq() || alignmentIn->isNullMapq())) {
            if (supportingHardAlignments.size() <= MAXPERMISSIBLEHARDCLIPS) {
                supportingHardAlignments.push_back(alignmentIn);
            }
        } else {
            if (totalLowMapqHardClips < MAXPERMISSIBLELOWMAPQHARDCLIPS) {
                supportingHardLowMapqAlignments.push_back(alignmentIn);
                ++totalLowMapqHardClips;
            } else {
                supportingHardLowMapqAlignments.clear();
            }
        }
    }
}

bool
Breakpoint::finalizeBreakpoint(
    const deque<MateInfo> &discordantAlignmentsPool,
    const deque<MateInfo> &discordantLowQualAlignmentsPool,
    const deque<MateInfo> &discordantAlignmentCandidatesPool) {
    auto overhangStr = string();
    auto eventTotal =
        unpairedBreaksSoft + unpairedBreaksHard + breaksShortIndel;
    auto artifactTotal =
        lowQualBreaksSoft + lowQualSpansSoft + lowQualSpansHard;
    if ((eventTotal + artifactTotal > 50) &&
        (artifactTotal / (0.0 + eventTotal + artifactTotal)) > 0.85) {
        ++bpindex;
        missingInfoBp = true;
    } else if (static_cast<int>(supportingSoftAlignments.size()) ==
                   MAXPERMISSIBLESOFTCLIPS &&
               eventTotal + normalSpans + artifactTotal >
                   MAXPERMISSIBLEHARDCLIPS * 20) {
        ++bpindex;
        missingInfoBp = true;
    } else {
        fillMatePool(discordantAlignmentsPool, discordantLowQualAlignmentsPool,
                     discordantAlignmentCandidatesPool);
        if (eventTotal < BPSUPPORTTHRESHOLD &&
            static_cast<int>(poolLeft.size()) < BPSUPPORTTHRESHOLD &&
            static_cast<int>(poolLowQualLeft.size()) < BPSUPPORTTHRESHOLD &&
            static_cast<int>(poolRight.size()) < BPSUPPORTTHRESHOLD &&
            static_cast<int>(poolLowQualRight.size()) < BPSUPPORTTHRESHOLD) {
            if (artifactTotal < normalSpans) {
                return false;
            } else {
                ++bpindex;
                missingInfoBp = true;
            }
        } else {
            overhangStr = finalizeOverhangs();
            detectDoubleSupportSupps();
            collectMateSupport();
        }
    }
    if (eventTotal + mateSupport + artifactTotal < BPSUPPORTTHRESHOLD ||
        (eventTotal + artifactTotal < BPSUPPORTTHRESHOLD &&
         doubleSidedMatches.empty() && supplementsPrimary.empty())) {
        return false;
    }
    if (missingInfoBp ||
        (doubleSidedMatches.empty() && supplementsPrimary.empty() &&
         overhangStr.size() < 2)) {
        auto artifactTotal2 = lowQualBreaksSoft + lowQualSpansHard +
                              lowQualSpansSoft + repetitiveOverhangBreaks;
        auto artifactTotal2Relaxed =
            lowQualBreaksSoft + lowQualSpansSoft + repetitiveOverhangBreaks;
        auto eventTotal2 = pairedBreaksSoft + pairedBreaksHard +
                           unpairedBreaksSoft + unpairedBreaksHard +
                           breaksShortIndel;
        auto eventTotal2Strict =
            pairedBreaksSoft + unpairedBreaksSoft + pairedBreaksHard;
        auto covCriterion = (eventTotal2 + artifactTotal2) > 10;
        if (!(covCriterion && eventTotal2Strict + artifactTotal2Relaxed > 0)) {
            return false;
        }
    }
    printBreakpointReport(overhangStr);
    return true;
}

void
Breakpoint::printBreakpointReport(const string &overhangStr) {
    string res{};
    res.reserve(350);
    res.append(ChrConverter::indexToChr[chrIndex]).append("\t");
    res.append(strtk::type_to_string<int>(pos)).append("\t");
    res.append(strtk::type_to_string<int>(pos + 1)).append("\t");

    res.append(strtk::type_to_string<int>(pairedBreaksSoft)).append(",");
    res.append(strtk::type_to_string<int>(pairedBreaksHard)).append(",");
    res.append(strtk::type_to_string<int>(mateSupport)).append(",");
    res.append(strtk::type_to_string<int>(unpairedBreaksSoft)).append(",");
    res.append(strtk::type_to_string<int>(unpairedBreaksHard)).append(",");
    res.append(strtk::type_to_string<int>(breaksShortIndel)).append(",");

    res.append(strtk::type_to_string<int>(normalSpans)).append(",");

    res.append(strtk::type_to_string<int>(lowQualSpansSoft)).append(",");
    res.append(strtk::type_to_string<int>(lowQualSpansHard)).append(",");
    res.append(strtk::type_to_string<int>(lowQualBreaksSoft)).append(",");
    res.append(strtk::type_to_string<int>(lowQualBreaksHard)).append(",");
    res.append(strtk::type_to_string<int>(repetitiveOverhangBreaks))
        .append("\t");

    res.append(strtk::type_to_string<int>(leftCoverage)).append(",");
    res.append(strtk::type_to_string<int>(rightCoverage)).append("\t");
    if (missingInfoBp) {
        res.append("#\t#\t#\n");
    } else {
        collapseSuppRange(res, doubleSidedMatches);
        res.append("\t");
        collapseSuppRange(res, supplementsPrimary);
        res.append("\t");
        if (overhangStr.empty()) {
            res.append(".").append("\n");
        } else {
            res.append(overhangStr).append("\n");
        }
    }
    cout << res;
}

void
Breakpoint::collapseSuppRange(string &res,
                              const vector<SuppAlignment> &vec) const {
    if (vec.empty()) {
        res.append(".");
    } else {
        for (const auto &suppAlignment : vec) {
            res.append(suppAlignment.print()).append(";");
        }
        res.pop_back();
    }
}

string
Breakpoint::finalizeOverhangs() {
    ++bpindex;
    for (auto i = 0u; i < supportingSoftAlignments.size(); ++i) {
        supportingSoftAlignments[i]->setChosenBp(pos, i);
        if (supportingSoftAlignments[i]->assessOutlierMateDistance()) {
            if (supportingSoftAlignments[i]->getMateChrIndex() < 1002) {
                if (supportingSoftAlignments[i]->isOverhangEncounteredM()) {
                    if (!(supportingSoftAlignments[i]->isNullMapq() ||
                          supportingSoftAlignments[i]->isLowMapq())) {
                        poolLeft.emplace_back(
                            supportingSoftAlignments[i]->getStartPos(),
                            supportingSoftAlignments[i]->getEndPos(),
                            supportingSoftAlignments[i]->getMateChrIndex(),
                            supportingSoftAlignments[i]->getMatePos(), 0,
                            supportingSoftAlignments[i]->isInvertedMate());
                    } else {
                        poolLowQualLeft.emplace_back(
                            supportingSoftAlignments[i]->getStartPos(),
                            supportingSoftAlignments[i]->getEndPos(),
                            supportingSoftAlignments[i]->getMateChrIndex(),
                            supportingSoftAlignments[i]->getMatePos(), 0,
                            supportingSoftAlignments[i]->isInvertedMate());
                    }
                } else {
                    if (!(supportingSoftAlignments[i]->isNullMapq() ||
                          supportingSoftAlignments[i]->isLowMapq())) {
                        poolRight.emplace_back(
                            supportingSoftAlignments[i]->getStartPos(),
                            supportingSoftAlignments[i]->getEndPos(),
                            supportingSoftAlignments[i]->getMateChrIndex(),
                            supportingSoftAlignments[i]->getMatePos(), 0,
                            supportingSoftAlignments[i]->isInvertedMate());
                    } else {
                        poolLowQualRight.emplace_back(
                            supportingSoftAlignments[i]->getStartPos(),
                            supportingSoftAlignments[i]->getEndPos(),
                            supportingSoftAlignments[i]->getMateChrIndex(),
                            supportingSoftAlignments[i]->getMatePos(), 0,
                            supportingSoftAlignments[i]->isInvertedMate());
                    }
                }
            }
        }
    }
    vector<SuppAlignment> supplementsPrimaryTmp{};
    sort(supportingSoftAlignments.begin(), supportingSoftAlignments.end(),
         [](const shared_ptr<Alignment> &a, const shared_ptr<Alignment> &b) {
             return a->getOverhangLength() < b->getOverhangLength();
         });
    vector<shared_ptr<Alignment>> supportingSoftParentAlignments{};
    while (!supportingSoftAlignments.empty()) {
        auto substrCheck = false;
        auto tmpSas = supportingSoftAlignments.back()->generateSuppAlignments(
            chrIndex, pos);
        for (auto overhangParent : supportingSoftParentAlignments) {
            if (matchDetector(overhangParent,
                              supportingSoftAlignments.back())) {
                substrCheck = true;
                overhangParent->addChildNode(
                    supportingSoftAlignments.back()->getOriginIndex());
                overhangParent->addSupplementaryAlignments(tmpSas);
            }
        }
        if (!substrCheck) {
            auto allDistant =
                !tmpSas.empty() &&
                all_of(cbegin(tmpSas), cend(tmpSas),
                       [](const SuppAlignment &sa) { return sa.isDistant(); });
            if (allDistant || supportingSoftAlignments.back()
                                      ->overhangComplexityMaskRatio() <= 0.5) {
                if (supportingSoftAlignments.back()->getOverhangLength() >=
                    20) {
                    supportingSoftAlignments.back()->addSupplementaryAlignments(
                        tmpSas);
                    supportingSoftParentAlignments.push_back(
                        supportingSoftAlignments.back());
                } else {
                    for (const auto &sa : tmpSas) {
                        auto it =
                            find_if(supplementsPrimaryTmp.begin(),
                                    supplementsPrimaryTmp.end(),
                                    [&](const SuppAlignment &suppAlignment) {
                                        return suppAlignment.saCloseness(sa, 5);
                                    });
                        if (it == supplementsPrimaryTmp.end()) {
                            supplementsPrimaryTmp.push_back(sa);
                        } else {
                            if (it->isFuzzy() && !sa.isFuzzy()) {
                                it->removeFuzziness(sa);
                            } else if (it->isFuzzy() && sa.isFuzzy()) {
                                it->extendSuppAlignment(sa.getPos(),
                                                        sa.getExtendedPos());
                            }
                            it->addSupportingIndices(
                                supportingSoftAlignments.back()
                                    ->getChildrenNodes());
                            if (it->isNullMapqSource() &&
                                !supportingSoftAlignments.back()
                                     ->isNullMapq()) {
                                it->setNullMapqSource(false);
                            }
                        }
                    }
                }
            } else {
                --unpairedBreaksSoft;
                ++repetitiveOverhangBreaks;
            }
        }
        supportingSoftAlignments.pop_back();
    }
    string consensusOverhangsTmp{};
    consensusOverhangsTmp.reserve(250);
    {
        auto i = 1;
        auto indexStr = strtk::type_to_string<int>(bpindex);
        for (const auto &overhangParent : supportingSoftParentAlignments) {
            if (static_cast<int>(overhangParent->getChildrenNodes().size()) >=
                BPSUPPORTTHRESHOLD) {
                consensusOverhangsTmp.append(">")
                    .append(indexStr)
                    .append("_")
                    .append(strtk::type_to_string<int>(i))
                    .append(":")
                    .append(overhangParent->printOverhang())
                    .append(";");
                ++i;
            }
            for (const auto &sa :
                 overhangParent->getSupplementaryAlignments()) {
                if (sa.isToRemove()) {
                    continue;
                }
                auto it = find_if(supplementsPrimary.begin(),
                                  supplementsPrimary.end(),
                                  [&](const SuppAlignment &suppAlignment) {
                                      return suppAlignment.saCloseness(sa, 5);
                                  });
                if (it == supplementsPrimary.end()) {
                    supplementsPrimary.push_back(sa);
                    supplementsPrimary.back().addSupportingIndices(
                        overhangParent->getChildrenNodes());
                } else {
                    if (it->isFuzzy() && !sa.isFuzzy()) {
                        it->removeFuzziness(sa);
                    } else if (it->isFuzzy() && sa.isFuzzy()) {
                        it->extendSuppAlignment(sa.getPos(),
                                                sa.getExtendedPos());
                    }
                    it->addSupportingIndices(
                        overhangParent->getChildrenNodes());
                    if (it->isNullMapqSource() &&
                        !supplementsPrimary.back().isNullMapqSource()) {
                        it->setNullMapqSource(false);
                    }
                }
            }
        }
        supportingSoftParentAlignments.clear();
        if (!consensusOverhangsTmp.empty()) {
            consensusOverhangsTmp.pop_back();
        } else {
            return string();
        }
    }
    for (const auto &sa : supplementsPrimaryTmp) {
        auto it = find_if(supplementsPrimary.begin(), supplementsPrimary.end(),
                          [&](const SuppAlignment &suppAlignment) {
                              return suppAlignment.saCloseness(sa, 5);
                          });
        if (it == supplementsPrimary.end()) {
            supplementsPrimary.push_back(sa);
        } else {
            if (it->isFuzzy() && !sa.isFuzzy()) {
                it->removeFuzziness(sa);
            } else if (it->isFuzzy() && sa.isFuzzy()) {
                it->extendSuppAlignment(sa.getPos(), sa.getExtendedPos());
            }
            it->addSupportingIndices(sa.getSupportingIndices());
            if (it->isNullMapqSource() && !sa.isNullMapqSource()) {
                it->setNullMapqSource(false);
            }
        }
    }
    return consensusOverhangsTmp;
}

bool
Breakpoint::matchDetector(const shared_ptr<Alignment> &longAlignment,
                          const shared_ptr<Alignment> &shortAlignment) const {
    if (longAlignment->isOverhangEncounteredM() !=
        shortAlignment->isOverhangEncounteredM()) {
        return false;
    }
    auto mismatches = 0;
    char cLong{};
    char cShort{};
    auto shortS = shortAlignment->getOverhangLength();
    auto longStart = longAlignment->getOverhangStartIndex();
    auto shortStart = shortAlignment->getOverhangStartIndex();
    const auto pointerToLongSeq = &longAlignment->getSamLine();
    const auto pointerToShortSeq = &shortAlignment->getSamLine();
    if (!longAlignment->isOverhangEncounteredM() &&
        !shortAlignment->isOverhangEncounteredM()) {
        auto lenDiff = longAlignment->getOverhangLength() - shortS;
        for (int i = shortS - 1; i >= 0; --i) {
            cLong = (*pointerToLongSeq)[longStart + i + lenDiff];
            if (cLong == 'N')
                continue;
            cShort = (*pointerToShortSeq)[shortStart + i];
            if (cShort == 'N')
                continue;
            if (cLong != cShort) {
                ++mismatches;
                if (mismatches > PERMISSIBLEMISMATCHES)
                    return false;
            }
        }
        return true;
    } else if (longAlignment->isOverhangEncounteredM() &&
               shortAlignment->isOverhangEncounteredM()) {
        for (auto i = 0; i < shortS; ++i) {
            cLong = (*pointerToLongSeq)[longStart + i];
            if (cLong == 'N')
                continue;
            cShort = (*pointerToShortSeq)[shortStart + i];
            if (cShort == 'N')
                continue;
            if (cLong != cShort) {
                ++mismatches;
                if (mismatches > PERMISSIBLEMISMATCHES)
                    return false;
            }
        }
        return true;
    }
    return false;
}

void
Breakpoint::detectDoubleSupportSupps() {
    vector<SuppAlignment> saHardTmpLowQual;
    {
        auto i = 0u;
        for (; i < supportingHardAlignments.size(); ++i) {
            supportingHardAlignments[i]->setChosenBp(pos, i);
        }
        for (auto hardAlignment : supportingHardAlignments) {
            for (const auto &sa :
                 hardAlignment->generateSuppAlignments(chrIndex, pos)) {
                if (!(sa.isInverted() && sa.getPos() == pos &&
                      sa.getChrIndex() == chrIndex)) {
                    supplementsSecondary.push_back(sa);
                    if (supplementsSecondary.back().isDistant()) {
                        if (supplementsSecondary.back().isEncounteredM()) {
                            if (!(hardAlignment->isNullMapq() ||
                                  hardAlignment->isLowMapq())) {
                                poolLeft.emplace_back(
                                    hardAlignment->getStartPos(),
                                    hardAlignment->getEndPos(),
                                    hardAlignment->getMateChrIndex(),
                                    hardAlignment->getMatePos(), 1,
                                    hardAlignment->isInvertedMate());
                            } else {
                                poolLowQualLeft.emplace_back(
                                    hardAlignment->getStartPos(),
                                    hardAlignment->getEndPos(),
                                    hardAlignment->getMateChrIndex(),
                                    hardAlignment->getMatePos(), 1,
                                    hardAlignment->isInvertedMate());
                            }

                        } else {
                            if (!(hardAlignment->isNullMapq() ||
                                  hardAlignment->isLowMapq())) {
                                poolRight.emplace_back(
                                    hardAlignment->getStartPos(),
                                    hardAlignment->getEndPos(),
                                    hardAlignment->getMateChrIndex(),
                                    hardAlignment->getMatePos(), 1,
                                    hardAlignment->isInvertedMate());
                            } else {
                                poolLowQualRight.emplace_back(
                                    hardAlignment->getStartPos(),
                                    hardAlignment->getEndPos(),
                                    hardAlignment->getMateChrIndex(),
                                    hardAlignment->getMatePos(), 1,
                                    hardAlignment->isInvertedMate());
                            }
                        }
                    }
                }
            }
        }
        supportingHardAlignments.clear();
        for (auto j = 0u; j < supportingHardLowMapqAlignments.size(); ++j) {
            supportingHardLowMapqAlignments[j]->setChosenBp(pos, i + j);
        }

        for (auto hardAlignment : supportingHardLowMapqAlignments) {
            for (const auto &sa :
                 hardAlignment->generateSuppAlignments(chrIndex, pos)) {
                if (!(sa.isInverted() && sa.getPos() == pos &&
                      sa.getChrIndex() == chrIndex)) {
                    saHardTmpLowQual.push_back(sa);
                    if (saHardTmpLowQual.back().isDistant()) {
                        if (saHardTmpLowQual.back().isEncounteredM()) {
                            poolLowQualLeft.emplace_back(
                                hardAlignment->getStartPos(),
                                hardAlignment->getEndPos(),
                                hardAlignment->getMateChrIndex(),
                                hardAlignment->getMatePos(), 1,
                                hardAlignment->isInvertedMate());
                        } else {
                            poolLowQualRight.emplace_back(
                                hardAlignment->getStartPos(),
                                hardAlignment->getEndPos(),
                                hardAlignment->getMateChrIndex(),
                                hardAlignment->getMatePos(), 1,
                                hardAlignment->isInvertedMate());
                        }
                    }
                }
            }
        }
        supportingHardLowMapqAlignments.clear();
    }
    {
        vector<int> lowMapqHardSupportWhitelistIndices{};
        for (auto primarySupptIt = supplementsPrimary.begin();
             primarySupptIt != supplementsPrimary.end(); ++primarySupptIt) {
            auto foundMatch = false;
            for (auto secondarySuppIt = supplementsSecondary.begin();
                 secondarySuppIt != supplementsSecondary.end();
                 ++secondarySuppIt) {
                if (primarySupptIt->saCloseness(*secondarySuppIt, 100)) {
                    if (primarySupptIt->isFuzzy() &&
                        !secondarySuppIt->isFuzzy()) {
                        primarySupptIt->removeFuzziness(*secondarySuppIt);
                    } else if (primarySupptIt->isFuzzy() &&
                               secondarySuppIt->isFuzzy()) {
                        primarySupptIt->extendSuppAlignment(
                            secondarySuppIt->getPos(),
                            secondarySuppIt->getExtendedPos());
                    }
                    foundMatch = true;
                    secondarySuppIt->setToRemove(true);
                    primarySupptIt->addSecondarySupportIndices(
                        secondarySuppIt->getSupportingIndicesSecondary());
                }
            }
            for (auto secondarySuppIt = saHardTmpLowQual.begin();
                 secondarySuppIt != saHardTmpLowQual.end(); ++secondarySuppIt) {
                if (primarySupptIt->saCloseness(*secondarySuppIt, 100)) {
                    if (primarySupptIt->isFuzzy() &&
                        !secondarySuppIt->isFuzzy()) {
                        primarySupptIt->removeFuzziness(*secondarySuppIt);
                    } else if (primarySupptIt->isFuzzy() &&
                               secondarySuppIt->isFuzzy()) {
                        primarySupptIt->extendSuppAlignment(
                            secondarySuppIt->getPos(),
                            secondarySuppIt->getExtendedPos());
                    }
                    foundMatch = true;
                    for (auto index :
                         secondarySuppIt->getSupportingIndicesSecondary()) {
                        lowMapqHardSupportWhitelistIndices.push_back(index);
                    }
                    secondarySuppIt->setToRemove(true);
                    primarySupptIt->addSecondarySupportIndices(
                        secondarySuppIt->getSupportingIndicesSecondary());
                }
            }
            if (foundMatch) {
                doubleSidedMatches.push_back(*primarySupptIt);
                primarySupptIt->setToRemove(true);
            }
        }
        cleanUpVector(supplementsPrimary);
        cleanUpVector(supplementsSecondary);
        cleanUpVector(saHardTmpLowQual);
        supplementsSecondary.insert(supplementsSecondary.end(),
                                    saHardTmpLowQual.begin(),
                                    saHardTmpLowQual.end());
        sort(lowMapqHardSupportWhitelistIndices.begin(),
             lowMapqHardSupportWhitelistIndices.end());
        auto whitelistSize =
            distance(lowMapqHardSupportWhitelistIndices.begin(),
                     unique(lowMapqHardSupportWhitelistIndices.begin(),
                            lowMapqHardSupportWhitelistIndices.end()));
        unpairedBreaksHard += whitelistSize;
        lowQualBreaksHard -= whitelistSize;
    }
    auto maxMapq = 0;
    for (const auto &sa : doubleSidedMatches) {
        if (sa.getMapq() > maxMapq) {
            maxMapq = sa.getMapq();
        }
    }
    for (const auto &sa : supplementsPrimary) {
        if (sa.getMapq() > maxMapq) {
            maxMapq = sa.getMapq();
        }
    }
    for (auto &sa : supplementsPrimary) {
        if (sa.getMapq() > 0) {
            if (sa.getMapq() < 13 && sa.getMapq() < maxMapq) {
                sa.setToRemove(true);
            }
        }
    }
    cleanUpVector(supplementsPrimary);
    for (auto &sa : supplementsPrimary) {
        sa.finalizeSupportingIndices();
    }
    vector<int> uniqueDoubleSupportPrimaryIndices{};
    vector<int> uniqueDoubleSupportSecondaryIndices{};
    for (auto &sa : doubleSidedMatches) {
        sa.finalizeSupportingIndices();
        uniqueDoubleSupportPrimaryIndices.insert(
            uniqueDoubleSupportPrimaryIndices.end(),
            sa.getSupportingIndices().cbegin(),
            sa.getSupportingIndices().cend());
        uniqueDoubleSupportSecondaryIndices.insert(
            uniqueDoubleSupportSecondaryIndices.end(),
            sa.getSupportingIndicesSecondary().cbegin(),
            sa.getSupportingIndicesSecondary().cend());
    }
    sort(uniqueDoubleSupportPrimaryIndices.begin(),
         uniqueDoubleSupportPrimaryIndices.end());
    sort(uniqueDoubleSupportSecondaryIndices.begin(),
         uniqueDoubleSupportSecondaryIndices.end());
    auto priCompensation =
        distance(uniqueDoubleSupportPrimaryIndices.begin(),
                 unique(uniqueDoubleSupportPrimaryIndices.begin(),
                        uniqueDoubleSupportPrimaryIndices.end()));
    auto secCompensation =
        distance(uniqueDoubleSupportSecondaryIndices.begin(),
                 unique(uniqueDoubleSupportSecondaryIndices.begin(),
                        uniqueDoubleSupportSecondaryIndices.end()));
    unpairedBreaksSoft -= priCompensation;
    unpairedBreaksHard -= secCompensation;
    pairedBreaksSoft += priCompensation;
    pairedBreaksHard += secCompensation;
}
void
Breakpoint::collectMateSupport() {
    sort(poolLeft.begin(), poolLeft.end());
    sort(poolRight.begin(), poolRight.end());
    compressMatePool(poolLeft);
    compressMatePool(poolRight);
    auto leftDiscordantsTotal = 0, rightDiscordantsTotal = 0;
    for (const auto &mateInfo : poolLeft) {
        leftDiscordantsTotal += mateInfo.matePower;
    }
    for (const auto &mateInfo : poolRight) {
        rightDiscordantsTotal += mateInfo.matePower;
    }
    auto leftSideExpectedErrors = 0.0;
    auto rightSideExpectedErrors = 0.0;
    if (PROPERPAIRCOMPENSATIONMODE) {
        leftSideExpectedErrors =
            leftSideDiscordantCandidates * IMPROPERPAIRRATIO;
        rightSideExpectedErrors =
            rightSideDiscordantCandidates * IMPROPERPAIRRATIO;
        leftDiscordantsTotal =
            max(0, static_cast<int>(
                       round(leftDiscordantsTotal - leftSideExpectedErrors)));
        rightDiscordantsTotal =
            max(0, static_cast<int>(
                       round(rightDiscordantsTotal - rightSideExpectedErrors)));
        leftSideExpectedErrors *= 0.5;
        rightSideExpectedErrors *= 0.5;
    }
    for (auto &sa : doubleSidedMatches) {
        if (sa.isDistant()) {
            if (sa.isEncounteredM()) {
                sa.setExpectedDiscordants(leftDiscordantsTotal);
                collectMateSupportHelper(sa, poolLeft, poolLowQualLeft);
            } else {
                sa.setExpectedDiscordants(rightDiscordantsTotal);
                collectMateSupportHelper(sa, poolRight, poolLowQualRight);
            }
        }
    }
    for (auto &sa : supplementsPrimary) {
        if (sa.isDistant()) {
            if (sa.isEncounteredM()) {
                sa.setExpectedDiscordants(leftDiscordantsTotal);
                collectMateSupportHelper(sa, poolLeft, poolLowQualLeft);
            } else {
                sa.setExpectedDiscordants(rightDiscordantsTotal);
                collectMateSupportHelper(sa, poolRight, poolLowQualRight);
            }
        } else if (sa.getSupport() < BPSUPPORTTHRESHOLD) {
            sa.setToRemove(true);
        }
    }
    for (auto &sa : supplementsPrimary) {
        if (sa.getMapq() == 0 && sa.getMateSupport() == 0 &&
            sa.getDistinctReads() < 2 * BPSUPPORTTHRESHOLD) {
            sa.setToRemove(true);
        }
    }
    vector<SuppAlignment> candidateSupplementsSecondary{};
    sort(supplementsSecondary.begin(), supplementsSecondary.end(),
         [](const SuppAlignment &a, const SuppAlignment &b) {
             return a.getMapq() < b.getMapq();
         });
    while (!supplementsSecondary.empty()) {
        auto foundMatch = false;
        for (auto &sa : candidateSupplementsSecondary) {
            if (sa.saCloseness(supplementsSecondary.back(), 100)) {
                if (sa.isFuzzy() && !supplementsSecondary.back().isFuzzy()) {
                    sa.removeFuzziness(supplementsSecondary.back());
                } else if (sa.isFuzzy() &&
                           supplementsSecondary.back().isFuzzy()) {
                    sa.extendSuppAlignment(
                        supplementsSecondary.back().getPos(),
                        supplementsSecondary.back().getExtendedPos());
                }
                if (supplementsSecondary.back().getMapq() > sa.getMapq()) {
                    sa.setMapq(supplementsSecondary.back().getMapq());
                }
                if (sa.isNullMapqSource() &&
                    !supplementsSecondary.back().isNullMapqSource()) {
                    sa.setNullMapqSource(false);
                }
                for (auto index : supplementsSecondary.back()
                                      .getSupportingIndicesSecondary()) {
                    sa.addSecondarySupportIndices(index);
                }
                foundMatch = true;
                break;
            }
        }
        if (!foundMatch) {
            candidateSupplementsSecondary.push_back(
                supplementsSecondary.back());
        }
        supplementsSecondary.pop_back();
    }
    for (auto &sa : candidateSupplementsSecondary) {
        sa.finalizeSupportingIndices();
    }
    vector<int> originIndices{};
    for (auto &sa : candidateSupplementsSecondary) {
        if (sa.isDistant()) {
            if (sa.isEncounteredM()) {
                sa.setExpectedDiscordants(leftDiscordantsTotal);
                collectMateSupportHelper(sa, poolLeft, poolLowQualLeft);
            } else {
                sa.setExpectedDiscordants(rightDiscordantsTotal);
                collectMateSupportHelper(sa, poolRight, poolLowQualRight);
            }
            if (sa.getMateSupport() > 0) {
                doubleSidedMatches.push_back(sa);
            } else {
                if (sa.isLowMapqSource() || sa.isNullMapqSource() ||
                    sa.getMapq() < 13) {
                    for (auto index : sa.getSupportingIndicesSecondary()) {
                        originIndices.push_back(index);
                    }
                }
            }
        }
    }
    sort(originIndices.begin(), originIndices.end());
    int uniqueCount = unique(originIndices.begin(), originIndices.end()) -
                      originIndices.begin();
    uniqueCount = min(lowQualBreaksHard, uniqueCount);
    lowQualBreaksHard -= uniqueCount;
    lowQualBreaksSoft += uniqueCount;

    for (const auto &mateInfo : poolLeft) {
        if (!mateInfo.saSupporter && mateInfo.evidenceLevel == 3 &&
            mateInfo.matePower / (0.0 + leftDiscordantsTotal) >= 0.33 &&
            (pos - mateInfo.readEndPos) < DEFAULTREADLENGTH / 2) {
            supplementsPrimary.emplace_back(
                mateInfo.mateChrIndex, mateInfo.mateStartPos,
                mateInfo.matePower, leftDiscordantsTotal, true,
                mateInfo.inversionSupport > mateInfo.straightSupport,
                mateInfo.mateEndPos, false, false, false, -1);
        }
    }
    for (const auto &mateInfo : poolRight) {
        if (!mateInfo.saSupporter && mateInfo.evidenceLevel == 3 &&
            mateInfo.matePower / (0.0 + rightDiscordantsTotal) >= 0.33 &&
            (mateInfo.readStartPos - pos) < DEFAULTREADLENGTH / 2) {
            supplementsPrimary.emplace_back(
                mateInfo.mateChrIndex, mateInfo.mateStartPos,
                mateInfo.matePower, rightDiscordantsTotal, false,
                mateInfo.inversionSupport > mateInfo.straightSupport,
                mateInfo.mateEndPos, false, false, false, -1);
        }
    }
    for (auto &sa : doubleSidedMatches) {
        if (sa.isFuzzy() && sa.getMateSupport() == 0) {
            sa.setToRemove(true);
        }
    }
    cleanUpVector(doubleSidedMatches);

    for (auto &sa : supplementsPrimary) {
        if (sa.isFuzzy() && sa.getMateSupport() == 0) {
            sa.setToRemove(true);
        }
    }
    cleanUpVector(supplementsPrimary);
    for (auto &sa : doubleSidedMatches) {
        for (auto &sa2 : supplementsPrimary) {
            if (sa.saCloseness(sa2, 5)) {
                sa.mergeSa(sa2);
                sa2.setToRemove(true);
            }
        }
    }
    cleanUpVector(supplementsPrimary);
    for (auto &sa : doubleSidedMatches) {
        if (!(sa.getSupport() > 0 && sa.getSecondarySupport() > 0) &&
            sa.getSupport() + sa.getSecondarySupport() + sa.getMateSupport() <
                BPSUPPORTTHRESHOLD) {
            sa.setToRemove(true);
        } else {
            if (sa.isDistant() && PROPERPAIRCOMPENSATIONMODE) {
                if (sa.isEncounteredM()) {
                    if (sa.getMateSupport() < leftSideExpectedErrors) {
                        sa.setProperPairErrorProne(true);
                    }
                } else {
                    if (sa.getMateSupport() < rightSideExpectedErrors) {
                        sa.setProperPairErrorProne(true);
                    }
                }
                if (sa.getMateSupport() > sa.getExpectedDiscordants()) {
                    sa.setExpectedDiscordants(sa.getMateSupport());
                }
            }
        }
    }
    cleanUpVector(doubleSidedMatches);
    for (auto &sa : supplementsPrimary) {
        if (!(sa.getSupport() > 0 && sa.getSecondarySupport() > 0) &&
            sa.getSupport() + sa.getSecondarySupport() + sa.getMateSupport() <
                BPSUPPORTTHRESHOLD) {
            sa.setToRemove(true);
        } else {
            if (sa.isDistant() && PROPERPAIRCOMPENSATIONMODE) {
                if (sa.isEncounteredM()) {
                    if (sa.getMateSupport() < leftSideExpectedErrors) {
                        sa.setProperPairErrorProne(true);
                    }
                } else {
                    if (sa.getMateSupport() < rightSideExpectedErrors) {
                        sa.setProperPairErrorProne(true);
                    }
                }
                if (sa.getMateSupport() > sa.getExpectedDiscordants()) {
                    sa.setExpectedDiscordants(sa.getMateSupport());
                }
            }
        }
    }
    cleanUpVector(supplementsPrimary);
}

void
Breakpoint::compressMatePool(vector<MateInfo> &discordantAlignmentsPool) {
    if (discordantAlignmentsPool.empty())
        return;
    auto lastIndex = 0;
    for (auto i = 1; i < static_cast<int>(discordantAlignmentsPool.size());
         ++i) {
        if (discordantAlignmentsPool[lastIndex].mateChrIndex !=
                discordantAlignmentsPool[i].mateChrIndex ||   //
            discordantAlignmentsPool[i].mateStartPos -
                    discordantAlignmentsPool[lastIndex].mateEndPos >
                3.5 * DEFAULTREADLENGTH) {
            lastIndex = i;
        } else {
            discordantAlignmentsPool[lastIndex].mateEndPos =
                max(discordantAlignmentsPool[lastIndex].mateEndPos,
                    discordantAlignmentsPool[i].mateEndPos);
            discordantAlignmentsPool[lastIndex].mateStartPos =
                min(discordantAlignmentsPool[lastIndex].mateStartPos,
                    discordantAlignmentsPool[i].mateStartPos);
            ++discordantAlignmentsPool[lastIndex].matePower;
            if (discordantAlignmentsPool[i].inverted) {
                ++discordantAlignmentsPool[lastIndex].inversionSupport;
            } else {
                ++discordantAlignmentsPool[lastIndex].straightSupport;
            }
            if (abs(pos - discordantAlignmentsPool[i].readStartPos) <=
                abs(pos - discordantAlignmentsPool[i].readEndPos)) {
                // left side
                if (abs(pos - discordantAlignmentsPool[lastIndex].readEndPos) >
                    abs(pos - discordantAlignmentsPool[i].readEndPos)) {
                    discordantAlignmentsPool[lastIndex].readStartPos =
                        discordantAlignmentsPool[i].readStartPos;
                    discordantAlignmentsPool[lastIndex].readEndPos =
                        discordantAlignmentsPool[i].readEndPos;
                }
            } else {
                // right side
                if (abs(pos -
                        discordantAlignmentsPool[lastIndex].readStartPos) >
                    abs(pos - discordantAlignmentsPool[i].readStartPos)) {
                    discordantAlignmentsPool[lastIndex].readStartPos =
                        discordantAlignmentsPool[i].readStartPos;
                    discordantAlignmentsPool[lastIndex].readEndPos =
                        discordantAlignmentsPool[i].readEndPos;
                }
            }
            if ((discordantAlignmentsPool[lastIndex].source == 0 &&
                 discordantAlignmentsPool[i].source == 1) ||
                (discordantAlignmentsPool[lastIndex].source == 1 &&
                 discordantAlignmentsPool[i].source == 0)) {
                discordantAlignmentsPool[lastIndex].evidenceLevel = 2;
                discordantAlignmentsPool[lastIndex].source = 2;
            }
            if (discordantAlignmentsPool[lastIndex].evidenceLevel != 3 &&
                discordantAlignmentsPool[i].evidenceLevel == 3) {
                discordantAlignmentsPool[lastIndex].evidenceLevel = 3;
                discordantAlignmentsPool[lastIndex].source = 2;
            }
            discordantAlignmentsPool[i].toRemove = true;
        }
    }

    for (auto &cluster : discordantAlignmentsPool) {
        if (cluster.evidenceLevel == 1 &&
            cluster.matePower < BPSUPPORTTHRESHOLD) {
            cluster.toRemove = true;
        }
    }
    cleanUpVector(discordantAlignmentsPool);
}

void
Breakpoint::fillMatePool(
    const deque<MateInfo> &discordantAlignmentsPool,
    const deque<MateInfo> &discordantLowQualAlignmentsPool,
    const deque<MateInfo> &discordantAlignmentCandidatesPool) {
    poolLeft.reserve(discordantAlignmentsPool.size());
    poolRight.reserve(discordantAlignmentsPool.size());
    {
        auto i = 0u;
        for (; i < discordantAlignmentsPool.size(); ++i) {
            if (discordantAlignmentsPool[i].readStartPos >= pos) {
                break;
            } else {
                if (discordantAlignmentsPool[i].readEndPos <= pos) {
                    poolLeft.push_back(discordantAlignmentsPool[i]);
                } else {
                    poolLeft.push_back(discordantAlignmentsPool[i]);
                    poolRight.push_back(discordantAlignmentsPool[i]);
                }
            }
        }
        for (; i < discordantAlignmentsPool.size(); ++i) {
            poolRight.push_back(discordantAlignmentsPool[i]);
        }
    }
    if (PROPERPAIRCOMPENSATIONMODE) {
        auto i = 0u;
        for (; i < discordantAlignmentCandidatesPool.size(); ++i) {
            if (discordantAlignmentCandidatesPool[i].readStartPos >= pos) {
                break;
            } else {
                if (discordantAlignmentCandidatesPool[i].readEndPos <= pos) {
                    ++leftSideDiscordantCandidates;
                } else {
                    ++leftSideDiscordantCandidates;
                    ++rightSideDiscordantCandidates;
                }
            }
        }
        for (; i < discordantAlignmentCandidatesPool.size(); ++i) {
            ++rightSideDiscordantCandidates;
        }
    }
    poolLowQualLeft.reserve(discordantLowQualAlignmentsPool.size());
    poolLowQualRight.reserve(discordantLowQualAlignmentsPool.size());
    {
        auto i = 0u;
        for (; i < discordantLowQualAlignmentsPool.size(); ++i) {
            if (discordantLowQualAlignmentsPool[i].readStartPos <
                pos - DISCORDANTLOWQUALLEFTRANGE) {
                continue;
            }
            if (discordantLowQualAlignmentsPool[i].readStartPos >= pos) {
                break;
            } else {
                if (discordantLowQualAlignmentsPool[i].readEndPos <= pos) {
                    poolLowQualLeft.push_back(
                        discordantLowQualAlignmentsPool[i]);
                } else {
                    poolLowQualLeft.push_back(
                        discordantLowQualAlignmentsPool[i]);
                    poolLowQualRight.push_back(
                        discordantLowQualAlignmentsPool[i]);
                }
            }
        }
        for (; i < discordantLowQualAlignmentsPool.size(); ++i) {
            if (discordantLowQualAlignmentsPool[i].readStartPos >
                pos + DISCORDANTLOWQUALRIGHTRANGE) {
                break;
            }
            poolLowQualRight.push_back(discordantLowQualAlignmentsPool[i]);
        }
    }
}

void
Breakpoint::collectMateSupportHelper(
    SuppAlignment &sa, vector<MateInfo> &discordantAlignmentsPool,
    vector<MateInfo> &discordantLowQualAlignmentsPool) {
    auto maxEvidenceLevel = 0;
    for (auto &mateInfo : discordantAlignmentsPool) {
        if (mateInfo.suppAlignmentFuzzyMatch(sa)) {
            if (!mateInfo.saSupporter) {
                mateSupport += mateInfo.matePower;
                mateInfo.saSupporter = true;
            }
            sa.incrementMateSupport(mateInfo.matePower);
            if (sa.isFuzzy()) {
                sa.extendSuppAlignment(mateInfo.mateStartPos,
                                       mateInfo.mateEndPos);
            }
            if (mateInfo.evidenceLevel > maxEvidenceLevel) {
                maxEvidenceLevel = mateInfo.evidenceLevel;
            }
        }
    }
    int lowQualSupports{0};
    for (auto &mateInfo : discordantLowQualAlignmentsPool) {
        if (mateInfo.suppAlignmentFuzzyMatch(sa)) {
            if (!mateInfo.saSupporter) {
                mateSupport += 1;
                mateInfo.saSupporter = true;
            }
            sa.incrementMateSupport(1);
            ++lowQualSupports;
            auto bpPosMatch = false;
            for (const auto bpPos : mateInfo.bpLocs) {
                if (bpPos == pos) {
                    bpPosMatch = true;
                    break;
                }
            }
            if (!bpPosMatch) {
                if (mateInfo.evidenceLevel > maxEvidenceLevel) {
                    maxEvidenceLevel = mateInfo.evidenceLevel;
                }
            }
        }
    }
    auto lowQualDiscordantSupports =
        lowQualSupports +
        min(lowQualSupports,
            static_cast<int>(discordantLowQualAlignmentsPool.size()) -
                lowQualSupports);
    sa.setExpectedDiscordants(sa.getExpectedDiscordants() +
                              lowQualDiscordantSupports);
    if (sa.getMateSupport() == 0) {
        if (sa.getSecondarySupport() == 0 &&
            sa.getSupport() < BPSUPPORTTHRESHOLD) {
            sa.setToRemove(true);
        } else {
            sa.setSuspicious(true);
        }
    } else {
        if (maxEvidenceLevel < 3) {
            sa.setSemiSuspicious(true);
            if (!((0.0 + sa.getMateSupport()) / (sa.getExpectedDiscordants()) >
                  0.33)) {
                if (sa.getSecondarySupport() == 0 &&
                    sa.getSupport() < BPSUPPORTTHRESHOLD) {
                    sa.setToRemove(true);
                }
            }
        } else {
            sa.setSemiSuspicious(false);
        }
    }
}

Breakpoint::Breakpoint(const string &bpIn, bool ignoreOverhang)
    : covFinalized{true}, missingInfoBp{false}, chrIndex{0}, pos{0},
      normalSpans{0}, lowQualSpansSoft{0}, lowQualSpansHard{0},
      unpairedBreaksSoft{0}, unpairedBreaksHard{0}, breaksShortIndel{0},
      lowQualBreaksSoft{0}, lowQualBreaksHard{0}, repetitiveOverhangBreaks{0},
      pairedBreaksSoft{0}, pairedBreaksHard{0}, mateSupport{0}, leftCoverage{0},
      rightCoverage{0}, hitsInMref{0}, germline{false} {
    auto index = 0;
    vector<int> bpChunkPositions{};
    bpChunkPositions.reserve(7);
    for (auto it = bpIn.cbegin(); it != bpIn.cend(); ++it) {
        if (*it == '\t') {
            bpChunkPositions.push_back(index);
        }
        ++index;
    }
    chrIndex = ChrConverter::readChromosomeIndex(bpIn.cbegin(), '\t');

    for (auto i = bpChunkPositions[0] + 1; i < bpChunkPositions[1]; ++i) {
        pos = pos * 10 + (bpIn[i] - '0');
    }
    auto mode = 0;
    for (auto i = bpChunkPositions[2] + 1; i < bpChunkPositions[3]; ++i) {
        if (bpIn[i] == ',') {
            ++mode;
        } else {
            switch (mode) {
            case 0:
                pairedBreaksSoft = 10 * pairedBreaksSoft + (bpIn[i] - '0');
                break;
            case 1:
                pairedBreaksHard = 10 * pairedBreaksHard + (bpIn[i] - '0');
                break;
            case 2:
                mateSupport = 10 * mateSupport + (bpIn[i] - '0');
                break;
            case 3:
                unpairedBreaksSoft = 10 * unpairedBreaksSoft + (bpIn[i] - '0');
                break;
            case 4:
                unpairedBreaksHard = 10 * unpairedBreaksHard + (bpIn[i] - '0');
                break;
            case 5:
                breaksShortIndel = 10 * breaksShortIndel + (bpIn[i] - '0');
                break;
            case 6:
                normalSpans = 10 * normalSpans + (bpIn[i] - '0');
                break;
            case 7:
                lowQualSpansSoft = 10 * lowQualSpansSoft + (bpIn[i] - '0');
                break;
            case 8:
                lowQualSpansHard = 10 * lowQualSpansHard + (bpIn[i] - '0');
                break;
            case 9:
                lowQualBreaksSoft = 10 * lowQualBreaksSoft + (bpIn[i] - '0');
                break;
            case 10:
                lowQualBreaksHard = 10 * lowQualBreaksHard + (bpIn[i] - '0');
                break;
            case 11:
                repetitiveOverhangBreaks =
                    10 * repetitiveOverhangBreaks + (bpIn[i] - '0');
                break;
            default:
                break;
            }
        }
    }
    mode = 0;
    for (auto i = bpChunkPositions[3] + 1; i < bpChunkPositions[4]; ++i) {
        if (bpIn[i] == ',') {
            ++mode;
        } else {
            switch (mode) {
            case 0:
                leftCoverage = 10 * leftCoverage + (bpIn[i] - '0');
                break;
            case 1:
                rightCoverage = 10 * rightCoverage + (bpIn[i] - '0');
                break;
            default:
                break;
            }
        }
    }

    auto shortClipTotal = normalSpans - min(leftCoverage, rightCoverage);
    if (shortClipTotal > 0) {
        normalSpans -= shortClipTotal;
        if (pairedBreaksSoft > 0) {
            pairedBreaksSoft += shortClipTotal;
        } else {
            unpairedBreaksSoft += shortClipTotal;
        }
    }

    if (bpIn[bpChunkPositions[4] + 1] == '#') {
        missingInfoBp = true;
    } else {
        if (bpIn[bpChunkPositions[4] + 1] != '.') {
            string saStr{};
            for (auto i = bpChunkPositions[4] + 1; i < bpChunkPositions[5];
                 ++i) {
                if (bpIn[i] == ';') {
                    doubleSidedMatches.emplace_back(saStr);
                    saStr.clear();
                } else {
                    saStr.push_back(bpIn[i]);
                }
            }
            SuppAlignment saTmp{saStr};
            if (saTmp.getChrIndex() < 1002) {
                doubleSidedMatches.push_back(saTmp);
            }
        }
        if (bpIn[bpChunkPositions[5] + 1] != '.') {
            string saStr{};
            for (auto i = bpChunkPositions[5] + 1; i < bpChunkPositions[6];
                 ++i) {
                if (bpIn[i] == ';') {
                    supplementsPrimary.emplace_back(saStr);
                    saStr.clear();
                } else {
                    saStr.push_back(bpIn[i]);
                }
            }
            SuppAlignment saTmp{saStr};
            if (saTmp.getChrIndex() < 1002) {
                supplementsPrimary.push_back(saTmp);
            }
        }
        cleanUpVector(supplementsPrimary);
        saHomologyClashSolver();
        if (!ignoreOverhang && bpIn[bpChunkPositions[6] + 1] != '.') {
            string overhang{};
            for (auto i = bpChunkPositions[6] + 1;
                 i < static_cast<int>(bpIn.length()); ++i) {
                if (bpIn[i] == ';') {
                    consensusOverhangs.emplace_back(overhang);
                    overhang.clear();
                } else {
                    overhang.push_back(bpIn[i]);
                }
            }
            consensusOverhangs.emplace_back(overhang);
        }
    }
}
void
Breakpoint::saHomologyClashSolver() {
    for (auto i = 0u; i < doubleSidedMatches.size(); ++i) {
        if (!doubleSidedMatches[i].isDistant() ||
            doubleSidedMatches[i].getMateSupport() == 0) {
            continue;
        }
        bool anyMatch{false};
        bool semiSuspiciousRescue{false};
        for (auto j = 0u; j < doubleSidedMatches.size(); ++j) {
            if (j == i) {
                continue;
            }
            if (doubleSidedMatches[i].saDistHomologyRescueCloseness(
                    doubleSidedMatches[j], 200000)) {
                if (!semiSuspiciousRescue &&
                    doubleSidedMatches[i].isSemiSuspicious() &&
                    !doubleSidedMatches[j].isSemiSuspicious()) {
                    semiSuspiciousRescue = true;
                }
                anyMatch = true;
                break;
            }
        }
        if (!anyMatch) {
            for (auto j = 0u; j < supplementsPrimary.size(); ++j) {
                if (doubleSidedMatches[i].saDistHomologyRescueCloseness(
                        supplementsPrimary[j], 200000)) {
                    if (!semiSuspiciousRescue &&
                        doubleSidedMatches[i].isSemiSuspicious() &&
                        !supplementsPrimary[j].isSemiSuspicious()) {
                        semiSuspiciousRescue = true;
                    }
                    anyMatch = true;
                    break;
                }
            }
        }
        if (anyMatch) {
            doubleSidedMatches[i].padMateSupportHomologyRescue();
            if (semiSuspiciousRescue) {
                doubleSidedMatches[i].setSemiSuspicious(false);
            }
        }
    }
    for (auto i = 0u; i < supplementsPrimary.size(); ++i) {
        if (!supplementsPrimary[i].isDistant() ||
            supplementsPrimary[i].getMateSupport() == 0) {
            continue;
        }
        bool anyMatch{false};
        bool semiSuspiciousRescue{false};
        for (auto j = 0u; j < supplementsPrimary.size(); ++j) {
            if (j == i) {
                continue;
            }
            if (supplementsPrimary[i].saDistHomologyRescueCloseness(
                    supplementsPrimary[j], 100000)) {
                if (!semiSuspiciousRescue &&
                    supplementsPrimary[i].isSemiSuspicious() &&
                    !supplementsPrimary[j].isSemiSuspicious()) {
                    semiSuspiciousRescue = true;
                }
                anyMatch = true;
                break;
            }
        }
        if (!anyMatch) {
            for (auto j = 0u; j < doubleSidedMatches.size(); ++j) {
                if (supplementsPrimary[i].saDistHomologyRescueCloseness(
                        doubleSidedMatches[j], 100000)) {
                    if (!semiSuspiciousRescue &&
                        supplementsPrimary[i].isSemiSuspicious() &&
                        !doubleSidedMatches[j].isSemiSuspicious()) {
                        semiSuspiciousRescue = true;
                    }
                    anyMatch = true;
                    break;
                }
            }
        }
        if (anyMatch) {
            supplementsPrimary[i].padMateSupportHomologyRescue();
            if (semiSuspiciousRescue) {
                supplementsPrimary[i].setSemiSuspicious(false);
            }
        }
    }
}

SuppAlignment *
Breakpoint::searchFuzzySa(const SuppAlignment &fuzzySa) {
    SuppAlignment *match = nullptr;
    for (auto &saDouble : doubleSidedMatches) {
        if (saDouble.saCloseness(fuzzySa, 1)) {
            match = &saDouble;
            return match;
        }
    }
    for (auto &saSingle : supplementsPrimary) {
        if (saSingle.saCloseness(fuzzySa, 1)) {
            match = &saSingle;
            return match;
        }
    }
    return nullptr;
}

}   // namespace sophia

/* namespace sophia */
