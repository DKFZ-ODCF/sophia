/*
 * DeFuzzier.cpp
 *
 *  Created on: 27 Nov 2016
 *      Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical Bioinformatics, Bioinformatics and Omics Data Analytics and currently Neuroblastoma Genomics)
 *      Copyright (C) 2018 Umut H. Toprak, Matthias Schlesner, Roland Eils and DKFZ Heidelberg
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

#include <DeFuzzier.h>
#include <algorithm>

namespace sophia {
    
    DeFuzzier::DeFuzzier(ChrSize maxDistanceIn, bool mrefModeIn)
        : MAX_DISTANCE { maxDistanceIn },
          MREF_MODE { mrefModeIn } {}

    void DeFuzzier::deFuzzyDb(std::vector<BreakpointReduced>& bps) const {
        for (auto it = bps.begin(); it != bps.end(); ++it) {
            for (auto &sa : it->getSupplementsPtr()) {
                if (sa->isFuzzy()) {
                    auto saTmp = sa;
                    processFuzzySa(bps, it, saTmp);
                }
            }
            it->removeMarkedFuzzies();
        }
        for (auto &bp : bps) {
            bp.removeMarkedFuzzies();
        }
    }

    void DeFuzzier::processFuzzySa(std::vector<BreakpointReduced>& bps,
                                   std::vector<BreakpointReduced>::iterator startingIt,
                                   SuppAlignmentAnno* startingSa) const {
        auto consensusSa = startingSa;
        std::vector<SuppAlignmentAnno*> processedSas { startingSa };
        if (!startingSa->isEncounteredM()) {
            dbSweep(bps, startingIt, 1, consensusSa, processedSas);
            dbSweep(bps, startingIt, -1, consensusSa, processedSas);
        } else {
            dbSweep(bps, startingIt, -1, consensusSa, processedSas);
            dbSweep(bps, startingIt, 1, consensusSa, processedSas);
        }
        selectBestSa(processedSas, consensusSa);
    }

    void DeFuzzier::dbSweep(std::vector<BreakpointReduced>& bps,
                            std::vector<BreakpointReduced>::iterator startingIt,
                            int increment,
                            SuppAlignmentAnno* consensusSa,
                            std::vector<SuppAlignmentAnno*>& processedSas) const {
        auto it = startingIt;
        if (it == bps.begin() || it == bps.end()) {
            return;
        }
        advance(it, increment);
        while (it != bps.begin() && it != bps.end()) {
            auto res = it->searchFuzzySa(*consensusSa);
            if (!res && ChrSize(abs(static_cast<long>(startingIt->getPos()) - static_cast<long>(it->getPos()))) > MAX_DISTANCE) {
                break;
            } else {
                if (res) {
                    processedSas.push_back(res);
                    if (res->isFuzzy()) {
                        consensusSa->extendSuppAlignment(std::min(res->getPos(), consensusSa->getPos()),
                                                         std::max(res->getExtendedPos(), consensusSa->getExtendedPos()));
                    }
                }
            }
            advance(it, increment);
        }
    }

    void DeFuzzier::selectBestSa(std::vector<SuppAlignmentAnno*>& processedSas,
                                 SuppAlignmentAnno* consensusSa) const {
        auto maxMateScore = -1;
        auto maxExpectedDiscordants = -1;
        auto index = 0;
        std::vector<int> nonFuzzyIndices { };
        for (auto &sa : processedSas) {
            if (sa->getMateSupport() > maxMateScore) {
                maxMateScore = sa->getMateSupport();
            }
            if (sa->getExpectedDiscordants() > maxExpectedDiscordants) {
                maxExpectedDiscordants = sa->getExpectedDiscordants();
            }
            if (!sa->isFuzzy()) {
                nonFuzzyIndices.push_back(index);
            }
            sa->setToRemove(true);
            ++index;
        }
        SuppAlignmentAnno* selectedSa = nullptr;
        if (!nonFuzzyIndices.empty()) {
            auto bestElement = max_element(
                nonFuzzyIndices.begin(),
                nonFuzzyIndices.end(), //
                [&](unsigned int a, unsigned int b) {
                    return (processedSas[a]->getSupport() + processedSas[a]->getSecondarySupport())
                                < (processedSas[b]->getSupport() + processedSas[b]->getSecondarySupport());
                });
            selectedSa = processedSas[(size_t) *bestElement];
        } else {
            auto bestElement = max_element(processedSas.begin(), processedSas.end(),
                    [&](SuppAlignmentAnno* a, SuppAlignmentAnno* b) {return a->getMateSupport() < b->getMateSupport();});
            selectedSa = *bestElement;
            selectedSa->extendSuppAlignment(consensusSa->getPos(), consensusSa->getExtendedPos());
        }
        selectedSa->setToRemove(false);
        selectedSa->setMateSupport(maxMateScore);
        selectedSa->setExpectedDiscordants(maxExpectedDiscordants);
    }

    void DeFuzzier::deFuzzyDb(std::vector<MrefEntry>& bps) const {
        for (auto it = bps.begin(); it != bps.end(); ++it) {
            if (!it->isValid()) {
                continue;
            }
            for (auto &sa : it->getSupplementsPtr()) {
                if (sa->isToRemove()) {
                    continue;
                }
                if (sa->isFuzzy() || sa->isStrictFuzzy()) {
                    auto saTmp = sa;
                    processFuzzySa(bps, it, saTmp);
                }
            }
            it->removeMarkedFuzzies();
            if (it->getValidityScore() == 0 && it->getSuppAlignments().empty()) {
                it->setAsInvalid();
            }
        }
        for (MrefEntry &bp : bps) {
            if (bp.getValidityScore() == -1 || !bp.isValid()) {
                bp.setAsInvalid();
            }
        }
    }

    void DeFuzzier::processFuzzySa(std::vector<MrefEntry>& bps,
                                   std::vector<MrefEntry>::iterator startingIt,
                                   SuppAlignmentAnno* startingSa) const {
        SuppAlignmentAnno* consensusSa = startingSa;
        std::unordered_set<unsigned short> fileIndices { };
        for (auto fileIndex : startingIt->getFileIndices()) {
            fileIndices.insert(fileIndex);
        }
        std::vector<SuppAlignmentAnno*> processedSas { startingSa };
        if (!startingSa->isEncounteredM()) {
            dbSweep(bps, startingIt, fileIndices, 1, consensusSa, processedSas);
            dbSweep(bps, startingIt, fileIndices, -1, consensusSa, processedSas);
        } else {
            dbSweep(bps, startingIt, fileIndices, -1, consensusSa, processedSas);
            dbSweep(bps, startingIt, fileIndices, 1, consensusSa, processedSas);
        }
        selectBestSa(processedSas, consensusSa, fileIndices);
    }

    void DeFuzzier::dbSweep(std::vector<MrefEntry>& bps,
                            std::vector<MrefEntry>::iterator startingIt,
                            std::unordered_set<unsigned short>& fileIndices,
                            int increment,
                            SuppAlignmentAnno* consensusSa,
                            std::vector<SuppAlignmentAnno*>& processedSas) const {
        auto it = startingIt;
        if (it == bps.begin() || it == bps.end()) {
            return;
        }
        advance(it, increment);
        while (it != bps.begin() && it != bps.end()) {
            if (it->isValid()) {
                auto res = it->searchFuzzySa(*consensusSa);
                if (!res && abs(static_cast<int>(startingIt->getPos()) - static_cast<int>(it->getPos())) > static_cast<int>(MAX_DISTANCE)) {
                    break;
                } else {
                    if (res) {
                        for (auto fileIndex : it->getFileIndices()) {
                            fileIndices.insert(fileIndex);
                        }
                        processedSas.push_back(res);
                        if (res->isFuzzy()) {
                            consensusSa->extendSuppAlignment(std::min(res->getPos(), consensusSa->getPos()),
                                                             std::max(res->getExtendedPos(), consensusSa->getExtendedPos()));
                        }
                    }
                }
            }
            advance(it, increment);
        }
    }

    void DeFuzzier::selectBestSa(std::vector<SuppAlignmentAnno*>& processedSas,
                                 SuppAlignmentAnno* consensusSa,
                                 const std::unordered_set<unsigned short>& fileIndices) const {
        auto maxMateScore = -1;
        auto maxExpectedDiscordants = -1;
        auto index = 0;
        std::vector<int> nonFuzzyIndices { };
        for (auto &sa : processedSas) {
            if (sa->getMateSupport() > maxMateScore) {
                maxMateScore = sa->getMateSupport();
            }
            if (sa->getExpectedDiscordants() > maxExpectedDiscordants) {
                maxExpectedDiscordants = sa->getExpectedDiscordants();
            }
            if (!sa->isFuzzy()) {
                nonFuzzyIndices.push_back(index);
            }
            sa->setToRemove(true);
            ++index;
        }
        SuppAlignmentAnno* selectedSa = nullptr;
        if (!nonFuzzyIndices.empty()) {
            auto bestElement = max_element(
                nonFuzzyIndices.begin(),
                nonFuzzyIndices.end(),
                [&](unsigned int a, unsigned int b) {
                    return processedSas[a]->getSupportingIndices().size() < processedSas[b]->getSupportingIndices().size();
                });
            selectedSa = processedSas[static_cast<unsigned long>(*bestElement)];
        } else {
            auto bestElement = max_element(processedSas.begin(), processedSas.end(), //
                    [&](SuppAlignmentAnno* a, SuppAlignmentAnno* b) {return a->getSupportingIndices().size() < b->getSupportingIndices().size();});
            if ((*bestElement)->getSupport() + (*bestElement)->getSecondarySupport() == 0) {
                if (consensusSa->isEncounteredM()) {
                    selectedSa = processedSas.back();
                } else {
                    selectedSa = processedSas[0];
                }
            } else {
                selectedSa = *bestElement;
            }
            selectedSa->extendSuppAlignment(consensusSa->getPos(), consensusSa->getExtendedPos());
        }
        selectedSa->setToRemove(false);
        selectedSa->setMateSupport(maxMateScore);
        selectedSa->setExpectedDiscordants(maxExpectedDiscordants);
        selectedSa->mrefSaConsensus(fileIndices);
    }

}

