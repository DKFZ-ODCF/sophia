/*
 * MrefEntry.cpp
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

#include "global.h"
#include "Breakpoint.h"
#include "strtk-wrap.h"
#include "GlobalAppConfig.h"
#include <boost/algorithm/string/join.hpp>
#include <MrefEntry.h>
#include <unordered_set>
#include "BreakpointReduced.h"

namespace sophia {

    boost::format MrefEntry::doubleFormatter { "%.5f" };

    unsigned int MrefEntry::NUM_PIDS { };

    ChrSize MrefEntry::DEFAULT_READ_LENGTH { };

    MrefEntry::MrefEntry() :
                    validity { -1 },
                    pos { std::numeric_limits<ChrSize>::max() },
                    fileIndices { },
                    fileIndicesWithArtifactRatios { },
                    artifactRatios { },
                    suppAlignments { } {

    }

    void MrefEntry::addEntry(BreakpointReduced& tmpBreakpoint,
                             int fileIndex) {
        pos = tmpBreakpoint.getPos();
        auto artifactBreakTotal =
            tmpBreakpoint.getLowQualBreaksSoft() +
            tmpBreakpoint.getLowQualBreaksHard() +
            tmpBreakpoint.getRepetitiveOverhangBreaks();
        auto eventTotal =
            tmpBreakpoint.getPairedBreaksSoft() +
            tmpBreakpoint.getPairedBreaksHard() +
            tmpBreakpoint.getUnpairedBreaksSoft() +
            tmpBreakpoint.getUnpairedBreaksHard() +
            tmpBreakpoint.getBreaksShortIndel();
        auto breakTotal =
            eventTotal +
            artifactBreakTotal;

        if (breakTotal < 200) {
            const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
            for (SuppAlignmentAnno *saPtr : tmpBreakpoint.getSupplementsPtr()) {
                // Original code from bitbucket repository
                // if (saPtr->isSuspicious()
                //     || saPtr->isToRemove()
                //     || (saPtr->getChrIndex() != 1001     // i.e. !mitochondrial
                //         && ChrConverter::indexConverter[saPtr->getChrIndex()] < 0))  // i.e. !compressedMref
                // `indexConverter` is now `indexToCompressedMrefIndex` and a mapping from ChrIndex to
                // CompressedMrefIndex, that contained `-2` values (now `NA` constant) for chromosomes that were not
                // in the compressed mref set.
                // Note that the mitochondrial chromosome itself is *not* among the compressed mrefs. So the
                // condition is somewhat redundant.
                if (saPtr->isSuspicious()
                    || saPtr->isToRemove()
                    || (!chrConverter.isExtrachromosomal(saPtr->getChrIndex())
                        && !chrConverter.isCompressedMref(saPtr->getChrIndex()))
                    ) {
                    continue;
                }
                auto qualCheck = false;
                auto splitTotal = saPtr->getSupport() + saPtr->getSecondarySupport();
                if (saPtr->isDistant()) {
                    auto clonalityCondition = (((0.0 + saPtr->getMateSupport()) / saPtr->getExpectedDiscordants()) >= 0.5);
                    if (!clonalityCondition) {
                        continue;
                    }
                    qualCheck = (splitTotal > 4 && saPtr->getMateSupport() > 2);
                    if (!qualCheck) {
                        qualCheck = (saPtr->getMateSupport() > 4);
                    }
                } else {
                    qualCheck = (splitTotal > 4) || (splitTotal > 2 && saPtr->getSupport() > 0 && saPtr->getSecondarySupport() > 0);
                }
                if (qualCheck) {
                    if (!saMatcher(saPtr)) {
                        auto saTmp = *saPtr;
                        saTmp.mrefSaTransform(fileIndex);
                        suppAlignments.push_back(saTmp);
                    }
                }
            }
        }
        auto covValidity = (tmpBreakpoint.getBreaksShortIndel() > 2 || breakTotal > 9);
        if (!covValidity) {
            if (breakTotal > 4) {
                auto clonality = ((breakTotal + 0.0) / (tmpBreakpoint.getNormalSpans() + breakTotal));
                if (clonality > 0.3) {
                    covValidity = true;
                } else if (clonality > 0.1) {
                    if (tmpBreakpoint.hasOverhang) {
                        covValidity = true;
                    }
                }
            }
        }
        if (eventTotal + artifactBreakTotal > 0) {
            if (covValidity) {
                auto eventTotalStrict = tmpBreakpoint.getPairedBreaksSoft() + tmpBreakpoint.getUnpairedBreaksSoft() + tmpBreakpoint.getPairedBreaksHard();
                auto artifactTotalRelaxed = tmpBreakpoint.getLowQualBreaksSoft() + tmpBreakpoint.getLowQualSpansSoft() + tmpBreakpoint.getRepetitiveOverhangBreaks();
                if ((eventTotalStrict + artifactTotalRelaxed) > 0) {
                    artifactRatios.push_back((0.0 + artifactTotalRelaxed) / (eventTotalStrict + artifactTotalRelaxed));
                    fileIndicesWithArtifactRatios.push_back(fileIndex);
                }
            }
        }
        if (covValidity) {
            fileIndices.push_back(fileIndex);
            validity = 1;
        } else if (!suppAlignments.empty()) {
            fileIndices.push_back(fileIndex);
            validity = 0;
        }
    }

    void MrefEntry::mergeMrefEntries(MrefEntry& entry2) {
        pos = entry2.getPos();
        for (auto artifactRatio : entry2.getArtifactRatios()) {
            artifactRatios.push_back(artifactRatio);
        }
        for (auto fileIndex : entry2.getFileIndicesWithArtifactRatios()) {
            fileIndicesWithArtifactRatios.push_back(fileIndex);
        }
        for (auto fileIndex : entry2.getFileIndices()) {
            fileIndices.push_back(fileIndex);
        }
        for (auto saPtr : entry2.getSupplementsPtr()) {
            if (!saMatcher(saPtr)) {
                suppAlignments.push_back(*saPtr);
            }
        }
        validity = std::max(validity, entry2.getValidityScore());
    }

    /** This prints the output of the `sophiaMref` tool. */
    std::string MrefEntry::printBpInfo(const std::string& chromosome) {
        finalizeFileIndices();
        std::vector<std::string> outputFields { };
        outputFields.emplace_back(chromosome);
        outputFields.emplace_back(strtk::type_to_string<int>(pos));
        outputFields.emplace_back(strtk::type_to_string<int>(pos + 1));
        outputFields.emplace_back(strtk::type_to_string<int>(fileIndices.size()));
        outputFields.emplace_back(strtk::type_to_string<int>(fileIndicesWithArtifactRatios.size()));
        outputFields.emplace_back(boost::str(doubleFormatter % ((fileIndices.size() + 0.0) / NUM_PIDS)));
        outputFields.emplace_back(boost::str(doubleFormatter % ((fileIndicesWithArtifactRatios.size() + 0.0) / NUM_PIDS)));
        if (!artifactRatios.empty()) {
            outputFields.emplace_back(boost::str(doubleFormatter % (accumulate(artifactRatios.cbegin(), artifactRatios.cend(), 0.0) / artifactRatios.size())));
        } else {
            outputFields.emplace_back("NA");
        }
        if (suppAlignments.empty()) {
            outputFields.emplace_back(".");
        } else {
            std::vector<std::string> saFields { };
            saFields.reserve(suppAlignments.size());
            for (auto &sa : suppAlignments) {
                sa.finalizeSupportingIndices();
                if (suppAlignments.size() < 11 || sa.getSupport() >= 0.2 * fileIndices.size()) {
                    saFields.emplace_back(sa.print());
                }
            }
            if (saFields.empty()) {
                outputFields.emplace_back(".");
            } else {
                outputFields.emplace_back(boost::join(saFields, ";"));
            }
        }
        std::vector<std::string> fileIndicesStr { };
        transform(fileIndices.begin(), fileIndices.end(), back_inserter(fileIndicesStr),
                  [](int fileIndex) {return strtk::type_to_string<int>(fileIndex);});
        outputFields.emplace_back(boost::join(fileIndicesStr, ","));
        return boost::join(outputFields, "\t").append("\n");
    }

    // Currently, not used.
    std::string MrefEntry::printArtifactRatios(const std::string& chromosome) {
        std::vector<std::string> outputFields { };
        outputFields.reserve(NUM_PIDS + 3u);
        outputFields.emplace_back(chromosome);
        outputFields.emplace_back(strtk::type_to_string<int>(pos));
        outputFields.emplace_back(strtk::type_to_string<int>(pos + 1));
        std::vector<std::string> artifactRatiosOutput(NUM_PIDS, ".");
        for (size_t i = 0; i < fileIndicesWithArtifactRatios.size(); ++i) {
            artifactRatiosOutput[fileIndicesWithArtifactRatios[i]] =
                boost::str(doubleFormatter % artifactRatios[i]);
        }
        for (const auto &artifactRatio : artifactRatiosOutput) {
            outputFields.push_back(artifactRatio);
        }
        return boost::join(outputFields, "\t").append("\n");
    }

    SuppAlignmentAnno* MrefEntry::searchFuzzySa(const SuppAlignmentAnno& fuzzySa) {
        SuppAlignmentAnno* match = nullptr;
        for (auto &sa : suppAlignments) {
            if (sa.isToRemove()) {
                continue;
            }
            if (sa.saClosenessDirectional(fuzzySa, 1)) {
                match = &sa;
                return match;
            }
        }
        return nullptr;
    }

    bool MrefEntry::saMatcher(SuppAlignmentAnno* saPtr) {
        if (saPtr->isToRemove()
               || saPtr->isSuspicious()
               || (saPtr->getExpectedDiscordants() > 0
                   && !(saPtr->getMateSupport() / (saPtr->getExpectedDiscordants() + 0.0) > 0.1))) {
            return true;
        }
        for (auto &sa : suppAlignments) {
            if (sa.saCloseness(*saPtr, 100)) {
                if (sa.isFuzzy()) {
                    if (saPtr->isFuzzy()) {
                        sa.extendSuppAlignment(saPtr->getPos(), saPtr->getExtendedPos());
                    } else {
                        sa.removeFuzziness(*saPtr);
                    }
                }
                sa.addSupportingIndices(saPtr->getSupportingIndices());
                sa.mergeMrefSa(*saPtr);
                saPtr->setToRemove(true);
                return true;
            }
        }
        return false;
    }

    void MrefEntry::finalizeFileIndices() {
        for (const auto &sa : suppAlignments) {
            auto tmpIndices = sa.getSupportingIndices();
            copy(tmpIndices.begin(), tmpIndices.end(), inserter(fileIndices, fileIndices.end()));
        }
        sort(fileIndices.begin(), fileIndices.end());
        auto last = unique(fileIndices.begin(), fileIndices.end());
        fileIndices.erase(last, fileIndices.end());
    }

} /* namespace sophia */

