/*
 * AnnotationProcessor.cpp
 *
 * Created on: 28 Apr 2016
 * Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical
 * Bioinformatics, Bioinformatics and Omics Data Analytics and currently
 * Neuroblastoma Genomics) Copyright (C) 2018 Umut H. Toprak, Matthias
 * Schlesner, Roland Eils and DKFZ Heidelberg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * LICENSE: GPL
 */

#include "Breakpoint.h"
#include "HelperFunctions.h"
#include "SuppAlignment.h"
#include "GlobalAppConfig.h"
#include <AnnotationProcessor.h>
#include <DeFuzzier.h>
#include <algorithm>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

namespace sophia {

    bool AnnotationProcessor::ABRIDGED_OUTPUT{false};

    AnnotationProcessor::AnnotationProcessor(const std::string &tumorResultsIn,
                                             std::vector<std::vector<MrefEntryAnno>> &mref,
                                             ChrSize defaultReadLengthTumorIn,
                                             bool controlCheckMode,
                                             int GERMLINE_DB_LIMIT)
        : NO_CONTROL_MODE{true},
          GERMLINE_DB_LIMIT{GERMLINE_DB_LIMIT},
          contaminationObserved{false},
          massiveInvFilteringLevel{0},
          filteredResults{},
          visitedLineIndices{} {

        CompressedMrefIndex nCompressedMrefChromosomes = GlobalAppConfig::getInstance().
            getChrConverter().nChromosomesCompressedMref();
        std::vector<std::vector<BreakpointReduced>>::size_type vectorSize =
            std::vector<std::vector<BreakpointReduced>>::size_type(nCompressedMrefChromosomes);

        tumorResults = std::vector<std::vector<BreakpointReduced>>
            { vectorSize, std::vector<BreakpointReduced>{} };
        controlResults = std::vector<std::vector<BreakpointReduced>>
            { vectorSize, std::vector<BreakpointReduced>{} };

        std::unique_ptr<std::ifstream> tumorInputHandle{
            std::make_unique<std::ifstream>(tumorResultsIn, std::ios_base::in | std::ios_base::binary)};
        std::unique_ptr<boost::iostreams::filtering_istream> tumorGzHandle{
            std::make_unique<boost::iostreams::filtering_istream>()};
        tumorGzHandle->push(boost::iostreams::gzip_decompressor());
        tumorGzHandle->push(*tumorInputHandle);
        std::string line;
        int lineIndex = 0;
        const ChrConverter& chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        while (error_terminating_getline(*tumorGzHandle, line)) {
            if (line.front() == '#') {
                continue;
            };

            Breakpoint tmpBp = Breakpoint::parse(line, true);
            CompressedMrefIndex compressedMrefChrIndex;
            if (!chrConverter.isCompressedMref(tmpBp.getChrIndex())) {
                continue;
            } else {
                compressedMrefChrIndex = chrConverter.indexToCompressedMrefIndex(tmpBp.getChrIndex());
            }
            auto hasOverhang = line.back() != '.' && line.back() != '#';
            tumorResults[static_cast<unsigned int>(compressedMrefChrIndex)].
                emplace_back(BreakpointReduced(tmpBp, lineIndex, hasOverhang));
            if (hasOverhang) {
                std::string overhang{};
                for (auto it = line.rbegin(); it != line.rend(); ++it) {
                    if (*it == '\t') {
                        break;
                    } else {
                        overhang.push_back(*it);
                    }
                }
                reverse(overhang.begin(), overhang.end());
                overhangs.emplace_back(lineIndex, overhang);
            } else {
                visitedLineIndices.push_back(lineIndex);
            }
            ++lineIndex;
        }
        for (auto &tres : tumorResults) {
            DeFuzzier deFuzzier { defaultReadLengthTumorIn * 6, false };
            deFuzzier.deFuzzyDb(tres);
        }
        searchMatches(mref);
        if (applyMassiveInversionFiltering(false, controlCheckMode)) {
            ++massiveInvFilteringLevel;
        }
        if (massiveInvFilteringLevel == 1) {
            if (applyMassiveInversionFiltering(true, controlCheckMode)) {
                ++massiveInvFilteringLevel;
            }
        }
        contaminationObserved = applyPathogenContaminationFiltering();
        if (!controlCheckMode && !contaminationObserved &&
            massiveInvFilteringLevel == 0) {
            printUnresolvedRareOverhangs(mref);
        }
    }

    AnnotationProcessor::AnnotationProcessor(
            const std::string &tumorResultsIn,
            std::vector<std::vector<MrefEntryAnno>> &mref,
            const std::string &controlResultsIn,
            ChrSize defaultReadLengthTumorIn,
            ChrSize defaultReadLengthControlIn,
            int germlineDbLimit,
            int lowQualControlIn,
            bool pathogenInControlIn)
        : NO_CONTROL_MODE{false},
          GERMLINE_DB_LIMIT{germlineDbLimit},
          contaminationObserved{false},
          massiveInvFilteringLevel{0},
          filteredResults{} {

        CompressedMrefIndex nCompressedMrefChromosomes = GlobalAppConfig::getInstance().
            getChrConverter().nChromosomesCompressedMref();
        std::vector<std::vector<BreakpointReduced>>::size_type vectorSize =
            std::vector<std::vector<BreakpointReduced>>::size_type(nCompressedMrefChromosomes);

        tumorResults = std::vector<std::vector<BreakpointReduced>>
            { vectorSize, std::vector<BreakpointReduced>{} };
        controlResults = std::vector<std::vector<BreakpointReduced>>
            { vectorSize, std::vector<BreakpointReduced>{} };

        std::unique_ptr<std::ifstream> controlInputHandle{
            std::make_unique<std::ifstream>(controlResultsIn, std::ios_base::in | std::ios_base::binary)};
        std::unique_ptr<boost::iostreams::filtering_istream> controlGzHandle{
            std::make_unique<boost::iostreams::filtering_istream>()};

        controlGzHandle->push(boost::iostreams::gzip_decompressor());
        controlGzHandle->push(*controlInputHandle);
        std::string line;
        auto lineIndex = 0;
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        while (error_terminating_getline(*controlGzHandle, line)) {
            if (line.front() == '#') {
                continue;
            };
            Breakpoint tmpBpPre = Breakpoint::parse(line, true);
            BreakpointReduced tmpBp{tmpBpPre, lineIndex, false};
            if (chrConverter.isTechnical(tmpBp.getChrIndex())) {
                continue;
            }
            if (pathogenInControlIn) {
                if ((tmpBp.getPairedBreaksSoft() + tmpBp.getUnpairedBreaksSoft()) >
                        19 &&
                    (tmpBp.getPairedBreaksHard() + tmpBp.getUnpairedBreaksHard() <
                     3)) {
                    if (line.back() != '.' && line.back() != '#') {
                        std::string overhang{};
                        ChrSize overhangLength = 0;
                        ChrSize maxOverhangLength = 0;
                        for (auto it = line.rbegin(); *it != '\t'; ++it) {
                            switch (*it) {
                            case '(':
                                overhangLength = 0;
                                break;
                            case ':':
                                maxOverhangLength = std::max(maxOverhangLength, overhangLength);
                                overhangLength = 0;
                                break;
                            default:
                                ++overhangLength;
                                break;
                            }
                        }
                        auto maxOverhangLengthRatio =
                            (maxOverhangLength + 0.0) / defaultReadLengthControlIn;
                        if (maxOverhangLengthRatio > 0.7) {
                            continue;
                        }
                    }
                }
            }
            if (lowQualControlIn > 0) {
                SuppAlignmentAnno *bestSa = nullptr;
                auto bestSaSupport = 0;
                if (!tmpBp.getSuppAlignments().empty()) {
                    for (const auto sa : tmpBp.getSupplementsPtr()) {
                        if (!sa->isSuspicious()) {
                            if (bestSa) {
                                auto saSupport = sa->getSupport() +
                                                 sa->getSecondarySupport() +
                                                 sa->getMateSupport();
                                if (saSupport > bestSaSupport) {
                                    bestSa = sa;
                                    bestSaSupport = saSupport;
                                }
                            } else {
                                bestSa = sa;
                            }
                        }
                    }
                }
                auto clipTotal =
                    tmpBp.getPairedBreaksHard() + tmpBp.getPairedBreaksSoft() +
                    tmpBp.getUnpairedBreaksHard() + tmpBp.getUnpairedBreaksSoft();
                if (!bestSa) {
                    if (clipTotal < 10) {
                        continue;
                    }
                } else {
                    if (bestSa->isInverted()) {
                        ChrDistance suppDist = tmpBp.distanceToSupp(*bestSa);
                        if (suppDist > 0) {
                            if (suppDist < 10000) {
                                if (lowQualControlIn == 2) {
                                    continue;
                                } else {
                                    if (bestSa->getSupport() < 5 ||
                                        bestSa->getSecondarySupport() < 5) {
                                        continue;
                                    }
                                }
                            } else {
                                if (bestSa->getSupport() < 5 &&
                                    bestSa->getSecondarySupport() < 5) {
                                    continue;
                                }
                            }
                        } else if (suppDist < 0) {
                            if (bestSa->getSupport() < 5 &&
                                bestSa->getSecondarySupport() < 5) {
                                continue;
                            }
                        }
                    }
                }
            }
            controlResults[static_cast<unsigned int>(chrConverter.indexToCompressedMrefIndex(tmpBp.getChrIndex()))]
                .push_back(tmpBp);
            ++lineIndex;
        }
        for (auto &cres : controlResults) {
            DeFuzzier deFuzzierControl{defaultReadLengthControlIn * 6, false};
            deFuzzierControl.deFuzzyDb(cres);
        }
        std::unique_ptr<std::ifstream> tumorInputHandle{
            std::make_unique<std::ifstream>(tumorResultsIn, std::ios_base::in | std::ios_base::binary)};
        std::unique_ptr<boost::iostreams::filtering_istream> tumorGzHandle{
            std::make_unique<boost::iostreams::filtering_istream>()};

        tumorGzHandle->push(boost::iostreams::gzip_decompressor());
        tumorGzHandle->push(*tumorInputHandle);
        lineIndex = 0;
        while (error_terminating_getline(*tumorGzHandle, line)) {
            if (line.front() == '#') {
                continue;
            };
            Breakpoint tmpBp = Breakpoint::parse(line, true);
            if (!chrConverter.isCompressedMref(tmpBp.getChrIndex())) {
                continue;
            }
            CompressedMrefIndex compressedMrefChrIndex =
                chrConverter.indexToCompressedMrefIndex(tmpBp.getChrIndex());
            auto hasOverhang = line.back() != '.' && line.back() != '#';
            tumorResults[static_cast<unsigned int>(compressedMrefChrIndex)].
                emplace_back(tmpBp, lineIndex, hasOverhang);
            if (line.back() != '.' && line.back() != '#') {
                std::string overhang{};
                for (auto it = line.rbegin(); it != line.rend(); ++it) {
                    if (*it == '\t') {
                        break;
                    } else {
                        overhang.push_back(*it);
                    }
                }
                reverse(overhang.begin(), overhang.end());
                overhangs.emplace_back(lineIndex, overhang);
            } else {
                visitedLineIndices.push_back(lineIndex);
            }
            ++lineIndex;
        }
        for (auto &tres : tumorResults) {
            DeFuzzier deFuzzierTumor{defaultReadLengthTumorIn * 6, false};
            deFuzzierTumor.deFuzzyDb(tres);
        }
        searchMatches(mref);
        if (applyMassiveInversionFiltering(false, false)) {
            ++massiveInvFilteringLevel;
        }
        if (massiveInvFilteringLevel == 1) {
            if (applyMassiveInversionFiltering(true, false)) {
                ++massiveInvFilteringLevel;
            }
        }
        contaminationObserved = applyPathogenContaminationFiltering();
        if (!contaminationObserved && massiveInvFilteringLevel == 0) {
            printUnresolvedRareOverhangs(mref);
        }
    }

    void
    AnnotationProcessor::searchMatches(std::vector<std::vector<MrefEntryAnno>> &mref) {
        CompressedMrefIndex nCompressedMrefChromosomes = GlobalAppConfig::getInstance().
            getChrConverter().nChromosomesCompressedMref();
        for (CompressedMrefIndex mrefIdx = 0; mrefIdx < nCompressedMrefChromosomes; ++mrefIdx) {
            for (size_t dbIdx = 0; dbIdx < tumorResults[static_cast<unsigned int>(mrefIdx)].size(); ++dbIdx) {
                for (const auto &sa : tumorResults[static_cast<unsigned int>(mrefIdx)][dbIdx].getSuppAlignments()) {
                    if (SvEvent::DEBUG_MODE || !sa.isSuspicious()) {
                        if (sa.getSecondarySupport() > 0 ||
                            (sa.getSupport() > 0 && sa.getMateSupport() > 0)) {
                            searchSa(mrefIdx, dbIdx, sa, true, mref);
                        } else {
                            searchSa(mrefIdx, dbIdx, sa, false, mref);
                        }
                    }
                }
            }
        }
    }

    void
    AnnotationProcessor::searchSa(CompressedMrefIndex compressedMrefIndex,
                                  size_t dbIndex,
                                  const SuppAlignmentAnno &sa,
                                  bool doubleSupportSa,
                                  std::vector<std::vector<MrefEntryAnno>> &mref) {
        if (sa.getSupport() + sa.getSecondarySupport() + sa.getMateSupport() < 3) {
            return;
        }
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        if (chrConverter.isExtrachromosomal(sa.getChrIndex())) {
            if (createUnknownMatchSvPreCheck(sa, doubleSupportSa)) {
                createUnknownMatchSv(tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex], sa, mref,
                                     doubleSupportSa);
            }
            return;
        }
        if (!chrConverter.isCompressedMref(sa.getChrIndex())) {
            return;
        }
        CompressedMrefIndex saChrIndex = chrConverter.indexToCompressedMrefIndex(sa.getChrIndex());
        ChrDistance fuzziness = ChrDistance(3 * SuppAlignmentAnno::DEFAULT_READ_LENGTH);
        std::vector<std::pair<int, std::vector<BreakpointReduced>::iterator>> dbHits{};

        if (!tumorResults[static_cast<unsigned int>(saChrIndex)].empty()) {
            auto itStart = lower_bound(tumorResults[static_cast<unsigned int>(saChrIndex)].begin(),
                                       tumorResults[static_cast<unsigned int>(saChrIndex)].end(), sa);
            if (itStart == tumorResults[static_cast<unsigned int>(saChrIndex)].end()) {
                --itStart;
            }
            if (itStart != tumorResults[static_cast<unsigned int>(saChrIndex)].begin() &&
                !itStart->closeToSupp(sa, fuzziness)) {
                --itStart;
                if (!itStart->closeToSupp(sa, fuzziness)) {
                    if (createUnknownMatchSvPreCheck(sa, doubleSupportSa)) {

                        createUnknownMatchSv(tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex],
                                             sa,
                                             mref,
                                             doubleSupportSa);
                    }
                    return;
                }
            }
            auto it = itStart;
            while (it != tumorResults[static_cast<unsigned int>(saChrIndex)].begin()) {
                ChrDistance distance = it->distanceToSupp(sa);
                if (distance <= fuzziness) {
                    dbHits.emplace_back(distance, it);
                    --it;
                } else {
                    break;
                }
            }
            if (itStart != tumorResults[static_cast<unsigned int>(saChrIndex)].end()) {
                auto it = next(itStart);
                while (it != tumorResults[static_cast<unsigned int>(saChrIndex)].end()) {
                    auto distance = it->distanceToSupp(sa);
                    if (distance <= fuzziness) {
                        dbHits.emplace_back(distance, it);
                        ++it;
                    } else {
                        break;
                    }
                }
            }
            sort(dbHits.begin(), dbHits.end());
        }
        if (dbHits.empty()) {
            if (createUnknownMatchSvPreCheck(sa, doubleSupportSa)) {
                createUnknownMatchSv(tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex],
                                     sa,
                                     mref,
                                     doubleSupportSa);
            }
            return;
        } else {
            auto createdMatch = false;
            for (auto &res : dbHits) {
                for (const auto &saMatch : res.second->getSuppAlignments()) {
                    if (tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex].closeToSupp(
                            saMatch, fuzziness)) {
                        if (createDoubleMatchSvPreCheck(saMatch)) {
                            createDoubleMatchSv(tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex],
                                                *res.second, sa, saMatch, true,
                                                mref);
                            createdMatch = true;
                        }
                    }
                }
            }
            if (!createdMatch) {
                auto res = dbHits[0];
                for (const auto &saMatch : res.second->getSuppAlignments()) {
                    if (tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex].closeToSupp(
                            saMatch, fuzziness * 3)) {
                        if (createDoubleMatchSvPreCheck(saMatch)) {
                            createDoubleMatchSv(tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex],
                                                *res.second, sa, saMatch, true,
                                                mref);
                            return;
                        }
                    }
                }
                createUnmatchedSaSv(tumorResults[static_cast<unsigned int>(compressedMrefIndex)][dbIndex], *res.second,
                                    sa, mref);
            }
        }
    }

    void
    AnnotationProcessor::createDoubleMatchSv(BreakpointReduced &sourceBp,
                                             BreakpointReduced &targetBp,
                                             const SuppAlignmentAnno &sa,
                                             const SuppAlignmentAnno &saMatch,
                                             bool checkOrder,
                                             std::vector<std::vector<MrefEntryAnno>> &mref) {
        if (checkOrder) {
            if (sourceBp.getMrefHits().getNumConsevativeHits() == -1) {
                auto germlineInfo = searchGermlineHitsNew(
                    sourceBp,
                    SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                    SvEvent::GERMLINE_OFFSET_THRESHOLD);
                sourceBp.setGermlineInfo(germlineInfo);
                auto mrefInfo = searchMrefHitsNew(
                    sourceBp,
                    SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                    SvEvent::GERMLINE_OFFSET_THRESHOLD,
                    mref);
                sourceBp.setMrefHits(mrefInfo);
            }
            if (targetBp.getMrefHits().getNumConsevativeHits() == -1) {
                auto germlineInfo = searchGermlineHitsNew(
                    targetBp,
                    SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                    SvEvent::GERMLINE_OFFSET_THRESHOLD);
                targetBp.setGermlineInfo(germlineInfo);
                auto mrefInfo = searchMrefHitsNew(
                    targetBp,
                    SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                    SvEvent::GERMLINE_OFFSET_THRESHOLD,
                    mref);
                targetBp.setMrefHits(mrefInfo);
            }
            if (targetBp.fullSmaller(sourceBp)) {
                createDoubleMatchSv(targetBp, sourceBp, saMatch, sa, false, mref);
            }
        }
        visitedLineIndices.push_back(sourceBp.getLineIndex());
        visitedLineIndices.push_back(targetBp.getLineIndex());
        filteredResults.emplace_back(sourceBp, targetBp, sa, saMatch, overhangs);
        checkSvQuality();
    }
    bool
    AnnotationProcessor::createDoubleMatchSvPreCheck(
        const SuppAlignmentAnno &saMatch) {
        if (SvEvent::DEBUG_MODE || !saMatch.isSuspicious()) {
            return true;
        }
        return false;
    }
    void
    AnnotationProcessor::createUnmatchedSaSv(BreakpointReduced &sourceBp,
                                             BreakpointReduced &targetBp,
                                             const SuppAlignmentAnno &sa,
                                             std::vector<std::vector<MrefEntryAnno>> &mref) {
        if (sourceBp.getMrefHits().getNumConsevativeHits() == -1) {
            auto germlineInfo = searchGermlineHitsNew(
                sourceBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                SvEvent::GERMLINE_OFFSET_THRESHOLD);
            sourceBp.setGermlineInfo(germlineInfo);
            auto mrefInfo = searchMrefHitsNew(
                sourceBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                SvEvent::GERMLINE_OFFSET_THRESHOLD, mref);
            sourceBp.setMrefHits(mrefInfo);
        }
        targetBp.addDummySa(sa, sourceBp);
        auto germlineInfo = searchGermlineHitsNew(
            targetBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
            SvEvent::GERMLINE_OFFSET_THRESHOLD);
        auto mrefHits =
            searchMrefHitsNew(targetBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                              SvEvent::GERMLINE_OFFSET_THRESHOLD, mref);
        targetBp.setMrefHits(mrefHits);
        targetBp.setGermlineInfo(germlineInfo);
        visitedLineIndices.push_back(sourceBp.getLineIndex());
        visitedLineIndices.push_back(targetBp.getLineIndex());
        filteredResults.emplace_back(sourceBp, targetBp, sa, overhangs,
                                     targetBp.getDummySa());
        checkSvQuality();
        targetBp.removeMarkedFuzzies();
    }
    bool
    AnnotationProcessor::createUnknownMatchSvPreCheck(const SuppAlignmentAnno &sa,
                                                      bool doubleSupportSa) {
        if (SvEvent::DEBUG_MODE || !sa.isSemiSuspicious()) {
            if (doubleSupportSa ||
                (!sa.isFuzzy() &&
                 (sa.getSupport() > 0 || sa.getSecondarySupport() > 0))) {
                return true;
            }
        }
        return false;
    }
    void
    AnnotationProcessor::createUnknownMatchSv(BreakpointReduced &sourceBp,
                                              const SuppAlignmentAnno &sa,
                                              std::vector<std::vector<MrefEntryAnno>> &mref,
                                              bool doubleSupportSa[[gnu::unused]] // TODO: remove
                                              ) {
        auto germlineInfo = searchGermlineHitsNew(
            sourceBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
            SvEvent::GERMLINE_OFFSET_THRESHOLD);
        sourceBp.setGermlineInfo(germlineInfo);
        auto mrefInfo =
            searchMrefHitsNew(sourceBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                              SvEvent::GERMLINE_OFFSET_THRESHOLD, mref);
        sourceBp.setMrefHits(mrefInfo);
        BreakpointReduced dummyBp{sa, sourceBp, false};
        auto dummyGermline =
            searchGermlineHitsNew(dummyBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                                  SvEvent::GERMLINE_OFFSET_THRESHOLD * 3);
        auto dummyMref =
            searchMrefHitsNew(dummyBp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                              SvEvent::GERMLINE_OFFSET_THRESHOLD * 3, mref);
        if (sa.isFuzzy() &&
            static_cast<long>(sa.getExtendedPos()) - static_cast<long>(sa.getPos()) > 3 * SvEvent::GERMLINE_OFFSET_THRESHOLD) {
            BreakpointReduced dummyBp2{sa, sourceBp, true};
            auto dummyGermline2 = searchGermlineHitsNew(
                dummyBp2, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                SvEvent::GERMLINE_OFFSET_THRESHOLD * 3);
            auto dummyMref2 = searchMrefHitsNew(
                dummyBp2, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                SvEvent::GERMLINE_OFFSET_THRESHOLD * 3, mref);
            if (dummyMref2.getNumConsevativeHits() >
                dummyMref.getNumConsevativeHits()) {
                dummyMref = dummyMref2;
                dummyGermline = dummyGermline2;
            }
        }
        visitedLineIndices.push_back(sourceBp.getLineIndex());
        filteredResults.emplace_back(sourceBp, sa, dummyGermline, dummyMref,
                                     overhangs, dummyBp.getDummySa());
        checkSvQuality();
    }
    void
    AnnotationProcessor::checkSvQuality() {
        auto key = filteredResults.back().getKey();
        if (!key.empty()) {
            if (filteredResultKeys.count(key) == 0) {
                filteredResultKeys.insert(key);
                return;
            }
        }
        filteredResults.pop_back();
    }

    MrefMatch
    AnnotationProcessor::searchMrefHitsNew(const BreakpointReduced &bpIn,
                                           int distanceThreshold,
                                           int conservativeDistanceThreshold[[gnu::unused]], // TODO: remove
                                           std::vector<std::vector<MrefEntryAnno>> &mref) {
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        std::vector<SuppAlignmentAnno> suppMatches{};
        if (!chrConverter.isCompressedMref(bpIn.getChrIndex())) {
            return MrefMatch{0, 0, 10000, suppMatches};
        }
        unsigned int compressedMrefIndex =
            static_cast<unsigned int>(chrConverter.indexToCompressedMrefIndex(bpIn.getChrIndex()));

        auto itStart = lower_bound(mref[compressedMrefIndex].begin(),
                                   mref[compressedMrefIndex].end(), bpIn);
        if (itStart == mref[compressedMrefIndex].end()) {
            return MrefMatch{0, 0, 10000, suppMatches};
        }
        if (itStart != mref[compressedMrefIndex].begin() &&
            !(itStart->distanceTo(bpIn) < SvEvent::GERMLINE_OFFSET_THRESHOLD) &&
            prev(itStart)->distanceTo( bpIn) < SvEvent::GERMLINE_OFFSET_THRESHOLD) {
            --itStart;
        }
        auto it = itStart;

        std::vector<std::vector<MrefEntryAnno>::iterator> dbHits{};
        std::vector<std::vector<MrefEntryAnno>::iterator> dbHitsConservative{};
        while (true) {
            auto tmpDistance = it->distanceTo(bpIn);
            if (tmpDistance < SvEvent::GERMLINE_OFFSET_THRESHOLD) {
                dbHitsConservative.push_back(it);
            }
            if (tmpDistance < distanceThreshold) {
                dbHits.push_back(it);
            } else {
                break;
            }
            if (it == mref[compressedMrefIndex].begin()) {
                break;
            }
            --it;
        }
        if (itStart != mref[compressedMrefIndex].end()) {
            it = next(itStart);
            while (true) {
                if (it == mref[compressedMrefIndex].end()) {
                    break;
                }
                auto tmpDistance = it->distanceTo(bpIn);
                if (tmpDistance < SvEvent::GERMLINE_OFFSET_THRESHOLD) {
                    dbHitsConservative.push_back(it);
                }
                if (tmpDistance < distanceThreshold) {
                    dbHits.push_back(it);
                } else {
                    break;
                }
                ++it;
            }
        }
        if (dbHits.empty()) {
            return MrefMatch{0, 0, 0, suppMatches};
        }
        short score{0};
        auto offset = 0;
        for (auto res : dbHits) {
            auto saMatch = false;
            for (const auto &saRef : res->getSuppAlignments()) {
                for (const auto &sa : bpIn.getSuppAlignments()) {
                    if (saRef.saCloseness(sa, SuppAlignmentAnno::DEFAULT_READ_LENGTH /
                                                  2)) {
                        saMatch = true;
                        auto previouslyRecorded = false;
                        for (auto &saTmp : suppMatches) {
                            if (saRef.saCloseness(
                                    saTmp,
                                    SuppAlignmentAnno::DEFAULT_READ_LENGTH / 2)) {
                                previouslyRecorded = true;
                                if (saTmp.isFuzzy()) {
                                    if (saRef.isFuzzy()) {
                                        saTmp.extendSuppAlignment(
                                            saRef.getPos(), saRef.getExtendedPos());
                                    } else {
                                        saTmp.removeFuzziness(saRef);
                                    }
                                }
                                if (saRef.getSupport() > saTmp.getSupport()) {
                                    saTmp.setSupport(saRef.getSupport());
                                }
                                if (saRef.getSecondarySupport() >
                                    saTmp.getSecondarySupport()) {
                                    saTmp.setSecondarySupport(
                                        saRef.getSecondarySupport());
                                }
                                if (saTmp.getSecondarySupport() <
                                    saTmp.getSupport()) {
                                    saTmp.setSecondarySupport(saTmp.getSupport());
                                }
                            }
                        }
                        if (!previouslyRecorded) {
                            suppMatches.push_back(saRef);
                        }
                    }
                }
            }
            if (saMatch) {
                auto tmpScore = res->getNumHits();
                if (tmpScore > score) {
                    score = tmpScore;
                    offset = res->distanceTo(bpIn);
                }
            }
        }
        short conservativeScore{0};
        for (const auto res : dbHitsConservative) {
            auto tmpScore = res->getNumHits();
            if (tmpScore < SvEvent::RELAXED_BP_FREQ_THRESHOLD) {
                auto sas = res->getSuppAlignments();
                if (sas.size() == 1) {
                    if ((sas[0].getSupport() + 0.0) / tmpScore > 0.8) {
                        continue;
                    }
                }
            }
            if (tmpScore > conservativeScore) {
                conservativeScore = tmpScore;
            }
        }
        return MrefMatch{std::max(score, conservativeScore),
                         conservativeScore,
                         offset,
                         suppMatches};
    }

    GermlineMatch
    AnnotationProcessor::searchGermlineHitsNew(const BreakpointReduced &bpIn,
                                               int distanceThreshold,
                                               int conservativeDistanceThreshold) {
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();

        GermlineMatch dummyMatchTrue{1.0, 1.0,
                                     std::vector<std::pair<SuppAlignmentAnno, double>>{}};
        GermlineMatch dummyMatchFalse{0.0, 0.0,
                                      std::vector<std::pair<SuppAlignmentAnno, double>>{}};
        if (NO_CONTROL_MODE) {
            return dummyMatchTrue;
        }

        if (!chrConverter.isCompressedMref(bpIn.getChrIndex())) {
            return dummyMatchFalse;
        }

        unsigned int compressedMrefIndex =
            static_cast<unsigned int>(GlobalAppConfig::getInstance().getChrConverter().
                                        indexToCompressedMrefIndex(bpIn.getChrIndex()));

        if (controlResults[compressedMrefIndex].empty()) {
            return dummyMatchFalse;
        }
        auto itStart = lower_bound(controlResults[compressedMrefIndex].begin(),
                                   controlResults[compressedMrefIndex].end(), bpIn);
        if (itStart == controlResults[compressedMrefIndex].end()) {
            return dummyMatchFalse;
        }
        if (itStart != controlResults[compressedMrefIndex].cbegin() &&
            !(itStart->distanceToBp(bpIn) < SvEvent::GERMLINE_OFFSET_THRESHOLD) &&
            prev(itStart)->distanceToBp(bpIn) < SvEvent::GERMLINE_OFFSET_THRESHOLD) {
            --itStart;
        }
        auto it = itStart;
        std::vector<std::vector<BreakpointReduced>::iterator> dbHits{};
        std::vector<std::vector<BreakpointReduced>::iterator> dbHitsConservative{};

        while (true) {
            auto tmpDistance = it->distanceToBp(bpIn);
            if (tmpDistance < 0) {
                break;
            }
            if (tmpDistance < conservativeDistanceThreshold) {
                dbHitsConservative.push_back(it);
            }
            if (tmpDistance < distanceThreshold) {
                dbHits.push_back(it);
            } else {
                break;
            }
            if (it == controlResults[compressedMrefIndex].begin()) {
                break;
            }
            --it;
        }

        if (itStart != controlResults[compressedMrefIndex].end()) {
            it = next(itStart);
            while (true) {
                if (it == controlResults[compressedMrefIndex].end()) {
                    break;
                }
                auto tmpDistance = it->distanceToBp(bpIn);
                if (tmpDistance < 0) {
                    break;
                }
                if (tmpDistance < conservativeDistanceThreshold) {
                    dbHitsConservative.push_back(it);
                }
                if (tmpDistance < distanceThreshold) {
                    dbHits.push_back(it);
                } else {
                    break;
                }
                ++it;
            }
        }
        if (dbHits.empty()) {
            return dummyMatchFalse;
        }
        auto conservativeClonality = 0.0;
        if (!dbHitsConservative.empty()) {
            auto maxSupport = 0.0;
            for (const auto res : dbHitsConservative) {
                auto mateSupSa = 0;
                for (const SuppAlignmentAnno &sa : res->getSuppAlignments()) {
                    if (sa.getMateSupport() > mateSupSa) {
                        mateSupSa = sa.getMateSupport();
                    }
                }
                if (res->getMateSupport() > mateSupSa) {
                    mateSupSa = res->getMateSupport();
                }
                auto support =
                    0.0 + res->getPairedBreaksSoft() + res->getBreaksShortIndel() +
                    res->getPairedBreaksHard() + res->getUnpairedBreaksSoft() +
                    res->getUnpairedBreaksHard() + mateSupSa;
                if (support > maxSupport) {
                    maxSupport = support;
                    if (support + res->getNormalSpans() > 0) {
                        conservativeClonality =
                            support / (support + res->getNormalSpans());
                    }
                }
            }
        }
        auto clonality = conservativeClonality;
        std::vector<std::pair<SuppAlignmentAnno, double>> suppMatches{};
        for (auto res : dbHits) {
            auto breakSupportSoft =
                res->getPairedBreaksSoft() + res->getUnpairedBreaksSoft();
            auto breakSupportHard =
                res->getPairedBreaksHard() + res->getUnpairedBreaksHard();
            for (const auto &saRef : res->getSuppAlignments()) {
                for (const auto &sa : bpIn.getSuppAlignments()) {
                    if (saRef.saCloseness(sa, SuppAlignmentAnno::DEFAULT_READ_LENGTH /
                                                  2)) {
                        auto previouslyRecorded = false;
                        for (auto &saTmp : suppMatches) {
                            if (saRef.saClosenessDirectional(
                                    saTmp.first,
                                    SuppAlignmentAnno::DEFAULT_READ_LENGTH / 2)) {
                                previouslyRecorded = true;
                                if (saTmp.first.isFuzzy()) {
                                    if (saRef.isFuzzy()) {
                                        saTmp.first.extendSuppAlignment(
                                            saRef.getPos(), saRef.getExtendedPos());
                                    } else {
                                        saTmp.first.removeFuzziness(saRef);
                                    }
                                }
                                if (saRef.getSupport() > saTmp.first.getSupport()) {
                                    saTmp.first.setSupport(saRef.getSupport());
                                }
                                if (saRef.getSecondarySupport() >
                                    saTmp.first.getSecondarySupport()) {
                                    saTmp.first.setSecondarySupport(
                                        saRef.getSecondarySupport());
                                }
                                if (saRef.getMateSupport() >
                                    saTmp.first.getMateSupport()) {
                                    saTmp.first.setMateSupport(
                                        saRef.getMateSupport());
                                    saTmp.first.setExpectedDiscordants(
                                        saRef.getExpectedDiscordants());
                                }
                                auto currentSoftSupport =
                                std::max(saTmp.first.getSupport(), breakSupportSoft);
                                auto currentHardSupport =
                                std::max(saTmp.first.getSecondarySupport(),
                                        breakSupportHard);
                                auto breakSupport =
                                    currentSoftSupport + currentHardSupport + 0.0;
                                auto currentClonality =
                                    (breakSupport + saTmp.first.getMateSupport()) /
                                    (breakSupport + saTmp.first.getMateSupport() +
                                     res->getNormalSpans());
                            saTmp.second = std::max(currentClonality, saTmp.second);
                            }
                        }
                        if (!previouslyRecorded) {
                            auto currentSoftSupport =
                            std::max(saRef.getSupport(), breakSupportSoft);
                            auto currentHardSupport =
                            std::max(saRef.getSecondarySupport(), breakSupportHard);
                            auto breakSupport =
                                currentSoftSupport + currentHardSupport + 0.0;
                            auto currentClonality =
                                (breakSupport + saRef.getMateSupport()) /
                                (breakSupport + saRef.getMateSupport() +
                                 res->getNormalSpans());
                            suppMatches.push_back({saRef, currentClonality});
                        }
                    }
                }
            }
        }
        return GermlineMatch{clonality, conservativeClonality, suppMatches};
    }

    bool
    AnnotationProcessor::applyMassiveInversionFiltering(bool stricterMode,
                                                        bool controlCheckMode) {
        auto deletionCandidateCount = 0;
        auto totalCount = 0;
        for (const auto &sv : filteredResults) {
            if (sv.isToRemove()) {
                continue;
            }
            if (sv.getContaminationCandidate() > 0) {
                continue;
            }
            if (sv.getSuspicious() == 0 && sv.getEventScore() > 2) {
                if (sv.isInverted() || sv.isSemiSuspicious()) {
                    ++deletionCandidateCount;
                } else if (stricterMode && sv.getEvidenceLevel2() == 0) {
                    ++deletionCandidateCount;
                }
                ++totalCount;
            }
        }
        auto totalCountThreshold = 300;
        if (stricterMode || controlCheckMode) {
            totalCountThreshold = 200;
        }
        auto invRatio = (deletionCandidateCount + 0.0) / totalCount;
        if (invRatio > 0.7 && totalCount > totalCountThreshold) {
            for (auto &sv : filteredResults) {
                if (sv.isToRemove()) {
                    continue;
                }
                if (sv.getContaminationCandidate() > 0) {
                    continue;
                }
                if (sv.getSuspicious() != 0) {
                    continue;
                }
                if (sv.isDistant()) {
                    if (sv.getSelectedSa1().isProperPairErrorProne() &&
                        sv.getSelectedSa1().getMateSupport() < 6) {
                        sv.setToRemove(true);
                        continue;
                    }
                    if (sv.getSelectedSa2().isProperPairErrorProne() &&
                        sv.getSelectedSa2().getMateSupport() < 6) {
                        sv.setToRemove(true);
                        continue;
                    }
                    if (sv.isSemiSuspicious()) {
                        if (sv.getEventScore() < 3) {
                            sv.setToRemove(true);
                            continue;
                        }
                        if (sv.getSelectedSa1().isSemiSuspicious() &&
                            sv.getSelectedSa2().isSemiSuspicious()) {
                            sv.setToRemove(true);
                            continue;
                        }
                        if (sv.getEvidenceLevel2() < 3 &&
                            sv.getSelectedSa1().isSemiSuspicious()) {
                            if (sv.getEvidenceLevel1() < 3 ||
                                sv.getSelectedSa1().getSupport() < 5) {
                                sv.setToRemove(true);
                                continue;
                            }
                        }
                        if (sv.getEvidenceLevel1() < 3 &&
                            sv.getSelectedSa2().isSemiSuspicious()) {
                            if (sv.getEvidenceLevel2() < 3 ||
                                sv.getSelectedSa2().getSupport() < 5) {
                                sv.setToRemove(true);
                                continue;
                            }
                        }
                    }
                } else {
                    if (sv.getEventScore() == 3 || sv.getEvidenceLevel2() == 0) {
                        sv.setToRemove(true);
                        continue;
                    }
                }
                if (sv.isInverted()) {
                    if (stricterMode) {
                        if (sv.getSelectedSa1().isStrictFuzzyCandidate() ||
                            sv.getSelectedSa2().isStrictFuzzyCandidate()) {
                            sv.setToRemove(true);
                            continue;
                        }
                        if (sv.isSemiSuspicious()) {
                            sv.setToRemove(true);
                            continue;
                        }
                    } else {
                        if (sv.getSelectedSa1().isStrictFuzzyCandidate() &&
                            sv.getSelectedSa2().isStrictFuzzyCandidate()) {
                            sv.setToRemove(true);
                            continue;
                        }
                    }
                    if (sv.getEventScore() < 3) {
                        sv.setToRemove(true);
                        continue;
                    }
                    if (sv.getMateRatio1() < 0.6 && sv.getMateRatio2() < 0.6) {
                        sv.setToRemove(true);
                        continue;
                    }
                    if (sv.getEventSize() > 0 && sv.getEventSize() < 10000) {
                        if (stricterMode || sv.getTotalEvidence1() < 5 ||
                            sv.getTotalEvidence2() < 5) {
                            sv.setToRemove(true);
                            continue;
                        }
                    } else {
                        if (sv.getTotalEvidence1() < 5 &&
                            sv.getTotalEvidence2() < 5) {
                            sv.setToRemove(true);
                            continue;
                        }
                    }
                }
                if (stricterMode) {
                    if (sv.getEvidenceLevel2() == 0) {
                        if (sv.getSelectedSa1().getSupport() < 10 ||
                            !sv.isOverhang1Compensation()) {
                            sv.setToRemove(true);
                            continue;
                        }
                    }
                    if (sv.getEventScore() == 3) {
                        if (sv.getSelectedSa1().getSupport() < 10) {
                            sv.setToRemove(true);
                            continue;
                        }
                        if (sv.getSelectedSa1().getSecondarySupport() < 5) {
                            sv.setToRemove(true);
                            continue;
                        }
                        if (sv.isDistant()) {
                            if (sv.getSelectedSa1().getMateSupport() < 10) {
                                sv.setToRemove(true);
                                continue;
                            }
                        }
                    }
                    if (sv.isDistant()) {
                        if (sv.getSelectedSa1().isProperPairErrorProne() &&
                            (sv.getSelectedSa2().isProperPairErrorProne() ||
                             sv.getEvidenceLevel2() == 0)) {
                            if (sv.getSelectedSa1().getSupport() < 5 ||
                                sv.getSelectedSa2().getSupport() < 5) {
                                sv.setToRemove(true);
                                continue;
                            }
                            if (sv.getSelectedSa1().getSecondarySupport() < 3 &&
                                sv.getSelectedSa2().getSecondarySupport() < 3) {
                                sv.setToRemove(true);
                                continue;
                            }
                        }
                        if (sv.getSelectedSa1().isSemiSuspicious()) {
                            if (sv.getEvidenceLevel1() < 3 ||
                                sv.getSelectedSa1().getSupport() < 5) {
                                sv.setToRemove(true);
                                continue;
                            }
                        }
                        if (sv.getSelectedSa2().isSemiSuspicious()) {
                            if (sv.getEvidenceLevel2() < 3 ||
                                sv.getSelectedSa2().getSupport() < 5) {
                                sv.setToRemove(true);
                                continue;
                            }
                        }

                        if (sv.getEvidenceLevel2() == 0) {
                            if (sv.getEvidenceLevel1() < 3 ||
                                sv.getSelectedSa1().getSupport() < 5) {
                                sv.setToRemove(true);
                                continue;
                            }
                        }
                        if (sv.getEvidenceLevel1() < 3 &&
                            sv.getEvidenceLevel2() < 3) {
                            sv.setToRemove(true);
                            continue;
                        }
                    } else {
                        if (sv.getEvidenceLevel1() < 2 ||
                            sv.getEvidenceLevel2() < 2) {
                            sv.setToRemove(true);
                            continue;
                        }
                        if (sv.getTotalEvidence1() < 5 &&
                            sv.getTotalEvidence2() < 5) {
                            sv.setToRemove(true);
                            continue;
                        }
                    }
                    if (sv.getTotalEvidence2() < 5) {
                        auto evidenceLevelThreshold = sv.isDistant() ? 3 : 2;
                        if (sv.getEvidenceLevel1() < evidenceLevelThreshold) {
                            sv.setToRemove(true);
                            continue;
                        }
                    }
                }
            }
            return true;
        }
        return false;
    }
    bool
    AnnotationProcessor::applyPathogenContaminationFiltering() {
        auto likelyContaminants = 0;
        for (auto &sv : filteredResults) {
            if (sv.isToRemove()) {
                continue;
            }
            if (sv.getContaminationCandidate() != 2) {
                continue;
            }
            if (sv.getSuspicious() != 0) {
                continue;
            }
            if (sv.getEventScore() < 3) {
                continue;
            }
            ++likelyContaminants;
            if (likelyContaminants > 9) {
                break;
            }
        }
        if (likelyContaminants > 9) {
            auto cleanedContaminants = 0;
            for (auto &sv : filteredResults) {
                if (sv.isToRemove()) {
                    continue;
                }
                if (sv.getSuspicious() != 0) {
                    continue;
                }
                if (sv.getEventScore() < 3) {
                    continue;
                }
                if (sv.getContaminationCandidate() == 2) {
                    ++cleanedContaminants;
                }
            }
            if (cleanedContaminants > 19) {
                for (auto &sv : filteredResults) {
                    if (sv.isToRemove()) {
                        continue;
                    }
                    if (sv.getSuspicious() != 0) {
                        continue;
                    }
                    if (sv.getEventScore() < 2) {
                        continue;
                    }
                    if (sv.getContaminationCandidate() > 0) {
                        sv.setEventScore(1);
                        sv.setEventType(5);
                        continue;
                    }
                    if (sv.getSelectedSa1().getSupport() > 19 &&
                        sv.getSelectedSa1().getSecondarySupport() < 3) {
                        sv.setEventScore(1);
                        sv.setEventType(5);
                        continue;
                    }
                    if (sv.getSelectedSa2().getSupport() > 19 &&
                        sv.getSelectedSa2().getSecondarySupport() < 3) {
                        sv.setEventScore(1);
                        sv.setEventType(5);
                        continue;
                    }
                    if (sv.getOverhang1lengthRatio() > 0.7 ||
                        sv.getOverhang2lengthRatio() > 0.7) {
                        sv.setEventScore(1);
                        sv.setEventType(5);
                        continue;
                    }
                }
            }
            return true;
        }
        return false;
    }
    void
    AnnotationProcessor::printFilteredResults(bool contaminationInControl,
                                              int controlPrefilteringLevel) const {
        if (controlPrefilteringLevel > 0) {
            std::cout << "#controlMassiveInvPrefilteringLevel\t"
                 << controlPrefilteringLevel << std::endl;
        }
        if (massiveInvFilteringLevel > 0) {
            std::cout << "#tumorMassiveInvFilteringLevel\t" << massiveInvFilteringLevel
                 << std::endl;
        }
        if (contaminationInControl) {
            std::cout << "#likelyPathogenInControl\tTRUE" << std::endl;
        }
        if (contaminationObserved) {
            std::cout << "#likelyPathogenInTumor\tTRUE" << std::endl;
        }
        for (const auto &sv : filteredResults) {
            if (!sv.isToRemove()) {
                std::cout << sv.printMatch(overhangs);
            }
        }
    }

    void
    AnnotationProcessor::printUnresolvedRareOverhangs(
        std::vector<std::vector<MrefEntryAnno>> &mref) {
        if (massiveInvFilteringLevel != 0) {
            return;
        }
        sort(visitedLineIndices.begin(), visitedLineIndices.end());
        std::unordered_set<int> visitedLineIndicesSet{visitedLineIndices.begin(),
                                                      visitedLineIndices.end()};
        for (const auto &tumorChromosome : tumorResults) {
            for (const auto &bp : tumorChromosome) {
                if (visitedLineIndicesSet.count(bp.getLineIndex())) {
                    continue;
                }
                if (bp.testOverhangBasedCandidacy()) {
                    auto mrefHits = searchMrefHitsNew(
                        bp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                        SvEvent::GERMLINE_OFFSET_THRESHOLD, mref);
                    if (mrefHits.getNumHits() < SvEvent::GERMLINE_DB_LIMIT) {
                        auto germlineClonality = 1.0;
                        if (!NO_CONTROL_MODE) {
                            germlineClonality =
                                searchGermlineHitsNew(
                                    bp, SuppAlignmentAnno::DEFAULT_READ_LENGTH * 6,
                                    SvEvent::GERMLINE_OFFSET_THRESHOLD)
                                    .getClonality();
                        }
                        std::string overhang{""};
                        {
                            std::pair<int, std::string> dummy{bp.getLineIndex(), ""};
                            auto lower = lower_bound(overhangs.cbegin(),
                                                     overhangs.cend(), dummy);
                            if (lower != overhangs.cend()) {
                                if (lower->first == bp.getLineIndex()) {
                                    overhang = lower->second;
                                } else if (next(lower) != overhangs.cend() &&
                                           next(lower)->first ==
                                               bp.getLineIndex()) {
                                    overhang = next(lower)->second;
                                } else if (lower != overhangs.cbegin() &&
                                           prev(lower)->first ==
                                               bp.getLineIndex()) {
                                    overhang = prev(lower)->second;
                                }
                            }
                        }
                        if (!overhang.empty()) {
                            std::cout << bp.printOverhang(
                                germlineClonality, mrefHits.getNumHits(), overhang);
                        }
                    }
                }
            }
        }
    }

} /* namespace sophia */
