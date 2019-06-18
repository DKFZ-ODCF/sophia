/*
 * AnnotationProcessor.cpp
 *
 *  Created on: 28 Apr 2016
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

#include <AnnotationProcessor.h>
#include "Breakpoint.h"
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <DeFuzzier.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "SuppAlignment.h"
#include "HelperFunctions.h"

namespace sophia {
bool AnnotationProcessor::ABRIDGEDOUTPUT { false };

AnnotationProcessor::AnnotationProcessor(const std::string& tumorResultsIn, std::vector<std::vector<MrefEntryAnno>> &mref, int defaultReadLengthTumorIn, bool controlCheckMode, int germlineDbLimit) :
				NOCONTROLMODE { true },
				GERMLINEDBLIMIT { germlineDbLimit },
				contaminationObserved { false },
				massiveInvFilteringLevel { 0 },
				filteredResults { },
				tumorResults { 85, std::vector<BreakpointReduced> { } },
				controlResults { 85, std::vector<BreakpointReduced> { } },
				visitedLineIndices { } {
	std::unique_ptr<std::ifstream> tumorInputHandle { std::make_unique<std::ifstream>(tumorResultsIn, std::ios_base::in | std::ios_base::binary) };
	std::unique_ptr<boost::iostreams::filtering_istream> tumorGzHandle { std::make_unique<boost::iostreams::filtering_istream>() };
	tumorGzHandle->push(boost::iostreams::gzip_decompressor());
	tumorGzHandle->push(*tumorInputHandle);
	std::string line;
	auto lineIndex = 0;
	while (error_terminating_getline(*tumorGzHandle, line)) {
		if (line.front() == '#') {
			continue;
		};
		Breakpoint tmpBp { line, true };
		auto chrIndex = ChrConverter::indexConverter[tmpBp.getChrIndex()];
		if (chrIndex < 0) {
			continue;
		}
		auto hasOverhang = line.back() != '.' && line.back() != '#';
		tumorResults[chrIndex].emplace_back(tmpBp, lineIndex, hasOverhang);
		if (hasOverhang) {
			std::string overhang { };
			for (auto it = line.rbegin(); it != line.rend(); ++it) {
				if (*it == '\t') {
					break;
				} else {
					overhang.push_back(*it);
				}
			}
			std::reverse(overhang.begin(), overhang.end());
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
	if (!controlCheckMode && !contaminationObserved && massiveInvFilteringLevel == 0) {
		printUnresolvedRareOverhangs(mref);
	}
}

AnnotationProcessor::AnnotationProcessor(const std::string& tumorResultsIn, std::vector<std::vector<MrefEntryAnno>> &mref, const std::string& controlResultsIn, int defaultReadLengthTumorIn, int defaultReadLengthControlIn, int germlineDbLimit, int lowQualControlIn, bool pathogenInControlIn) :
				NOCONTROLMODE { false },
				GERMLINEDBLIMIT { germlineDbLimit },
				contaminationObserved { false },
				massiveInvFilteringLevel { 0 },
				filteredResults { },
				tumorResults { 85, std::vector<BreakpointReduced> { } },
				controlResults { 85, std::vector<BreakpointReduced> { } } {
	std::unique_ptr<std::ifstream> controlInputHandle { std::make_unique<std::ifstream>(controlResultsIn, std::ios_base::in | std::ios_base::binary) };
	std::unique_ptr<boost::iostreams::filtering_istream> controlGzHandle { std::make_unique<boost::iostreams::filtering_istream>() };
	controlGzHandle->push(boost::iostreams::gzip_decompressor());
	controlGzHandle->push(*controlInputHandle);
	std::string line;
	auto lineIndex = 0;
	while (error_terminating_getline(*controlGzHandle, line)) {
		if (line.front() == '#') {
			continue;
		};
		Breakpoint tmpBpPre { line, true };
		BreakpointReduced tmpBp { tmpBpPre, lineIndex, false };
		if (tmpBp.getChrIndex() > 1001) {
			continue;
		}
		if (pathogenInControlIn) {
			if ((tmpBp.getPairedBreaksSoft() + tmpBp.getUnpairedBreaksSoft()) > 19 && (tmpBp.getPairedBreaksHard() + tmpBp.getUnpairedBreaksHard() < 3)) {
				if (line.back() != '.' && line.back() != '#') {
					std::string overhang { };
					auto overhangLength = 0;
					auto maxOverhangLength = 0;
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
					auto maxOverhangLengthRatio = (maxOverhangLength + 0.0) / defaultReadLengthControlIn;
					if (maxOverhangLengthRatio > 0.7) {
						continue;
					}
				}
			}
		}
		if (lowQualControlIn > 0) {
			SuppAlignmentAnno* bestSa = nullptr;
			auto bestSaSupport = 0;
			if (!tmpBp.getSuppAlignments().empty()) {
				for (const auto sa : tmpBp.getSupplementsPtr()) {
					if (!sa->isSuspicious()) {
						if (bestSa) {
							auto saSupport = sa->getSupport() + sa->getSecondarySupport() + sa->getMateSupport();
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
			auto clipTotal = tmpBp.getPairedBreaksHard() + tmpBp.getPairedBreaksSoft() + tmpBp.getUnpairedBreaksHard() + tmpBp.getUnpairedBreaksSoft();
			if (!bestSa) {
				if (clipTotal < 10) {
					continue;
				}
			} else {
				if (bestSa->isInverted()) {
					auto suppDist = tmpBp.distanceToSupp(*bestSa);
					if (suppDist > 0) {
						if (suppDist < 10000) {
							if (lowQualControlIn == 2) {
								continue;
							} else {
								if (bestSa->getSupport() < 5 || bestSa->getSecondarySupport() < 5) {
									continue;
								}
							}
						} else {
							if (bestSa->getSupport() < 5 && bestSa->getSecondarySupport() < 5) {
								continue;
							}
						}
					} else if (suppDist < 0) {
						if (bestSa->getSupport() < 5 && bestSa->getSecondarySupport() < 5) {
							continue;
						}
					}
				}
			}
		}
		auto chrIndex = ChrConverter::indexConverter[tmpBp.getChrIndex()];
		controlResults[chrIndex].push_back(tmpBp);
		++lineIndex;
	}
	for (auto &cres : controlResults) {
		DeFuzzier deFuzzierControl { defaultReadLengthControlIn * 6, false };
		deFuzzierControl.deFuzzyDb(cres);
	}
	std::unique_ptr<std::ifstream> tumorInputHandle { std::make_unique<std::ifstream>(tumorResultsIn, std::ios_base::in | std::ios_base::binary) };
	std::unique_ptr<boost::iostreams::filtering_istream> tumorGzHandle { std::make_unique<boost::iostreams::filtering_istream>() };
	tumorGzHandle->push(boost::iostreams::gzip_decompressor());
	tumorGzHandle->push(*tumorInputHandle);
	lineIndex = 0;
	while (error_terminating_getline(*tumorGzHandle, line)) {
		if (line.front() == '#') {
			continue;
		};
		Breakpoint tmpBp { line, true };
		auto chrIndex = ChrConverter::indexConverter[tmpBp.getChrIndex()];
		if (chrIndex < 0) {
			continue;
		}
		auto hasOverhang = line.back() != '.' && line.back() != '#';
		tumorResults[chrIndex].emplace_back(tmpBp, lineIndex, hasOverhang);
		if (line.back() != '.' && line.back() != '#') {
			std::string overhang { };
			for (auto it = line.rbegin(); it != line.rend(); ++it) {
				if (*it == '\t') {
					break;
				} else {
					overhang.push_back(*it);
				}
			}
			std::reverse(overhang.begin(), overhang.end());
			overhangs.emplace_back(lineIndex, overhang);
		} else {
			visitedLineIndices.push_back(lineIndex);
		}
		++lineIndex;
	}
	for (auto &tres : tumorResults) {
		DeFuzzier deFuzzierTumor { defaultReadLengthTumorIn * 6, false };
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

void AnnotationProcessor::searchMatches(std::vector<std::vector<MrefEntryAnno>> &mref) {
	for (auto j = 0; j < 85; ++j) {
		for (auto i = 0u; i < tumorResults[j].size(); ++i) {
			for (const auto &sa : tumorResults[j][i].getSuppAlignments()) {
				if (SvEvent::DEBUGMODE || !sa.isSuspicious()) {
					if (sa.getSecondarySupport() > 0 || (sa.getSupport() > 0 && sa.getMateSupport() > 0)) {
						searchSa(j, i, sa, true, mref);
					} else {
						searchSa(j, i, sa, false, mref);
					}
				}
			}
		}
	}

}

void AnnotationProcessor::searchSa(int chrIndex, int dbIndex, const SuppAlignmentAnno& sa, bool doubleSupportSa, std::vector<std::vector<MrefEntryAnno>> &mref) {
	if (sa.getSupport() + sa.getSecondarySupport() + sa.getMateSupport() < 3) {
		return;
	}
	if (sa.getChrIndex() == 1001) {
		if (createUnknownMatchSvPreCheck(sa, doubleSupportSa)) {
			createUnknownMatchSv(tumorResults[chrIndex][dbIndex], sa, mref, doubleSupportSa);
		}
		return;
	}
	auto saChrIndex = ChrConverter::indexConverter[sa.getChrIndex()];
	if (saChrIndex < 0) {
		return;
	}
	auto fuzziness = 3 * SuppAlignmentAnno::DEFAULTREADLENGTH;
	std::vector<std::pair<int, std::vector<BreakpointReduced>::iterator>> dbHits { };

	if (!tumorResults[saChrIndex].empty()) {
		auto itStart = std::lower_bound(tumorResults[saChrIndex].begin(), tumorResults[saChrIndex].end(), sa);
		if (itStart == tumorResults[saChrIndex].end()) {
			--itStart;
		}
		if (itStart != tumorResults[saChrIndex].begin() && !itStart->closeToSupp(sa, fuzziness)) {
			--itStart;
			if (!itStart->closeToSupp(sa, fuzziness)) {
				if (createUnknownMatchSvPreCheck(sa, doubleSupportSa)) {

					createUnknownMatchSv(tumorResults[chrIndex][dbIndex], sa, mref, doubleSupportSa);
				}
				return;
			}
		}
		auto it = itStart;
		while (it != tumorResults[saChrIndex].begin()) {
			auto distance = it->distanceToSupp(sa);
			if (distance <= fuzziness) {
				dbHits.emplace_back(distance, it);
				--it;
			} else {
				break;
			}
		}
		if (itStart != tumorResults[saChrIndex].end()) {
			auto it = std::next(itStart);
			while (it != tumorResults[saChrIndex].end()) {
				auto distance = it->distanceToSupp(sa);
				if (distance <= fuzziness) {
					dbHits.emplace_back(distance, it);
					++it;
				} else {
					break;
				}
			}
		}
		std::sort(dbHits.begin(), dbHits.end());
	}
	if (dbHits.empty()) {
		if (createUnknownMatchSvPreCheck(sa, doubleSupportSa)) {
			createUnknownMatchSv(tumorResults[chrIndex][dbIndex], sa, mref, doubleSupportSa);
		}
		return;
	} else {
		auto createdMatch = false;
		for (auto &res : dbHits) {
			for (const auto &saMatch : res.second->getSuppAlignments()) {
				if (tumorResults[chrIndex][dbIndex].closeToSupp(saMatch, fuzziness)) {
					if (createDoubleMatchSvPreCheck(saMatch)) {
						createDoubleMatchSv(tumorResults[chrIndex][dbIndex], *res.second, sa, saMatch, true, mref);
						createdMatch = true;
					}
				}
			}
		}
		if (!createdMatch) {
			auto res = dbHits[0];
			for (const auto &saMatch : res.second->getSuppAlignments()) {
				if (tumorResults[chrIndex][dbIndex].closeToSupp(saMatch, fuzziness * 3)) {
					if (createDoubleMatchSvPreCheck(saMatch)) {
						createDoubleMatchSv(tumorResults[chrIndex][dbIndex], *res.second, sa, saMatch, true, mref);
						return;
					}
				}
			}
			createUnmatchedSaSv(tumorResults[chrIndex][dbIndex], *res.second, sa, mref);
		}
	}
}

void AnnotationProcessor::createDoubleMatchSv(BreakpointReduced& sourceBp, BreakpointReduced& targetBp, const SuppAlignmentAnno& sa, const SuppAlignmentAnno& saMatch, bool checkOrder, std::vector<std::vector<MrefEntryAnno>>& mref) {
	if (checkOrder) {
		if (sourceBp.getMrefHits().getNumConsevativeHits() == -1) {
			auto germlineInfo = searchGermlineHitsNew(sourceBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD);
			sourceBp.setGermlineInfo(germlineInfo);
			auto mrefInfo = searchMrefHitsNew(sourceBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD, mref);
			sourceBp.setMrefHits(mrefInfo);
		}
		if (targetBp.getMrefHits().getNumConsevativeHits() == -1) {
			auto germlineInfo = searchGermlineHitsNew(targetBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD);
			targetBp.setGermlineInfo(germlineInfo);
			auto mrefInfo = searchMrefHitsNew(targetBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD, mref);
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
bool AnnotationProcessor::createDoubleMatchSvPreCheck(const SuppAlignmentAnno& saMatch) {
	if (SvEvent::DEBUGMODE || !saMatch.isSuspicious()) {
		return true;
	}
	return false;
}
void AnnotationProcessor::createUnmatchedSaSv(BreakpointReduced& sourceBp, BreakpointReduced& targetBp, const SuppAlignmentAnno& sa, std::vector<std::vector<MrefEntryAnno>>& mref) {
	if (sourceBp.getMrefHits().getNumConsevativeHits() == -1) {
		auto germlineInfo = searchGermlineHitsNew(sourceBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD);
		sourceBp.setGermlineInfo(germlineInfo);
		auto mrefInfo = searchMrefHitsNew(sourceBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD, mref);
		sourceBp.setMrefHits(mrefInfo);
	}
	targetBp.addDummySa(sa, sourceBp);
	auto germlineInfo = searchGermlineHitsNew(targetBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD);
	auto mrefHits = searchMrefHitsNew(targetBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD, mref);
	targetBp.setMrefHits(mrefHits);
	targetBp.setGermlineInfo(germlineInfo);
	visitedLineIndices.push_back(sourceBp.getLineIndex());
	visitedLineIndices.push_back(targetBp.getLineIndex());
	filteredResults.emplace_back(sourceBp, targetBp, sa, overhangs, targetBp.getDummySa());
	checkSvQuality();
	targetBp.removeMarkedFuzzies();
}
bool AnnotationProcessor::createUnknownMatchSvPreCheck(const SuppAlignmentAnno& sa, bool doubleSupportSa) {
	if (SvEvent::DEBUGMODE || !sa.isSemiSuspicious()) {
		if (doubleSupportSa || (!sa.isFuzzy() && (sa.getSupport() > 0 || sa.getSecondarySupport() > 0))) {
			return true;
		}
	}
	return false;
}
void AnnotationProcessor::createUnknownMatchSv(BreakpointReduced& sourceBp, const SuppAlignmentAnno& sa, std::vector<std::vector<MrefEntryAnno>>& mref, bool doubleSupportSa) {
	auto germlineInfo = searchGermlineHitsNew(sourceBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD);
	sourceBp.setGermlineInfo(germlineInfo);
	auto mrefInfo = searchMrefHitsNew(sourceBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD, mref);
	sourceBp.setMrefHits(mrefInfo);
	BreakpointReduced dummyBp { sa, sourceBp, false };
	auto dummyGermline = searchGermlineHitsNew(dummyBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD * 3);
	auto dummyMref = searchMrefHitsNew(dummyBp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD * 3, mref);
	if (sa.isFuzzy() && sa.getExtendedPos() - sa.getPos() > 3 * SvEvent::GERMLINEOFFSETTHRESHOLD) {
		BreakpointReduced dummyBp2 { sa, sourceBp, true };
		auto dummyGermline2 = searchGermlineHitsNew(dummyBp2, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD * 3);
		auto dummyMref2 = searchMrefHitsNew(dummyBp2, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD * 3, mref);
		if (dummyMref2.getNumConsevativeHits() > dummyMref.getNumConsevativeHits()) {
			dummyMref = dummyMref2;
			dummyGermline = dummyGermline2;
		}
	}
	visitedLineIndices.push_back(sourceBp.getLineIndex());
	filteredResults.emplace_back(sourceBp, sa, dummyGermline, dummyMref, overhangs, dummyBp.getDummySa());
	checkSvQuality();
}
void AnnotationProcessor::checkSvQuality() {
	auto key = filteredResults.back().getKey();
	if (!key.empty()) {
		if (filteredResultKeys.count(key) == 0) {
			filteredResultKeys.insert(key);
			return;
		}
	}
	filteredResults.pop_back();
}

MrefMatch AnnotationProcessor::searchMrefHitsNew(const BreakpointReduced& bpIn, int distanceThreshold, int conservativeDistanceThreshold, std::vector<std::vector<MrefEntryAnno> >& mref) {
	auto convertedChrIndex = ChrConverter::indexConverter[bpIn.getChrIndex()];
	std::vector<SuppAlignmentAnno> suppMatches { };
	if (convertedChrIndex < 0) {
		return MrefMatch { 0, 0, 10000, suppMatches };
	}
	auto itStart = std::lower_bound(mref[convertedChrIndex].begin(), mref[convertedChrIndex].end(), bpIn);
	if (itStart == mref[convertedChrIndex].end()) {
		return MrefMatch { 0, 0, 10000, suppMatches };
	}
	if (itStart != mref[convertedChrIndex].begin() && !(itStart->distanceTo(bpIn) < SvEvent::GERMLINEOFFSETTHRESHOLD) && std::prev(itStart)->distanceTo(bpIn) < SvEvent::GERMLINEOFFSETTHRESHOLD) {
		--itStart;
	}
	auto it = itStart;

	std::vector<std::vector<MrefEntryAnno>::iterator> dbHits { };
	std::vector<std::vector<MrefEntryAnno>::iterator> dbHitsConservative { };
	while (true) {
		auto tmpDistance = it->distanceTo(bpIn);
		if (tmpDistance < SvEvent::GERMLINEOFFSETTHRESHOLD) {
			dbHitsConservative.push_back(it);
		}
		if (tmpDistance < distanceThreshold) {
			dbHits.push_back(it);
		} else {
			break;
		}
		if (it == mref[convertedChrIndex].begin()) {
			break;
		}
		--it;
	}
	if (itStart != mref[convertedChrIndex].end()) {
		it = std::next(itStart);
		while (true) {
			if (it == mref[convertedChrIndex].end()) {
				break;
			}
			auto tmpDistance = it->distanceTo(bpIn);
			if (tmpDistance < SvEvent::GERMLINEOFFSETTHRESHOLD) {
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
		return MrefMatch { 0, 0, 0, suppMatches };
	}
	short score { 0 };
	auto offset = 0;
	for (auto res : dbHits) {
		auto saMatch = false;
		for (const auto &saRef : res->getSuppAlignments()) {
			for (const auto &sa : bpIn.getSuppAlignments()) {
				if (saRef.saCloseness(sa, SuppAlignmentAnno::DEFAULTREADLENGTH / 2)) {
					saMatch = true;
					auto previouslyRecorded = false;
					for (auto &saTmp : suppMatches) {
						if (saRef.saCloseness(saTmp, SuppAlignmentAnno::DEFAULTREADLENGTH / 2)) {
							previouslyRecorded = true;
							if (saTmp.isFuzzy()) {
								if (saRef.isFuzzy()) {
									saTmp.extendSuppAlignment(saRef.getPos(), saRef.getExtendedPos());
								} else {
									saTmp.removeFuzziness(saRef);
								}
							}
							if (saRef.getSupport() > saTmp.getSupport()) {
								saTmp.setSupport(saRef.getSupport());
							}
							if (saRef.getSecondarySupport() > saTmp.getSecondarySupport()) {
								saTmp.setSecondarySupport(saRef.getSecondarySupport());
							}
							if (saTmp.getSecondarySupport() < saTmp.getSupport()) {
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
	short conservativeScore { 0 };
	for (const auto res : dbHitsConservative) {
		auto tmpScore = res->getNumHits();
		if (tmpScore < SvEvent::RELAXEDBPFREQTHRESHOLD) {
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
	return MrefMatch { std::max(score, conservativeScore), conservativeScore, offset, suppMatches };
}

GermlineMatch AnnotationProcessor::searchGermlineHitsNew(const BreakpointReduced& bpIn, int distanceThreshold, int conservativeDistanceThreshold) {
	GermlineMatch dummyMatchTrue { 1.0, 1.0, std::vector<std::pair<SuppAlignmentAnno, double>> { } };
	GermlineMatch dummyMatchFalse { 0.0, 0.0, std::vector<std::pair<SuppAlignmentAnno, double>> { } };
	if (NOCONTROLMODE) {
		return dummyMatchTrue;
	}
	auto convertedChrIndex = ChrConverter::indexConverter[bpIn.getChrIndex()];
	if (convertedChrIndex < 0) {
		return dummyMatchFalse;
	}
	if (controlResults[convertedChrIndex].empty()) {
		return dummyMatchFalse;
	}
	auto itStart = std::lower_bound(controlResults[convertedChrIndex].begin(), controlResults[convertedChrIndex].end(), bpIn);
	if (itStart == controlResults[convertedChrIndex].end()) {
		return dummyMatchFalse;
	}
	if (itStart != controlResults[convertedChrIndex].cbegin() && !(itStart->distanceToBp(bpIn) < SvEvent::GERMLINEOFFSETTHRESHOLD) && std::prev(itStart)->distanceToBp(bpIn) < SvEvent::GERMLINEOFFSETTHRESHOLD) {
		--itStart;
	}
	auto it = itStart;
	std::vector<std::vector<BreakpointReduced>::iterator> dbHits { };
	std::vector<std::vector<BreakpointReduced>::iterator> dbHitsConservative { };

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
		if (it == controlResults[convertedChrIndex].begin()) {
			break;
		}
		--it;
	}
	if (itStart != controlResults[convertedChrIndex].end()) {
		it = std::next(itStart);
		while (true) {
			if (it == controlResults[convertedChrIndex].end()) {
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
			for (const auto sa : res->getSuppAlignments()) {
				if (sa.getMateSupport() > mateSupSa) {
					mateSupSa = sa.getMateSupport();
				}
			}
			if (res->getMateSupport() > mateSupSa) {
				mateSupSa = res->getMateSupport();
			}
			auto support = 0.0 + res->getPairedBreaksSoft() + res->getBreaksShortIndel() + res->getPairedBreaksHard() + res->getUnpairedBreaksSoft() + res->getUnpairedBreaksHard() + mateSupSa;
			if (support > maxSupport) {
				maxSupport = support;
				if (support + res->getNormalSpans() > 0) {
					conservativeClonality = support / (support + res->getNormalSpans());
				}
			}
		}
	}
	auto clonality = conservativeClonality;
	std::vector<std::pair<SuppAlignmentAnno, double>> suppMatches { };
	for (auto res : dbHits) {
		auto breakSupportSoft = res->getPairedBreaksSoft() + res->getUnpairedBreaksSoft();
		auto breakSupportHard = res->getPairedBreaksHard() + res->getUnpairedBreaksHard();
		for (const auto &saRef : res->getSuppAlignments()) {
			for (const auto &sa : bpIn.getSuppAlignments()) {
				if (saRef.saCloseness(sa, SuppAlignmentAnno::DEFAULTREADLENGTH / 2)) {
					auto previouslyRecorded = false;
					for (auto &saTmp : suppMatches) {
						if (saRef.saClosenessDirectional(saTmp.first, SuppAlignmentAnno::DEFAULTREADLENGTH / 2)) {
							previouslyRecorded = true;
							if (saTmp.first.isFuzzy()) {
								if (saRef.isFuzzy()) {
									saTmp.first.extendSuppAlignment(saRef.getPos(), saRef.getExtendedPos());
								} else {
									saTmp.first.removeFuzziness(saRef);
								}
							}
							if (saRef.getSupport() > saTmp.first.getSupport()) {
								saTmp.first.setSupport(saRef.getSupport());
							}
							if (saRef.getSecondarySupport() > saTmp.first.getSecondarySupport()) {
								saTmp.first.setSecondarySupport(saRef.getSecondarySupport());
							}
							if (saRef.getMateSupport() > saTmp.first.getMateSupport()) {
								saTmp.first.setMateSupport(saRef.getMateSupport());
								saTmp.first.setExpectedDiscordants(saRef.getExpectedDiscordants());
							}
							auto currentSoftSupport = std::max(saTmp.first.getSupport(), breakSupportSoft);
							auto currentHardSupport = std::max(saTmp.first.getSecondarySupport(), breakSupportHard);
							auto breakSupport = currentSoftSupport + currentHardSupport + 0.0;
							auto currentClonality = (breakSupport + saTmp.first.getMateSupport()) / (breakSupport + saTmp.first.getMateSupport() + res->getNormalSpans());
							saTmp.second = std::max(currentClonality, saTmp.second);
						}
					}
					if (!previouslyRecorded) {
						auto currentSoftSupport = std::max(saRef.getSupport(), breakSupportSoft);
						auto currentHardSupport = std::max(saRef.getSecondarySupport(), breakSupportHard);
						auto breakSupport = currentSoftSupport + currentHardSupport + 0.0;
						auto currentClonality = (breakSupport + saRef.getMateSupport()) / (breakSupport + saRef.getMateSupport() + res->getNormalSpans());
						suppMatches.push_back( { saRef, currentClonality });
					}
				}
			}
		}
	}
	return GermlineMatch { clonality, conservativeClonality, suppMatches };
}

bool AnnotationProcessor::applyMassiveInversionFiltering(bool stricterMode, bool controlCheckMode) {
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
				if (sv.getSelectedSa1().isProperPairErrorProne() && sv.getSelectedSa1().getMateSupport() < 6) {
					sv.setToRemove(true);
					continue;
				}
				if (sv.getSelectedSa2().isProperPairErrorProne() && sv.getSelectedSa2().getMateSupport() < 6) {
					sv.setToRemove(true);
					continue;
				}
				if (sv.isSemiSuspicious()) {
					if (sv.getEventScore() < 3) {
						sv.setToRemove(true);
						continue;
					}
					if (sv.getSelectedSa1().isSemiSuspicious() && sv.getSelectedSa2().isSemiSuspicious()) {
						sv.setToRemove(true);
						continue;
					}
					if (sv.getEvidenceLevel2() < 3 && sv.getSelectedSa1().isSemiSuspicious()) {
						if (sv.getEvidenceLevel1() < 3 || sv.getSelectedSa1().getSupport() < 5) {
							sv.setToRemove(true);
							continue;
						}
					}
					if (sv.getEvidenceLevel1() < 3 && sv.getSelectedSa2().isSemiSuspicious()) {
						if (sv.getEvidenceLevel2() < 3 || sv.getSelectedSa2().getSupport() < 5) {
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
					if (sv.getSelectedSa1().isStrictFuzzyCandidate() || sv.getSelectedSa2().isStrictFuzzyCandidate()) {
						sv.setToRemove(true);
						continue;
					}
					if (sv.isSemiSuspicious()) {
						sv.setToRemove(true);
						continue;
					}
				} else {
					if (sv.getSelectedSa1().isStrictFuzzyCandidate() && sv.getSelectedSa2().isStrictFuzzyCandidate()) {
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
					if (stricterMode || sv.getTotalEvidence1() < 5 || sv.getTotalEvidence2() < 5) {
						sv.setToRemove(true);
						continue;
					}
				} else {
					if (sv.getTotalEvidence1() < 5 && sv.getTotalEvidence2() < 5) {
						sv.setToRemove(true);
						continue;
					}
				}
			}
			if (stricterMode) {
				if (sv.getEvidenceLevel2() == 0) {
					if (sv.getSelectedSa1().getSupport() < 10 || !sv.isOverhang1Compensation()) {
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
					if (sv.getSelectedSa1().isProperPairErrorProne() && (sv.getSelectedSa2().isProperPairErrorProne() || sv.getEvidenceLevel2() == 0)) {
						if (sv.getSelectedSa1().getSupport() < 5 || sv.getSelectedSa2().getSupport() < 5) {
							sv.setToRemove(true);
							continue;
						}
						if (sv.getSelectedSa1().getSecondarySupport() < 3 && sv.getSelectedSa2().getSecondarySupport() < 3) {
							sv.setToRemove(true);
							continue;
						}
					}
					if (sv.getSelectedSa1().isSemiSuspicious()) {
						if (sv.getEvidenceLevel1() < 3 || sv.getSelectedSa1().getSupport() < 5) {
							sv.setToRemove(true);
							continue;
						}
					}
					if (sv.getSelectedSa2().isSemiSuspicious()) {
						if (sv.getEvidenceLevel2() < 3 || sv.getSelectedSa2().getSupport() < 5) {
							sv.setToRemove(true);
							continue;
						}
					}

					if (sv.getEvidenceLevel2() == 0) {
						if (sv.getEvidenceLevel1() < 3 || sv.getSelectedSa1().getSupport() < 5) {
							sv.setToRemove(true);
							continue;
						}
					}
					if (sv.getEvidenceLevel1() < 3 && sv.getEvidenceLevel2() < 3) {
						sv.setToRemove(true);
						continue;
					}
				} else {
					if (sv.getEvidenceLevel1() < 2 || sv.getEvidenceLevel2() < 2) {
						sv.setToRemove(true);
						continue;
					}
					if (sv.getTotalEvidence1() < 5 && sv.getTotalEvidence2() < 5) {
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
bool AnnotationProcessor::applyPathogenContaminationFiltering() {
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
				if (sv.getSelectedSa1().getSupport() > 19 && sv.getSelectedSa1().getSecondarySupport() < 3) {
					sv.setEventScore(1);
					sv.setEventType(5);
					continue;
				}
				if (sv.getSelectedSa2().getSupport() > 19 && sv.getSelectedSa2().getSecondarySupport() < 3) {
					sv.setEventScore(1);
					sv.setEventType(5);
					continue;
				}
				if (sv.getOverhang1lengthRatio() > 0.7 || sv.getOverhang2lengthRatio() > 0.7) {
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
void AnnotationProcessor::printFilteredResults(bool contaminationInControl, int controlPrefilteringLevel) const {
	if (controlPrefilteringLevel > 0) {
		std::cout << "#controlMassiveInvPrefilteringLevel\t" << controlPrefilteringLevel << std::endl;
	}
	if (massiveInvFilteringLevel > 0) {
		std::cout << "#tumorMassiveInvFilteringLevel\t" << massiveInvFilteringLevel << std::endl;
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

void AnnotationProcessor::printUnresolvedRareOverhangs(std::vector<std::vector<MrefEntryAnno>> &mref) {
	if (massiveInvFilteringLevel != 0) {
		return;
	}
	std::sort(visitedLineIndices.begin(), visitedLineIndices.end());
	std::unordered_set<int> visitedLineIndicesSet { visitedLineIndices.begin(), visitedLineIndices.end() };
	for (const auto &tumorChromosome : tumorResults) {
		for (const auto &bp : tumorChromosome) {
			if (visitedLineIndicesSet.count(bp.getLineIndex())) {
				continue;
			}
			if (bp.testOverhangBasedCandidacy()) {
				auto mrefHits = searchMrefHitsNew(bp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD, mref);
				if (mrefHits.getNumHits() < SvEvent::GERMLINEDBLIMIT) {
					auto germlineClonality = 1.0;
					if (!NOCONTROLMODE) {
						germlineClonality = searchGermlineHitsNew(bp, SuppAlignmentAnno::DEFAULTREADLENGTH * 6, SvEvent::GERMLINEOFFSETTHRESHOLD).getClonality();
					}
					std::string overhang { "" };
					{
						std::pair<int, std::string> dummy { bp.getLineIndex(), "" };
						auto lower = std::lower_bound(overhangs.cbegin(), overhangs.cend(), dummy);
						if (lower != overhangs.cend()) {
							if (lower->first == bp.getLineIndex()) {
								overhang = lower->second;
							} else if (std::next(lower) != overhangs.cend() && std::next(lower)->first == bp.getLineIndex()) {
								overhang = std::next(lower)->second;
							} else if (lower != overhangs.cbegin() && std::prev(lower)->first == bp.getLineIndex()) {
								overhang = std::prev(lower)->second;
							}
						}
					}
					if (!overhang.empty()) {
						std::cout << bp.printOverhang(germlineClonality, mrefHits.getNumHits(), overhang);
					}
				}
			}
		}
	}
}

} /* namespace sophia */
