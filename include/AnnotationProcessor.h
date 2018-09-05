/*
 * AnnotationProcessor.h
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

#ifndef ANNOTATIONPROCESSOR_H_
#define ANNOTATIONPROCESSOR_H_
#include <BreakpointReduced.h>
#include <MrefEntryAnno.h>
#include <SvEvent.h>
#include <string>
#include <vector>
#include <deque>
#include <utility>
#include <unordered_set>
#include "SuppAlignmentAnno.h"
#include "ChrConverter.h"
#include "MrefMatch.h"
#include "GermlineMatch.h"
#include <memory>
//
//struct VectorHash {
//	size_t operator()(const std::vector<int>& v) const {
//		std::hash<int> hasher;
//		size_t seed = 0;
//		for (int i : v) {
//			seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//		}
//		return seed;
//	}
//};

namespace sophia {
class AnnotationProcessor {
public:
	static bool ABRIDGEDOUTPUT;
	AnnotationProcessor(const std::string &tumorResultsIn, std::vector<std::vector<MrefEntryAnno>> &mref, int defaultReadLengthTumorIn, bool controlCheckModeIn, int germlineDbLimit);
	AnnotationProcessor(const std::string &tumorResultsIn, std::vector<std::vector<MrefEntryAnno>> &mref, const std::string &controlResultsIn, int defaultReadLengthTumorIn, int defaultReadLengthControlIn, int germlineDbLimit, int lowQualControlIn, bool pathogenInControlIn);
	void printFilteredResults(bool contaminationInControl, int controlPrefilteringLevel) const;
	int getMassiveInvFilteringLevel() const {
		return massiveInvFilteringLevel;
	}

	bool isContaminationObserved() const {
		return contaminationObserved;
	}

private:
	void searchMatches(std::vector<std::vector<MrefEntryAnno>>& mref);
	void createDoubleMatchSv(BreakpointReduced& sourceBp, BreakpointReduced& targetBp, const SuppAlignmentAnno& sa, const SuppAlignmentAnno& saMatch, bool checkOrder, std::vector<std::vector<MrefEntryAnno>>& mref);
	bool createDoubleMatchSvPreCheck(const SuppAlignmentAnno& saMatch);
	void createUnmatchedSaSv(BreakpointReduced& sourceBp, BreakpointReduced& targetBp, const SuppAlignmentAnno& sa, std::vector<std::vector<MrefEntryAnno>>& mref);
	void createUnknownMatchSv(BreakpointReduced& sourceBp, const SuppAlignmentAnno& sa, std::vector<std::vector<MrefEntryAnno>>& mref, bool doubleSupportSa);
	bool createUnknownMatchSvPreCheck(const SuppAlignmentAnno& sa, bool doubleSupportSa);
	void checkSvQuality();
	MrefMatch searchMrefHitsNew(const BreakpointReduced& bpIn, int distanceThreshold, int conservativeDistanceThreshold, std::vector<std::vector<MrefEntryAnno>> &mref);
	GermlineMatch searchGermlineHitsNew(const BreakpointReduced& bpIn, int distanceThreshold, int conservativeDistanceThreshold);

	void searchSa(int chrIndex, int dbIndex, const SuppAlignmentAnno& sa, bool doubleSupportSa, std::vector<std::vector<MrefEntryAnno>> &mref);
	bool applyMassiveInversionFiltering(bool stricterMode, bool controlCheckMode);
	bool applyPathogenContaminationFiltering();
	void printUnresolvedRareOverhangs(std::vector<std::vector<MrefEntryAnno>> &mref);
	const bool NOCONTROLMODE;
	const int GERMLINEDBLIMIT;
	bool contaminationObserved;
	int massiveInvFilteringLevel;
//	std::unordered_set<std::vector<int>, VectorHash> filteredResultKeys;
	std::unordered_set<std::string> filteredResultKeys;
	std::vector<SvEvent> filteredResults;
	std::vector<std::vector<BreakpointReduced>> tumorResults;
	std::vector<std::vector<BreakpointReduced>> controlResults;
	std::vector<std::pair<int, std::string>> overhangs;
	std::vector<int> visitedLineIndices;
};
} /* namespace sophia */

#endif /* ANNOTATIONPROCESSOR_H_ */
