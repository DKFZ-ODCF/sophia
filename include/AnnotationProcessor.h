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
//	size_t operator()(const vector<int>& v) const {
//		hash<int> hasher;
//		size_t seed = 0;
//		for (int i : v) {
//			seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//		}
//		return seed;
//	}
//};

namespace sophia {

    using namespace std;
    
    class AnnotationProcessor {
        public:
            static bool ABRIDGEDOUTPUT;
            AnnotationProcessor(const string &tumorResultsIn, vector<vector<MrefEntryAnno>> &mref, int defaultReadLengthTumorIn, bool controlCheckModeIn, int germlineDbLimit);
            AnnotationProcessor(const string &tumorResultsIn, vector<vector<MrefEntryAnno>> &mref, const string &controlResultsIn, int defaultReadLengthTumorIn, int defaultReadLengthControlIn, int germlineDbLimit, int lowQualControlIn, bool pathogenInControlIn);
            void printFilteredResults(bool contaminationInControl, int controlPrefilteringLevel) const;
            int getMassiveInvFilteringLevel() const {
                return massiveInvFilteringLevel;
            }
        
            bool isContaminationObserved() const {
                return contaminationObserved;
            }
        
        private:
            void searchMatches(vector<vector<MrefEntryAnno>>& mref);
            void createDoubleMatchSv(BreakpointReduced& sourceBp, BreakpointReduced& targetBp, const SuppAlignmentAnno& sa, const SuppAlignmentAnno& saMatch, bool checkOrder, vector<vector<MrefEntryAnno>>& mref);
            bool createDoubleMatchSvPreCheck(const SuppAlignmentAnno& saMatch);
            void createUnmatchedSaSv(BreakpointReduced& sourceBp, BreakpointReduced& targetBp, const SuppAlignmentAnno& sa, vector<vector<MrefEntryAnno>>& mref);
            void createUnknownMatchSv(BreakpointReduced& sourceBp, const SuppAlignmentAnno& sa, vector<vector<MrefEntryAnno>>& mref, bool doubleSupportSa);
            bool createUnknownMatchSvPreCheck(const SuppAlignmentAnno& sa, bool doubleSupportSa);
            void checkSvQuality();
            MrefMatch searchMrefHitsNew(const BreakpointReduced& bpIn, int distanceThreshold, int conservativeDistanceThreshold, vector<vector<MrefEntryAnno>> &mref);
            GermlineMatch searchGermlineHitsNew(const BreakpointReduced& bpIn, int distanceThreshold, int conservativeDistanceThreshold);
        
            void searchSa(int chrIndex, int dbIndex, const SuppAlignmentAnno& sa, bool doubleSupportSa, vector<vector<MrefEntryAnno>> &mref);
            bool applyMassiveInversionFiltering(bool stricterMode, bool controlCheckMode);
            bool applyPathogenContaminationFiltering();
            void printUnresolvedRareOverhangs(vector<vector<MrefEntryAnno>> &mref);
            const bool NOCONTROLMODE;
            const int GERMLINEDBLIMIT;
            bool contaminationObserved;
            int massiveInvFilteringLevel;
        //	unordered_set<vector<int>, VectorHash> filteredResultKeys;
            unordered_set<string> filteredResultKeys;
            vector<SvEvent> filteredResults;
            vector<vector<BreakpointReduced>> tumorResults;
            vector<vector<BreakpointReduced>> controlResults;
            vector<pair<int, string>> overhangs;
            vector<int> visitedLineIndices;
    };
} /* namespace sophia */

#endif /* ANNOTATIONPROCESSOR_H_ */
