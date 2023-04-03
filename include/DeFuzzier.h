/*
 * FuzzyCompression.h
 *
 *  Created on: 24 Nov 2016
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

#ifndef DEFUZZIER_H_
#define DEFUZZIER_H_

#include <BreakpointReduced.h>
#include <MrefEntry.h>
#include "SuppAlignment.h"
#include "SuppAlignmentAnno.h"
#include "Breakpoint.h"
#include <vector>
#include <map>
#include <unordered_set>

namespace sophia {

    using namespace std;
    
class DeFuzzier {
public:
	DeFuzzier(int maxDistanceIn, bool mrefModeIn);
	void deFuzzyDb(vector<BreakpointReduced>& bps) const;
	void deFuzzyDb(vector<MrefEntry>& bps) const;
private:
	void processFuzzySa(vector<BreakpointReduced>& bps, vector<BreakpointReduced>::iterator startingIt, SuppAlignmentAnno* startingSa) const;
	void dbSweep(vector<BreakpointReduced>& bps, vector<BreakpointReduced>::iterator startingIt, int increment, SuppAlignmentAnno* consensusSa, vector<SuppAlignmentAnno*>& processedSas) const;
	void selectBestSa(vector<SuppAlignmentAnno*>& processedSas, SuppAlignmentAnno* consensusSa) const;

	void processFuzzySa(vector<MrefEntry>& bps, vector<MrefEntry>::iterator startingIt, SuppAlignmentAnno* startingSa) const;
	void dbSweep(vector<MrefEntry>& bps, vector<MrefEntry>::iterator startingIt, unordered_set<short>& fileIndices, int increment, SuppAlignmentAnno* consensusSa, vector<SuppAlignmentAnno*>& processedSas) const;
	void selectBestSa(vector<SuppAlignmentAnno*>& processedSas, SuppAlignmentAnno* consensusSa, const unordered_set<short>& fileIndices) const;


	const int MAXDISTANCE;
	const bool MREFMODE;
};

}

#endif /* DEFUZZIER_H_ */
