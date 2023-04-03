/*
 * MrefEntry.h
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

#ifndef MREFENTRY_H_
#define MREFENTRY_H_

#include <string>
#include <boost/format.hpp>
#include "SuppAlignment.h"
#include "BreakpointReduced.h"
namespace sophia {

    using namespace std;

class MrefEntry {
public:
	static int NUMPIDS;
	static int DEFAULTREADLENGTH;
	static boost::format doubleFormatter;
	MrefEntry();
	void addEntry(Breakpoint& tmpBreakpoint, int fileIndex);
	void addEntry(BreakpointReduced& tmpBreakpoint, int fileIndex);
	void mergeMrefEntries(MrefEntry &entry2);

	int getPos() const {
		return pos;
	}

	const vector<float>& getArtifactRatios() const {
		return artifactRatios;
	}

	const vector<short>& getFileIndices() const {
		return fileIndices;
	}

	short getValidityScore() const {
		return validity;
	}
	void removeMarkedFuzzies() {
		suppAlignments.erase(remove_if(suppAlignments.begin(), suppAlignments.end(), [](const SuppAlignmentAnno& sa) {return sa.isToRemove();}), suppAlignments.end());
	}
	string printBpInfo(const string& chromosome);
	string printArtifactRatios(const string& chromosome);
	SuppAlignmentAnno* searchFuzzySa(const SuppAlignmentAnno& fuzzySa);
	vector<SuppAlignmentAnno*> getSupplementsPtr() {
		vector<SuppAlignmentAnno*> res { };
		for (auto &sa : suppAlignments) {
			res.push_back(&sa);
		}
		return res;
	}
	const vector<short>& getFileIndicesWithArtifactRatios() const {
		return fileIndicesWithArtifactRatios;
	}
	const vector<SuppAlignmentAnno>& getSuppAlignments() const {
		return suppAlignments;
	}

	void setAsInvalid() {
		pos = -1;
		validity = -1;
	}

private:
	bool saMatcher(SuppAlignmentAnno* saPtr);
	void finalizeFileIndices();
	short validity; //-1 nothing, 0 only sa, 1 sa and support
	int pos;
	vector<short> fileIndices;
	vector<short> fileIndicesWithArtifactRatios;
	vector<float> artifactRatios;
	vector<SuppAlignmentAnno> suppAlignments;
};

} /* namespace sophia */

#endif /* MREFENTRY_H_ */
