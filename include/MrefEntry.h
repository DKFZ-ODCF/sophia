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

	const std::vector<float>& getArtifactRatios() const {
		return artifactRatios;
	}

	const std::vector<short>& getFileIndices() const {
		return fileIndices;
	}

	short getValidityScore() const {
		return validity;
	}
	void removeMarkedFuzzies() {
		suppAlignments.erase(std::remove_if(suppAlignments.begin(), suppAlignments.end(), [](const SuppAlignmentAnno& sa) {return sa.isToRemove();}), suppAlignments.end());
	}
	std::string printBpInfo(const std::string& chromosome);
	std::string printArtifactRatios(const std::string& chromosome);
	SuppAlignmentAnno* searchFuzzySa(const SuppAlignmentAnno& fuzzySa);
	std::vector<SuppAlignmentAnno*> getSupplementsPtr() {
		std::vector<SuppAlignmentAnno*> res { };
		for (auto &sa : suppAlignments) {
			res.push_back(&sa);
		}
		return res;
	}
	const std::vector<short>& getFileIndicesWithArtifactRatios() const {
		return fileIndicesWithArtifactRatios;
	}
	const std::vector<SuppAlignmentAnno>& getSuppAlignments() const {
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
	std::vector<short> fileIndices;
	std::vector<short> fileIndicesWithArtifactRatios;
	std::vector<float> artifactRatios;
	std::vector<SuppAlignmentAnno> suppAlignments;
};

} /* namespace sophia */

#endif /* MREFENTRY_H_ */
