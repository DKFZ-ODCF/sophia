/*
 * Alignment.h
 *
 *  Created on: 16 Apr 2016
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

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_
#include <OverhangRange.h>
#include <vector>
#include <deque>
#include <string>
#include <algorithm>
#include <memory>
#include "CigarChunk.h"
#include "SuppAlignment.h"
#include "ChosenBp.h"
#include "CoverageAtBase.h"

namespace sophia {
class Alignment {
public:
	Alignment();
	void continueConstruction();
	static int LOWQUALCLIPTHRESHOLD, BASEQUALITYTHRESHOLD, BASEQUALITYTHRESHOLDLOW, CLIPPEDNUCLEOTIDECOUNTTHRESHOLD, INDELNUCLEOTIDECOUNTTHRESHOLD;
	static double ISIZEMAX;
	int getStartPos() const {
		return startPos;
	}
	int getEndPos() const {
		return endPos;
	}
	int getReadType() const {
		return readType;
	}
	const std::vector<int>& getReadBreakpoints() const {
		return readBreakpoints;
	}
	bool isValidLine() const {
		return validLine;
	}
	const std::string& getSamLine() const {
		return samLine;
	}
	const std::vector<int>& getSamChunkPositions() const {
		return samChunkPositions;
	}
	bool assessOutlierMateDistance();
	int getMateChrIndex() const {
		return mateChrIndex;
	}
	int getMatePos() const {
		return matePos;
	}
	const std::vector<char>& getReadBreakpointTypes() const {
		return readBreakpointTypes;
	}
	void setChosenBp(int chosenBpLoc, int alignmentIndex);
	bool isOverhangEncounteredM() const {
		return chosenBp->bpEncounteredM;
	}
	int getOverhangLength() const {
		return chosenBp->overhangLength;
	}
	int getOverhangStartIndex() const {
		return chosenBp->overhangStartIndex;
	}
	std::vector<SuppAlignment> generateSuppAlignments(int bpChrIndex, int bpPos);
	const std::vector<SuppAlignment>& getSupplementaryAlignments() const {
		return chosenBp->supplementaryAlignments;
	}
	int getChrIndex() const {
		return chrIndex;
	}
	const std::vector<int>& getReadBreakpointsSizes() const {
		return readBreakpointSizes;
	}
	bool isLowMapq() const {
		return lowMapq;
	}
	bool isNullMapq() const {
		return nullMapq;
	}
	bool isSupplementary() const {
		return supplementary;
	}
	void addChildNode(int indexIn) {
		chosenBp->addChildNode(indexIn);
	}
	void addSupplementaryAlignments(const std::vector<SuppAlignment>& suppAlignments) {
		chosenBp->addSupplementaryAlignments(suppAlignments);
	}
	const std::vector<int>& getChildrenNodes() const {
		return chosenBp->childrenNodes;
	}
	int getOriginIndex() const {
		return chosenBp->selfNodeIndex;
	}
	std::string printOverhang() const;
	double overhangComplexityMaskRatio() const;

	bool isInvertedMate() const {
		return invertedMate;
	}
	bool isDistantMate() const {
		return distantMate == 1;
	}

private:
	void mappingQualityCheck();
	bool isEventCandidate() const;
	void createCigarChunks();
	void assignBreakpointsAndOverhangs();
	void qualityCheckCascade();
	bool clipCountCheck();
	bool uniqueSuppCheck();
	double overhangMedianQuality(const CigarChunk& cigarChunk) const;
	template<typename Iterator>
	void fullMedianQuality(Iterator qualBegin, Iterator qualEnd, std::vector<int>& overhangPerBaseQuality) const;
	template<typename Iterator>
	double getMedian(Iterator begin, Iterator end) const;
	void assessReadType();
	bool lowMapq;
	bool nullMapq;
	int distantMate;
	std::unique_ptr<ChosenBp> chosenBp;
	int chrIndex;
	int readType;
	int startPos, endPos;
	int mateChrIndex, matePos;
	std::string samLine;
	bool validLine;
	std::vector<int> samChunkPositions;
	std::string::const_iterator saCbegin, saCend;
	bool hasSa;
	bool supplementary;
	bool fwdStrand;
	bool invertedMate;
	bool qualChecked;
	std::vector<CigarChunk> cigarChunks;
	std::vector<int> readBreakpoints;
	std::vector<char> readBreakpointTypes;
	std::vector<int> readBreakpointSizes;
	std::vector<double> readBreakpointComplexityMaskRatios;
	std::deque<bool> readBreakpointsEncounteredM;
	std::vector<OverhangRange> readOverhangCoords;
};

template<typename Iterator>
void Alignment::fullMedianQuality(Iterator qualBegin, Iterator qualEnd, std::vector<int>& overhangPerBaseQuality) const {
	overhangPerBaseQuality.reserve(std::distance(qualBegin, qualEnd));
	auto consecutiveLowQuals = 0;
	for (auto cit = qualBegin; cit != qualEnd; ++cit) {
		if (*cit < BASEQUALITYTHRESHOLDLOW) { //33 + phred 11
			if (consecutiveLowQuals == 5) {
				overhangPerBaseQuality.clear();
				return;
			}
			++consecutiveLowQuals;
		} else {
			consecutiveLowQuals = 0;
		}
		overhangPerBaseQuality.push_back(*cit);
	}
}

//Median Code taken from http://rosettacode.org/wiki/Averages/Median#C.2B.2B
template<typename Iterator>
double Alignment::getMedian(Iterator begin, Iterator end) const {
	// this is middle for odd-length, and "upper-middle" for even length
	Iterator middle = begin + (end - begin) / 2;
	// This function runs in O(n) on average, according to the standard
	std::nth_element(begin, middle, end);
	if ((end - begin) % 2 != 0) { // odd length
		return *middle;
	} else { // even length
		// the "lower middle" is the max of the lower half
		Iterator lower_middle = std::max_element(begin, middle);
		return (*middle + *lower_middle) / 2.0;
	}
}

} /* namespace sophia */

#endif /* ALIGNMENT_H_ */
