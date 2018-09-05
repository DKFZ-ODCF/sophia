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

#ifndef BREAKPOINTREDUCED_H_
#define BREAKPOINTREDUCED_H_

#include <string>
#include <boost/format.hpp>
#include "SuppAlignmentAnno.h"
#include "Breakpoint.h"
#include "MrefMatch.h"
#include "GermlineMatch.h"
#include <iostream>
namespace sophia {
class BreakpointReduced {
public:
	static int DEFAULTREADLENGTH;
	static double CLONALITYSTRICTLOWTHRESHOLD;
	static double ARTIFACTFREQHIGHTHRESHOLD;
	static std::string PIDSINMREFSTR;
	static boost::format doubleFormatter;
	BreakpointReduced(const Breakpoint& tmpBp, int lineIndexIn, bool hasOverhangIn);
	BreakpointReduced(const SuppAlignmentAnno& sa, const BreakpointReduced& emittingBp, bool fuzzySecondary);

	template<typename T>
	bool operator<(const T& rhs) const {
		return pos < rhs.getPos();
	}
	bool fullSmaller(const BreakpointReduced& rhs) const {
		if (chrIndex < rhs.getChrIndex()) {
			return true;
		}
		if (chrIndex > rhs.getChrIndex()) {
			return false;
		}
		return pos < rhs.getPos();
	}
	template<typename T>
	int distanceTo(const T& rhs) const {
		if (chrIndex != rhs.getChrIndex()) {
			return 1000000;
		} else {
			return std::abs(pos - rhs.getPos());
		}
	}
	template<typename T>
	int distanceToBp(const T &compIn) const {
		if (chrIndex == compIn.getChrIndex()) {
			return std::abs(pos - compIn.getPos());
		} else {
			return -1;
		}
	}
	int getChrIndex() const {
		return chrIndex;
	}

	int getPos() const {
		return pos;
	}

	bool isToRemove() const {
		return toRemove;
	}
	void setToRemove(bool toRemove) {
		this->toRemove = toRemove;
	}
	void removeMarkedFuzzies();
	SuppAlignmentAnno* searchFuzzySa(const SuppAlignmentAnno& fuzzySa);

	int getBreaksShortIndel() const {
		return breaksShortIndel;
	}

	int getLeftCoverage() const {
		return leftCoverage;
	}

	int getLowQualBreaksHard() const {
		return lowQualBreaksHard;
	}

	int getLowQualBreaksSoft() const {
		return lowQualBreaksSoft;
	}

	int getLowQualSpansHard() const {
		return lowQualSpansHard;
	}

	int getLowQualSpansSoft() const {
		return lowQualSpansSoft;
	}

	int getMateSupport() const {
		return mateSupport;
	}
	int getNormalSpans() const {
		return normalSpans;
	}
	int getPairedBreaksHard() const {
		return pairedBreaksHard;
	}
	int getPairedBreaksSoft() const {
		return pairedBreaksSoft;
	}
	int getRepetitiveOverhangBreaks() const {
		return repetitiveOverhangBreaks;
	}
	int getRightCoverage() const {
		return rightCoverage;
	}
	const std::vector<SuppAlignmentAnno>& getSuppAlignments() const {
		return suppAlignments;
	}

	int getUnpairedBreaksHard() const {
		return unpairedBreaksHard;
	}
	int getUnpairedBreaksSoft() const {
		return unpairedBreaksSoft;
	}
	int getLineIndex() const {
		return lineIndex;
	}
	std::vector<SuppAlignmentAnno*> getSupplementsPtr() {
		std::vector<SuppAlignmentAnno*> res { };
		for (auto &sa : suppAlignments) {
			res.push_back(&sa);
		}
		return res;
	}
	bool closeToSupp(const SuppAlignmentAnno &compIn, int fuzziness) const {
		if (chrIndex == compIn.getChrIndex()) {
			if (compIn.isFuzzy()) {
				fuzziness = 2.5 * DEFAULTREADLENGTH;
				return (pos - fuzziness) <= (compIn.getExtendedPos() + fuzziness) && (compIn.getPos() - fuzziness) <= (pos + fuzziness);
			} else {
				return std::abs(pos - compIn.getPos()) <= fuzziness;
			}
		} else {
			return false;
		}
	}
	int distanceToSupp(const SuppAlignmentAnno &compIn) const {
		if (chrIndex == compIn.getChrIndex()) {
			if (compIn.isFuzzy()) {
				if (compIn.getPos() <= pos && pos <= compIn.getExtendedPos()) {
					return 0;
				} else {
					if (pos < compIn.getPos()) {
						return compIn.getPos() - pos;
					} else {
						return pos - compIn.getExtendedPos();
					}
				}
			} else {
				return std::abs(pos - compIn.getPos());
			}
		} else {
			return 1000000;
		}
	}
	const MrefMatch& getMrefHits() const {
		return mrefHits;
	}
	void setMrefHits(MrefMatch mrefHits) {
		this->mrefHits = mrefHits;
	}
	void setGermlineInfo(GermlineMatch germlineInfo) {
		this->germlineInfo = germlineInfo;
	}

	bool testOverhangBasedCandidacy() const;
	std::string printOverhang(double germlineClonality, int numHits, const std::string &overhang) const;

	const GermlineMatch& getGermlineInfo() const {
		return germlineInfo;
	}
	void addFileIndex(int fileIndex) {
		for (auto &sa : suppAlignments) {
			sa.addFileIndex(fileIndex);
		}
	}
	void complexRearrangementMateRatioRescue(bool encounteredM);
	bool hasOverhang;
	void addDummySa(const SuppAlignmentAnno& sa, const BreakpointReduced& emittingBp);
	const SuppAlignmentAnno& getDummySa();
private:
	bool toRemove;
	int lineIndex;
	int chrIndex;
	int pos;
	int normalSpans, lowQualSpansSoft, lowQualSpansHard, unpairedBreaksSoft, unpairedBreaksHard, breaksShortIndel, lowQualBreaksSoft, lowQualBreaksHard, repetitiveOverhangBreaks;
	int pairedBreaksSoft, pairedBreaksHard;
	int mateSupport;
	int leftCoverage, rightCoverage;
	MrefMatch mrefHits;
	GermlineMatch germlineInfo;
	std::vector<SuppAlignmentAnno> suppAlignments;
};

} /* namespace sophia */

#endif /* BREAKPOINTREDUCED_H_ */
