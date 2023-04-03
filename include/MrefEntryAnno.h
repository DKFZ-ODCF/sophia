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

#ifndef MREFENTRYANNO_H_
#define MREFENTRYANNO_H_

#include <string>
#include <boost/format.hpp>
#include <BreakpointReduced.h>
#include "SuppAlignmentAnno.h"

namespace sophia {

    using namespace std;

class MrefEntryAnno {

public:
	static int PIDSINMREF;
	static int DEFAULTREADLENGTH;
	static boost::format doubleFormatter;
	MrefEntryAnno(const string &mrefEntryIn);
	template<typename T>
	bool operator<(const T& rhs) const {
		return pos < rhs.getPos();
	}
	template<typename T>
	int distanceTo(const T& rhs) const {
		return abs(pos - rhs.getPos());
	}
	template<typename T>
	int distanceToBp(const T &compIn) const {
		return abs(pos - compIn.getPos());
	}
	bool operator==(const MrefEntryAnno& rhs) const {
		return pos == rhs.getPos();
	}

	int getPos() const {
		return pos;
	}
	vector<SuppAlignmentAnno*> getSuppAlignmentsPtr() {
		vector<SuppAlignmentAnno*> res { };
		for (auto &sa : suppAlignments) {
			res.push_back(&sa);
		}
		return res;
	}

	void removeMarkedFuzzies() {
		while (!suppAlignments.empty() && suppAlignments.back().isToRemove()) {
			suppAlignments.pop_back();
		}
		for (auto saIt = suppAlignments.begin(); saIt != suppAlignments.end(); ++saIt) {
			if (saIt->isToRemove()) {
				swap(*saIt, suppAlignments.back());
			}
			while (!suppAlignments.empty() && suppAlignments.back().isToRemove()) {
				suppAlignments.pop_back();
			}
		}
	}
//	SuppAlignmentAnno* searchFuzzySa(const SuppAlignmentAnno& fuzzySa);

	const vector<SuppAlignmentAnno>& getSuppAlignments() const {
		return suppAlignments;
	}

	vector<SuppAlignmentAnno*> getSupplementsPtr() {
		vector<SuppAlignmentAnno*> res { };
		for (auto &sa : suppAlignments) {
			res.push_back(&sa);
		}
		return res;
	}
	bool closeToSupp(const SuppAlignmentAnno &compIn, int fuzziness) const {
		if (compIn.isFuzzy()) {
			fuzziness = 2.5 * DEFAULTREADLENGTH;
			return (pos - fuzziness) <= (compIn.getExtendedPos() + fuzziness) && (compIn.getPos() - fuzziness) <= (pos + fuzziness);
		} else {
			return abs(pos - compIn.getPos()) <= fuzziness;
		}
	}
	int distanceToSupp(const SuppAlignmentAnno &compIn) const {
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
			return abs(pos - compIn.getPos());
		}
	}
	short getNumHits() const {
		return numHits;
	}

	void setNumHits(short numHits) {
		this->numHits = numHits;
	}
private:
	int pos;
	short numHits;
	vector<SuppAlignmentAnno> suppAlignments;
};

} /* namespace sophia */

#endif /* MREFENTRYANNO_H_ */
