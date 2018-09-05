/*
 * CoverageAtBase.h
 *
 *  Created on: 18 Apr 2016
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

#ifndef COVERAGEATBASE_H_
#define COVERAGEATBASE_H_
#include <memory>
namespace sophia {

class CoverageAtBase {
public:
	CoverageAtBase() :
					coverage { 0 },
					normalBpsSoft { 0 },
					normalBpsHard { 0 },
					normalBpsShortIndel { 0 },
					normalSpans { 0 },
					lowQualSpansSoft { 0 },
					lowQualSpansHard { 0 },
					lowQualBpsSoft { 0 },
					lowQualBpsHard { 0 },
					decoyBlacklisted { false } {
	}
	~CoverageAtBase() = default;
	int getCoverage() const {
		return coverage;
	}
	int getLowQualBpsSoft() const {
		return lowQualBpsSoft;
	}
	int getLowQualBpsHard() const {
		return lowQualBpsHard;
	}
	int getLowQualSpansSoft() const {
		return lowQualSpansSoft;
	}
	int getLowQualSpansHard() const {
		return lowQualSpansHard;
	}
	int getNormalBpsHard() const {
		return normalBpsHard;
	}
	int getNormalBpsSoft() const {
		return normalBpsSoft;
	}
	int getNormalSpans() const {
		return normalSpans;
	}
	void incrementCoverage() {
		++coverage;
	}
	void incrementLowQualBpsSoft() {
		++lowQualBpsSoft;
	}
	void incrementLowQualBpsHard() {
		++lowQualBpsHard;
	}
	void incrementLowQualSpansSoft() {
		++lowQualSpansSoft;
	}
	void incrementLowQualSpansHard() {
		++lowQualSpansHard;
	}
	void incrementNormalBpsHard() {
		++normalBpsHard;
	}
	void incrementNormalBpsSoft() {
		++normalBpsSoft;
	}
	void incrementNormalBpsShortIndel() {
		++normalBpsShortIndel;
	}
	void incrementNormalSpans() {
		++normalSpans;
	}
	void decrementLowQualSpansHard() {
		--lowQualSpansHard;
	}
	void decrementLowQualSpansSoft() {
		--lowQualSpansSoft;
	}
	void decrementNormalSpans() {
		--normalSpans;
	}
	int getNormalBpsShortIndel() const {
		return normalBpsShortIndel;
	}
	bool isDecoyBlacklisted() const {
		return decoyBlacklisted;
	}
	void setDecoyBlacklisted(bool decoyBlacklisted) {
		this->decoyBlacklisted = decoyBlacklisted;
	}

private:
	int coverage;
	int normalBpsSoft;
	int normalBpsHard;
	int normalBpsShortIndel;
	int normalSpans;
	int lowQualSpansSoft;
	int lowQualSpansHard;
	int lowQualBpsSoft;
	int lowQualBpsHard;
	bool decoyBlacklisted;
};

} /* namespace sophia */

#endif /* COVERAGEATBASE_H_ */
