/*
 * Sdust.h
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

#ifndef SDUST_H_
#define SDUST_H_
#include <vector>
#include <deque>
#include <set>
#include <algorithm>
namespace sophia {
struct PerfectInterval;
class Sdust {
public:
	Sdust(const std::vector<int>& overhangIn);
	~Sdust() = default;
	const std::vector<PerfectInterval>& getRes() const {
		return res;
	}
private:
	static const int SCORETHRESHOLD = 20;
	static const int WINDOWSIZE = 64;
	std::vector<PerfectInterval> res;
	std::set<PerfectInterval> P;
	std::deque<int> w;
	int L;
	int rW;
	int rV;
	std::vector<int> cW;
	std::vector<int> cV;
	void saveMaskedRegions(int wStart);
	int triplet(const std::vector<int>&overhangIn, int indexPos);
	void shiftWindow(int t);
	void addTripletInfo(int& r, std::vector<int>& c, int t);
	void removeTripletInfo(int& r, std::vector<int>& c, int t);
	void findPerfectRegions(int wStart, int r, std::vector<int> c);
};
struct PerfectInterval {
	int startIndex;
	int endIndex;
	double score;
	bool operator<(const PerfectInterval & rhs) const {
		if (startIndex > rhs.startIndex) return true;
		if (startIndex < rhs.startIndex) return false;
		if (endIndex > rhs.endIndex) return false;
		if (endIndex < rhs.endIndex) return true;
		return false;
	}
};
} /* namespace sophia */

#endif /* SDUST_H_ */
