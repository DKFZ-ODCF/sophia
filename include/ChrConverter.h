/*
 * ChrConverter.h
 *
 *  Created on: 28 Dec 2017
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

#ifndef CHRCONVERTER_H_
#define CHRCONVERTER_H_
#include <string>
#include <array>

#include <iterator>
namespace sophia {

    using namespace std;

class ChrConverter {
public:
	static inline int readChromosomeIndex(string::const_iterator startIt, char stopChar) {
		int chrIndex { 0 };
		if (isdigit(*startIt)) {
			for (auto chr_cit = startIt; *chr_cit != stopChar; ++chr_cit) {
				chrIndex = chrIndex * 10 + (*chr_cit - '0');
			}
			return chrIndex;
		} else {
			switch (*startIt) {
			case 'h':
				return 999;
			case 'X':
				return 40;
			case 'G':
				for (auto cit = next(startIt, 2); *cit != '.'; ++cit) {
					chrIndex = 10 * chrIndex + *cit - '0';
				}
				return chrIndex;
			case 'Y':
				return 41;
			case 'M':
				++startIt;
				if (*startIt == 'T') {
					return 1001;
				} else {
					return 1003;
				}
			case 'N':
				return 1000;
			case 'p':
				return 1002;
			default:
				return 1003;
			}
		}
		return 0;
	}
	static const array<string, 1004> indexToChr;
	static const array<int, 1004> indexConverter;
	static const array<string, 85> indexToChrCompressedMref;

};

} /* namespace sophia */

#endif /* CHRCONVERTER_H_ */
