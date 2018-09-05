/*
 * MrefMatch.cpp
 *
 *  Created on: 13 Jan 2018
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

#include <MrefMatch.h>
#include <iostream>

namespace sophia {

MrefMatch::MrefMatch(short numHitsIn, short numConsevativeHitsIn, int offsetDistanceIn, const std::vector<SuppAlignmentAnno>& suppMatchesIn) :
				numHits { numHitsIn },
				numConsevativeHits { numConsevativeHitsIn },
				offsetDistance { offsetDistanceIn },
				suppMatches { suppMatchesIn } {
}

}
/* namespace sophia */
