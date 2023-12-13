/*
 *     Author: Philip R. Kensche, DKFZ Heidelberg (Omics IT and Data Management Core Facility)
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
 *     LICENSE: GPL
 */

#include "ChrConverter.h"
#include <stdexcept>
#include <vector>
#include <string>

namespace sophia {

    using namespace std;

    ChrConverter::ChrConverter(const vector<string> &indexToChr,
                               const vector<string> &indexToChrCompressedMref,
                               const vector<int> &chrSizesCompressedMref,
                               const vector<int> &indexConverter) :
                    indexToChr(indexToChr),
                    indexToChrCompressedMref(indexToChrCompressedMref),
                    chrSizesCompressedMref(chrSizesCompressedMref),
                    indexConverter(indexConverter) {
            if (indexToChr.size() != indexConverter.size())
                throw invalid_argument(
                    "indexToChr and indexConverter must have the same size");
            if (indexToChrCompressedMref.size() != chrSizesCompressedMref.size())
                throw invalid_argument(
                    "indexToChrCompressedMref and chrSizesCompressedMref must have the same size");
        }

    ChrConverter::~ChrConverter() {}

}