/*
 * MateInfo.h
 *
 *  Created on: 21 Apr 2016
 *      Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical
 * Bioinformatics, Bioinformatics and Omics Data Analytics and currently
 * Neuroblastoma Genomics) Copyright (C) 2018 Umut H. Toprak, Matthias
 * Schlesner, Roland Eils and DKFZ Heidelberg
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

#ifndef MATEINFO_H_
#define MATEINFO_H_
#include "SuppAlignment.h"
#include <cmath>

namespace sophia {

struct MateInfo {
    int readStartPos;
    int readEndPos;
    int mateChrIndex;
    int mateStartPos;
    int mateEndPos;
    bool inverted;
    int source;
    int evidenceLevel;
    int matePower;
    int inversionSupport;
    int straightSupport;
    std::vector<int> bpLocs;
    bool saSupporter;
    bool toRemove;
    bool operator<(const MateInfo &rhs) const {
        if (mateChrIndex < rhs.mateChrIndex)
            return true;
        if (mateChrIndex > rhs.mateChrIndex)
            return false;
        if (mateStartPos < rhs.mateStartPos)
            return true;
        return false;
    }
    bool suppAlignmentFuzzyMatch(const SuppAlignment &sa) const {
        if (mateChrIndex != sa.getChrIndex()) {
            return false;
        } else {
            if (!sa.isFuzzy()) {
                return sa.getPos() >= (mateStartPos - sa.getMatchFuzziness()) &&
                       sa.getPos() <= (mateEndPos + sa.getMatchFuzziness());
            } else {
                return (mateStartPos - sa.getMatchFuzziness()) <=
                           sa.getExtendedPos() &&
                       sa.getPos() <= (mateEndPos + sa.getMatchFuzziness());
            }
        }
    }
    MateInfo(int readStartPosIn, int readEndPosIn, int mateChrIndexIn,
             int mateStartPosIn, int sourceType, bool invertedIn)
        : readStartPos{readStartPosIn}, readEndPos{readEndPosIn},
          mateChrIndex{mateChrIndexIn}, mateStartPos{mateStartPosIn},
          mateEndPos{mateStartPosIn}, inverted{invertedIn}, source{sourceType},
          evidenceLevel{sourceType == 2 ? 3 : 1}, matePower{1},
          inversionSupport{invertedIn}, straightSupport{!invertedIn}, bpLocs{},
          saSupporter{false}, toRemove{false} {}
    MateInfo(int readStartPosIn, int readEndPosIn, int mateChrIndexIn,
             int mateStartPosIn, int sourceType, bool invertedIn,
             const std::vector<int> &bpLocsIn)
        : readStartPos{readStartPosIn}, readEndPos{readEndPosIn},
          mateChrIndex{mateChrIndexIn}, mateStartPos{mateStartPosIn},
          mateEndPos{mateStartPosIn}, inverted{invertedIn}, source{sourceType},
          evidenceLevel{sourceType == 2 ? 3 : 1}, matePower{1},
          inversionSupport{invertedIn}, straightSupport{!invertedIn},
          bpLocs{bpLocsIn}, saSupporter{false}, toRemove{false} {}

    bool isToRemove() const { return toRemove; }
};

}   // namespace sophia
#endif /* MATEINFO_H_ */
