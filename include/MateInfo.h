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
#include "global.h"
#include <cmath>

namespace sophia {

    struct MateInfo {

        ChrSize readStartPos;
        ChrSize readEndPos;
        ChrIndex mateChrIndex;
        ChrSize mateStartPos;
        ChrSize mateEndPos;
        bool inverted;
        int source;
        int evidenceLevel;
        int matePower;
        int inversionSupport;
        int straightSupport;
        std::vector<ChrSize> bpLocs;
        bool saSupporter;
        bool toRemove;

        bool operator<(const MateInfo &rhs) const;

        bool suppAlignmentFuzzyMatch(const SuppAlignment &sa) const;

        MateInfo(ChrSize readStartPosIn,
                 ChrSize readEndPosIn,
                 ChrIndex mateChrIndexIn,
                 ChrSize mateStartPosIn,
                 int sourceType,
                 bool invertedIn);

        MateInfo(ChrSize readStartPosIn,
                 ChrSize readEndPosIn,
                 ChrIndex mateChrIndexIn,
                 ChrSize mateStartPosIn,
                 int sourceType,
                 bool invertedIn,
                 const std::vector<ChrSize> &bpLocsIn);

        bool isToRemove() const;
        
    };

}   // namespace sophia
#endif /* MATEINFO_H_ */
