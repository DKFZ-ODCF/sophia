/*
 * ChosenBp.h
 *
 *  Created on: 18 Apr 2016
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

#ifndef CHOSENBP_H_
#define CHOSENBP_H_
#include "SuppAlignment.h"
#include <string>
#include <vector>

namespace sophia {

using namespace std;

class ChosenBp {
    friend class Alignment;

  public:
    ChosenBp(char bpTypeIn, int bpSizeIn, bool bpEncounteredMIn,
             int overhangStartIndexIn, int overhangLengthIn,
             int selfNodeIndexIn)
        : bpType{bpTypeIn}, bpSize{bpSizeIn}, bpEncounteredM{bpEncounteredMIn},
          overhangStartIndex{overhangStartIndexIn},
          overhangLength{overhangLengthIn}, supplementaryAlignments{},
          childrenNodes{{selfNodeIndexIn}}, selfNodeIndex{selfNodeIndexIn} {}
    ~ChosenBp() = default;
    static int BPSUPPORTTHRESHOLD;

  private:
    char bpType;
    int bpSize;
    bool bpEncounteredM;
    int overhangStartIndex, overhangLength;
    vector<SuppAlignment> supplementaryAlignments;
    vector<int> childrenNodes;
    int selfNodeIndex;
    void addChildNode(int indexIn);
    void
    addSupplementaryAlignments(const vector<SuppAlignment> &suppAlignments);
};

}   // namespace sophia

#endif /* CHOSENBP_H_ */
