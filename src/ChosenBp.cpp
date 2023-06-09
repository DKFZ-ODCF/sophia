/*
 * ChosenBp.cpp
 *
 *  Created on: 7 May 2016
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

#include "ChosenBp.h"
#include <algorithm>

namespace sophia {

using namespace std;

int ChosenBp::BPSUPPORTTHRESHOLD{};

void
ChosenBp::addChildNode(int indexIn) {
    childrenNodes.push_back(indexIn);
}

void
ChosenBp::addSupplementaryAlignments(
    const vector<SuppAlignment> &suppAlignments) {
    for (const auto &sa : suppAlignments) {
        auto it = find_if(supplementaryAlignments.begin(),
                          supplementaryAlignments.end(),
                          [&](const SuppAlignment &suppAlignment) {
                              return suppAlignment.saCloseness(sa, 5);
                          });
        if (it == supplementaryAlignments.end()) {
            supplementaryAlignments.push_back(sa);
        } else {
            if (it->isFuzzy() && !sa.isFuzzy()) {
                it->removeFuzziness(sa);
            } else if (it->isFuzzy() && sa.isFuzzy()) {
                it->extendSuppAlignment(sa.getPos(), sa.getExtendedPos());
            }
            // it->addSupportingIndices(sa.getSupportingIndices());
            if (sa.getMapq() > it->getMapq()) {
                it->setMapq(sa.getMapq());
            }
            it->incrementDistinctReads();
        }
    }
}

}   // namespace sophia
