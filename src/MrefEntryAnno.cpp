/*
 * MrefEntry.cpp
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

#include "Breakpoint.h"
#include "strtk.hpp"
#include <boost/algorithm/string/join.hpp>
#include <MrefEntryAnno.h>
#include <unordered_set>

namespace sophia {
    boost::format MrefEntryAnno::doubleFormatter { "%.5f" };
    int MrefEntryAnno::DEFAULTREADLENGTH { };
    int MrefEntryAnno::PIDSINMREF { };

    MrefEntryAnno::MrefEntryAnno(const string& mrefEntryIn) :
                    pos { 0 },
                    numHits { 0 },
                    suppAlignments { } {
        auto index = 0;
        vector<int> bpChunkPositions { };
        bpChunkPositions.reserve(7);
        auto cit = mrefEntryIn.cbegin();
        if (mrefEntryIn.back() != '.') {
            while (bpChunkPositions.size() < 8) {
                if (*cit == '\t') {
                    bpChunkPositions.push_back(index);
                }
                ++index;
                ++cit;
            }
            string saStr { };
            for (auto i = bpChunkPositions[7] + 1; i < static_cast<int>(mrefEntryIn.length()); ++i) {
                if (mrefEntryIn[i] == ';') {
                    suppAlignments.emplace_back(saStr);
                    saStr.clear();
                } else {
                    saStr.push_back(mrefEntryIn[i]);
                }
            }
            suppAlignments.emplace_back(saStr);
        } else {
            while (bpChunkPositions.size() < 5) {
                if (*cit == '\t') {
                    bpChunkPositions.push_back(index);
                }
                ++index;
                ++cit;
            }
        }

        for (auto i = bpChunkPositions[0] + 1; i < bpChunkPositions[1]; ++i) {
            pos = pos * 10 + (mrefEntryIn[i] - '0');
        }
        for (auto i = bpChunkPositions[2] + 1; i < bpChunkPositions[3]; ++i) {
            numHits = numHits * 10 + (mrefEntryIn[i] - '0');
        }
        if (mrefEntryIn[0] == 'Y') {
            numHits = min(PIDSINMREF, 2 * numHits);
        }
        for (auto &sa : suppAlignments) {
            sa.setSecondarySupport(numHits);
        }

    }

    //SuppAlignmentAnno* MrefEntryAnno::searchFuzzySa(const SuppAlignmentAnno& fuzzySa) {
    //	SuppAlignmentAnno* match = nullptr;
    //	for (auto &sa : suppAlignments) {
    //		if (sa.saCloseness(fuzzySa, 1)) {
    //			match = &sa;
    //			return match;
    //		}
    //	}
    //	return nullptr;
    //}

} /* namespace sophia */

