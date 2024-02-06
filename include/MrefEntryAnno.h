/*
 * MrefEntry.h
 *
 *  Created on: 27 Nov 2016
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

#ifndef MREFENTRYANNO_H_
#define MREFENTRYANNO_H_

#include "SuppAlignmentAnno.h"
#include <BreakpointReduced.h>
#include <boost/format.hpp>
#include <string>

namespace sophia {

using namespace std;

class MrefEntryAnno {

  public:

    static int PIDS_IN_MREF;

    static ChrSize DEFAULT_READ_LENGTH;

    static boost::format doubleFormatter;

    MrefEntryAnno(const string &mrefEntryIn);

    template <typename T> bool operator<(const T &rhs) const {
        return static_cast<int>(pos) < static_cast<int>(rhs.getPos());
    }

    template <typename T> int distanceTo(const T &rhs) const {
        return abs(static_cast<int>(pos) - static_cast<int>(rhs.getPos()));
    }
    template <typename T> int distanceToBp(const T &compIn) const {
        return abs(static_cast<int>(pos) - static_cast<int>(compIn.getPos()));
    }

    bool operator==(const MrefEntryAnno &rhs) const {
        return pos == rhs.getPos();
    }

    ChrSize getPos() const { return pos; }

    vector<SuppAlignmentAnno *> getSuppAlignmentsPtr() {
        vector<SuppAlignmentAnno *> res{};
        for (auto &sa : suppAlignments) {
            res.push_back(&sa);
        }
        return res;
    }

    void removeMarkedFuzzies() {
        while (!suppAlignments.empty() && suppAlignments.back().isToRemove()) {
            suppAlignments.pop_back();
        }
        for (auto saIt = suppAlignments.begin(); saIt != suppAlignments.end();
             ++saIt) {
            if (saIt->isToRemove()) {
                swap(*saIt, suppAlignments.back());
            }
            while (!suppAlignments.empty() &&
                   suppAlignments.back().isToRemove()) {
                suppAlignments.pop_back();
            }
        }
    }
    //	SuppAlignmentAnno* searchFuzzySa(const SuppAlignmentAnno& fuzzySa);

    const vector<SuppAlignmentAnno> &getSuppAlignments() const {
        return suppAlignments;
    }

    vector<SuppAlignmentAnno *> getSupplementsPtr() {
        vector<SuppAlignmentAnno *> res{};
        for (auto &sa : suppAlignments) {
            res.push_back(&sa);
        }
        return res;
    }

    bool closeToSupp(const SuppAlignmentAnno &compIn, ChrDistance fuzziness) const {
        if (compIn.isFuzzy()) {
            fuzziness = int(2.5 * DEFAULT_READ_LENGTH);   /* truncate */
            return (static_cast<long>(pos) - fuzziness) <= (static_cast<long>(compIn.getExtendedPos()) + fuzziness) &&
                   (static_cast<long>(compIn.getPos()) - fuzziness) <= (static_cast<long>(pos) + fuzziness);
        } else {
            return ChrDistance(abs(static_cast<long>(pos) - static_cast<long>(compIn.getPos()))) <= fuzziness;
        }
    }

    ChrDistance distanceToSupp(const SuppAlignmentAnno &compIn) const {
        ChrDistance result;
        if (compIn.isFuzzy()) {
            if (compIn.getPos() <= pos && pos <= compIn.getExtendedPos()) {
                result = 0;
            } else {
                if (pos < compIn.getPos()) {
                    result = ChrDistance(compIn.getPos() - pos);
                } else {
                    result = ChrDistance(pos - compIn.getExtendedPos());
                }
            }
        } else {
            result = ChrDistance(abs(static_cast<long>(pos) - static_cast<long>(compIn.getPos())));
        }
        return result;
    }

    short getNumHits() const { return numHits; }

    void setNumHits(short numHits) { this->numHits = numHits; }

  private:
    ChrSize pos;
    short numHits;
    vector<SuppAlignmentAnno> suppAlignments;
};

} /* namespace sophia */

#endif /* MREFENTRYANNO_H_ */
