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
        return (int) pos < (int) rhs.getPos();
    }

    template <typename T> int distanceTo(const T &rhs) const {
        return abs((int) pos - (int) rhs.getPos());
    }
    template <typename T> int distanceToBp(const T &compIn) const {
        return abs((int) pos - (int) compIn.getPos());
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

    bool closeToSupp(const SuppAlignmentAnno &compIn, int fuzziness) const {
        if (compIn.isFuzzy()) {
            fuzziness = 2.5 * DEFAULT_READ_LENGTH;
            return ((long) pos - (long) fuzziness) <= ((long) compIn.getExtendedPos() + (long) fuzziness) &&
                   ((long) compIn.getPos() - (long) fuzziness) <= ((long) pos + (long) fuzziness);
        } else {
            return abs((long) pos - (long) compIn.getPos()) <= (long) fuzziness;
        }
    }

    ChrSize distanceToSupp(const SuppAlignmentAnno &compIn) const {
        ChrSize result;
        if (compIn.isFuzzy()) {
            if (compIn.getPos() <= pos && pos <= compIn.getExtendedPos()) {
                result = 0;
            } else {
                if (pos < compIn.getPos()) {
                    result = ChrSize(compIn.getPos() - pos);
                } else {
                    // TODO Why here getExtendedPos, and getPos() above?
                    result = ChrSize(pos - compIn.getExtendedPos());
                }
            }
        } else {
            result = ChrSize(abs((long) pos - (long) compIn.getPos()));
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
