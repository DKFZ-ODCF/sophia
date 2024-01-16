/*
 * SamSegmentMapper.h
 *
 *  Created on: 16 Apr 2016
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

#ifndef SAMSEGMENTMAPPER_H_
#define SAMSEGMENTMAPPER_H_
#include "global.h"
#include "Breakpoint.h"
#include "CoverageAtBase.h"
#include "MateInfo.h"
#include <ctime>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace sophia {

using namespace std;

    class SamSegmentMapper {
      public:

        SamSegmentMapper(int defaultReadLengthIn);

        ~SamSegmentMapper() = default;

        void parseSamStream();

      private:

        // Does not print anything by itself, but lets via via another call to
        // Breakpoint::finalizeBreakpoint and then Breakpoint::printBreakpointReport it
        // prints to stdout.
        void printBps(int alignmentStart);

        void switchChromosome(const Alignment &alignment);

        void incrementCoverages(const Alignment &alignment);

        void assignBps(shared_ptr<Alignment> &alignment);

        const time_t STARTTIME;

        const bool PROPERPARIRCOMPENSATIONMODE;

        const int DISCORDANTLEFTRANGE;

        const int DISCORDANTRIGHTRANGE;

        unsigned int printedBps;

        ChrIndex chrIndexCurrent;
        int minPos, maxPos;
        map<int, Breakpoint> breakpointsCurrent;
        deque<CoverageAtBase> coverageProfiles;
        deque<MateInfo> discordantAlignmentsPool;
        deque<MateInfo> discordantAlignmentCandidatesPool;
        deque<MateInfo> discordantLowQualAlignmentsPool;
    };

} /* namespace sophia */

#endif /* SAMSEGMENTMAPPER_H_ */
