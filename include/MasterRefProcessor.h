/*
 * MasterRefProcessor.h
 *
 *  Created on: 27 Apr 2016
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

#ifndef MASTERREFPROCESSOR_H_
#define MASTERREFPROCESSOR_H_

#include "SuppAlignment.h"
#include "global.h"
#include <BreakpointReduced.h>
#include <MrefEntry.h>
#include <array>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace sophia {

    class MasterRefProcessor {
      public:
        MasterRefProcessor(const std::vector<std::string> &filesIn,
                           const std::string &outputRootName,
                           const std::string &version,
                           const ChrSize defaultReadLengthIn);

        ~MasterRefProcessor() = default;

      private:
        // Note that all methods and fields are private.
        // The MasterRefProcessor does all the work during construction time.

        unsigned long long processFile(const std::string &gzPath, short fileIndex);
        bool processBp(BreakpointReduced &bp, ChrIndex chrIndex, short fileIndex);

        const int NUM_PIDS;
        const ChrSize DEFAULT_READ_LENGTH;
        std::unique_ptr<std::ofstream> mergedBpsOutput;

        /** This will be a huge data structure, that contains one MrefEntry per position in the
         *  master reference chromosomes.
         **/
        std::vector<std::vector<MrefEntry>> mrefDb;
    };

}   // namespace sophia

#endif /* MASTERREFPROCESSOR_H_ */
