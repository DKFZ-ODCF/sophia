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

#ifndef MREFENTRY_H_
#define MREFENTRY_H_

#include "BreakpointReduced.h"
#include "SuppAlignment.h"
#include <boost/format.hpp>
#include <string>

namespace sophia {

    using namespace std;

    /**
     * @brief MrefEntry class is a container for the mref entries
     * It contains the position, the file indices, the artifact ratios and the support alignments.
     * This class should be memory optimized, as it will be instantiated in the billions, once
     * for every chromosome position by the MasterRefProcessor.
     */
    class MrefEntry {
      public:

        using ValidityScore = signed char;

        static unsigned int NUM_PIDS;

        static ChrSize DEFAULT_READ_LENGTH;

        static boost::format doubleFormatter;

        MrefEntry();

        void addEntry(Breakpoint &tmpBreakpoint, int fileIndex);

        void addEntry(BreakpointReduced &tmpBreakpoint, int fileIndex);

        void mergeMrefEntries(MrefEntry &entry2);

        ChrSize getPos() const {
            if (!isValid()) {
                throw_with_trace(logic_error("MrefEntry is invalid"));
            }
            return pos;
        }

        const vector<float> &getArtifactRatios() const { return artifactRatios; }

        const vector<unsigned short> &getFileIndices() const { return fileIndices; }

        ValidityScore getValidityScore() const { return validity; }

        void removeMarkedFuzzies() {
            suppAlignments.erase(remove_if(suppAlignments.begin(),
                                           suppAlignments.end(),
                                           [](const SuppAlignmentAnno &sa) {
                                               return sa.isToRemove();
                                           }),
                                 suppAlignments.end());
        }

        string printBpInfo(const string &chromosome);

        string printArtifactRatios(const string &chromosome);

        SuppAlignmentAnno *searchFuzzySa(const SuppAlignmentAnno &fuzzySa);

        vector<SuppAlignmentAnno *> getSupplementsPtr() {
            vector<SuppAlignmentAnno *> res{};
            for (auto &sa : suppAlignments) {
                res.push_back(&sa);
            }
            return res;
        }

        const vector<unsigned short> &getFileIndicesWithArtifactRatios() const {
            return fileIndicesWithArtifactRatios;
        }

        const vector<SuppAlignmentAnno> &getSuppAlignments() const {
            return suppAlignments;
        }

        void setAsInvalid() {
            pos = std::numeric_limits<ChrSize>::max();
            validity = -1;
        }

        bool isValid() const {
            return pos != std::numeric_limits<ChrSize>::max();
        }

      private:

        bool saMatcher(SuppAlignmentAnno *saPtr);

        void finalizeFileIndices();

        ValidityScore validity;   // -1 nothing, 0 only sa, 1 sa and support

        ChrSize pos;

        vector<unsigned short> fileIndices;

        vector<unsigned short> fileIndicesWithArtifactRatios;

        vector<float> artifactRatios;

        vector<SuppAlignmentAnno> suppAlignments;
    };

} /* namespace sophia */

#endif /* MREFENTRY_H_ */
