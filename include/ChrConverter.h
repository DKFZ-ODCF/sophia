/*
 *     Author: Philip R. Kensche, DKFZ Heidelberg (Omics IT and Data Management Core Facility)
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
 *     LICENSE: GPL
 */

#ifndef _CHRCONVERTER_H_
#define _CHRCONVERTER_H_

#include <stdexcept>
#include <string>
#include <vector>


namespace sophia {

    using namespace std;

    /** ChrConverter contains information the names of chromosomes in an assembly. */
    class ChrConverter {
      protected:

        /** The constructor should be used to initialize the fields from subclasses. It does
            additional checks of the dimensions of the input vectors. */
        ChrConverter(const vector<string>& indexToChr,
                     const vector<string>& indexToChrCompressedMref,
                     const vector<int>& chrSizesCompressedMref,
                     const vector<int>& indexConverter);


      public:

        virtual ~ChrConverter();

        /** The name of the assembly. */
        static const string assembly_name;

        /** Mapping indices to chromosome names. */
        const vector<string> indexToChr;

        /** Mapping indices to chromosome names for compressed mref files. */
        const vector<string> indexToChrCompressedMref;

        /** Chromosome sizes in base pairs. */
        const vector<int> chrSizesCompressedMref;

        /** Mapping chromosome names to indices. */
        const vector<int> indexConverter;

        /** Parse chromosome index.  It takes a position in a character stream, and translates the
            following character(s) into index positions (using ChrConverter::indexToChr). */
        virtual int readChromosomeIndex(string::const_iterator startIt, char stopChar) const = 0;

        size_t n_chromosomes() {
            return indexToChr.size();
        };

        size_t n_chromosomes_compressed_mref() {
            return indexToChrCompressedMref.size();
        };


    };

}

#endif /* _CHRCONVERTER_H_ */