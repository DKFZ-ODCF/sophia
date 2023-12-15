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

#ifndef HG38CHRCONVERTER_H_
#define HG38CHRCONVERTER_H_

#include "ChrConverter.h"
#include <string>
#include <vector>
#include <boost/unordered/unordered_map.hpp>


namespace sophia {

    using namespace boost::unordered;

    class Hg38ChrConverter: public ChrConverter {

      protected:

        /* This makes a search structure to map chromosome names to index positions in an array.
           This is useful mostly during the initial parsing of chromosome names from the input.
           Afterwards, we continue working with the positions. */
        unordered_map<std::string, std::vector<std::string>::size_type> makeChrToIndex(
            const std::vector<std::string>& chromosomes
            ) const;

        /** Given the index-maps for all chromosomes, and the selected "compressed mref" chromosomes,
            create a mapping of mref indices to all chromosome indices. */
        std::vector<ChrIndex> makeCompressedMrefIndexToIndexConverter(
            const unordered_map<std::string, ChrIndex> &allChromosomeLookup,
            const unordered_map<std::string, CompressedMrefIndex> &compressedMrefChromosomeLookup
            ) const;

        const unordered_map<std::string, ChrIndex> allChromosomeLookup;

        const unordered_map<std::string, CompressedMrefIndex> compressedMrefChromosomeLookup;

        const std::vector<int> allChromosomeLengths;

        const std::vector<ChrIndex> compressedMrefIndexToIndexLookup;
        const std::vector<std::string> compressedMrefChromosomes;

        Hg38ChrConverter(const std::vector<std::string> &allChromosomes,
                         const std::vector<int> &allChromosomeLengths,
                         const std::vector<std::string> &compressedMrefChromosomes);

      public:

        static const std::string assemblyName;

        /** This default constructor only makes sense, as long as the hg38 chromosome names are
            hard-coded. */
        Hg38ChrConverter();

        /** Number of chromosomes. */
        int nChromosomes() const;

        /** Number of compressed mref chromosomes. */
        int nChromosomesCompressedMref() const;

        /** Map an index position to a chromosome name. */
        std::string indexToChrName(ChrIndex index) const;

        /** Map an index position to a chromosome name for compressed mref files. */
        std::string indexToChrNameCompressedMref(CompressedMrefIndex index) const;

        /** Map the compressed mref index to the uncompressed mref index. */
        ChrIndex compressedMrefIndexToIndex(CompressedMrefIndex index) const;

        /** Map compressed mref index to chromosome size. */
        int chrSizeCompressedMref(CompressedMrefIndex index) const;

        /** Map a chromosome name to an index position. */
        ChrIndex chrNameToIndex(std::string chrName) const;

        /** Parse chromosome index. It takes a position in a character stream, and translates the
            following character(s) into index positions (using ChrConverter::indexToChr).
            If the name cannot be parsed, throws a domain_error exception. */
        ChrIndex parseChrAndReturnIndex(std::string::const_iterator startIt,
                                        char stopChar) const;

    };

}

#endif /* HG38CHRCONVERTER_H_ */