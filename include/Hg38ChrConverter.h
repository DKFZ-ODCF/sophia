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
#include "global.h"
#include <string>
#include <vector>
#include <optional>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered_set.hpp>


namespace sophia {

    class Hg38ChrConverter: public ChrConverter {

      protected:

        /* This makes a search structure to map chromosome names to index positions in an array.
           This is useful mostly during the initial parsing of chromosome names from the input.
           Afterwards, we continue working with the positions. */
        boost::unordered::unordered_map<std::string, std::vector<std::string>::size_type>
        makeChrToIndex(const std::vector<std::string>& chromosomes) const;

        /** Given the index-maps for all chromosomes, and the selected "compressed mref" chromosomes,
            create a mapping of mref indices to all chromosome indices. */
        std::vector<ChrIndex> makeCompressedMrefIndexToIndexConverter(
            const boost::unordered::unordered_map<std::string, ChrIndex>
            &allChromosomeLookup,
            const boost::unordered::unordered_map<std::string, CompressedMrefIndex>
            &compressedMrefChromosomeLookup
            ) const;

        // All chromosomes mappings.
        const std::vector<std::string> allChromosomes;
        const boost::unordered::unordered_map<std::string, ChrIndex>
            allChromosomeLookup;
        const std::vector<ChrSize> allChromosomeLengths;

        // Mapping compressed Mref chromosome indices to all chromosome indices.
        const boost::unordered::unordered_map<std::string, CompressedMrefIndex>
            compressedMrefChromosomeLookup;

        // Compressed Mref chromosome mappings.
        const std::vector<ChrIndex> compressedMrefIndexToIndexLookup;
        const std::vector<std::string> compressedMrefChromosomes;

        // Ignored chromosomes.
        const boost::unordered_set<std::string> ignoredChromosomes;

        Hg38ChrConverter(const std::vector<std::string> &allChromosomes,
                         const std::vector<ChrSize> &allChromosomeLengths,
                         const std::vector<std::string> &compressedMrefChromosomes,
                         const std::vector<std::string> &ignoredChromosomes);

      public:

        static const std::string assemblyName;

        /** This default constructor only makes sense, as long as the hg38 chromosome names are
            hard-coded. */
        Hg38ChrConverter();

        /** Number of chromosomes. */
        ChrIndex nChromosomes() const;

        /** Number of compressed mref chromosomes. */
        CompressedMrefIndex nChromosomesCompressedMref() const;

        /** Map an index position to a chromosome name. */
        std::string indexToChrName(ChrIndex index) const;

        /** Map an index position to a chromosome name for compressed mref files. */
        std::string indexToChrNameCompressedMref(CompressedMrefIndex index) const;

        /** Whether the chromosome index is that of an ignored chromosome. Ignored chromosomes
          * are not the same as the ones that are not among the compressedMref chromosomes.
          * This should include phiX. */
        bool isIgnoredChromosome(ChrIndex index) const;

        /** Map the compressed mref index to the uncompressed mref index. */
        std::optional<ChrIndex> compressedMrefIndexToIndex(CompressedMrefIndex index) const;

        /** Map compressed mref index to chromosome size. */
        ChrSize chrSizeCompressedMref(CompressedMrefIndex index) const;

        /** Map a chromosome name to an index position. */
        ChrIndex chrNameToIndex(std::string chrName) const;

        /** Parse chromosome index. It takes a position in a character stream, and translates the
          * following character(s) into index positions (using ChrConverter::indexToChr).
          * If the name cannot be parsed, throws a domain_error exception.
          *
          * This method parses up to the first occurrence of the `stopCharExt`. Then within the
          * identified start and end range, parses up to the last occurrence of `stopChar`. This
          * allows to parse a chromosome name "HLA-DRB1*13:01:01" from a string
          * "HLA-DRB1*13:01:01:2914|(4,0,0?/0)" by first separating out the `|` separator
          * (stopCharExt), and then finding the last `:` separator (stopChar).
          **/
        ChrIndex parseChrAndReturnIndex(std::string::const_iterator startIt,
                                        std::string::const_iterator endIt,
                                        char stopChar,
                                        const std::string &stopCharExt = "\0") const;

    };

}

#endif /* HG38CHRCONVERTER_H_ */