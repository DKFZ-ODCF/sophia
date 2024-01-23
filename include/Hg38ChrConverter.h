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
#include "IndexRange.h"
#include "ChrCategory.h"
#include "ChrInfo.h"
#include "ChrInfoTable.h"
#include <string>
#include <vector>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered_set.hpp>


namespace sophia {

    /** This converter provides fast by-index access to the ChrInfoTable. */
    class Hg38ChrConverter: public ChrConverter {

      public:

        using ChrToIndexMap =
            boost::unordered::unordered_map<ChrName, ChrIndex>;
        using CompressedMrefChrToIndexMap =
            boost::unordered::unordered_map<ChrName, CompressedMrefIndex>;

      protected:

        /** The ChrInfoTable provides some access to chromosomes and chromosome categories and
          * guarantees a consistent order of chromosomes. Thus, we use the ChrInfoTable as the
          * to access chromosome names and sizes using their global index. */
        const ChrInfoTable chrInfoTable;

        /** ChrInfoTable, is *not* actually for the index-based access. Therefore, as we need
          * a ChrName -> ChrIndex mapping, we manage this mapping here. */
        const ChrToIndexMap allChromosomeLookup;
        static ChrToIndexMap buildAllChromosomeLookup(const ChrInfoTable::ChrNames &chr_info);

        /** A mapping table to convert the compressed mref indices into the global index space */
        const std::vector<ChrIndex> compressedToAllMapping;
        static std::vector<ChrIndex> buildCompressedMrefToAllMapping(ChrInfoTable chrInfoIn);

        // Helper functions

        // ... for parsing
        static
        ChrName
        parseChrBreakPoint(std::string::const_iterator startIt,
                           std::string::const_iterator endIt,
                           char stopChar,
                           const std::string &stopCharExt);

        static
        ChrName
        parseChrSimple(std::string::const_iterator startIt,
                       std::string::const_iterator endIt,
                       char stopChar);

      public:

        const std::string assemblyName;

        /** Initialize the hg38 chromosome converter with different types of contig/chromosome
          * names and the sizes of the corresponding chromosomes.
          **/
        Hg38ChrConverter(std::string assemblyName,
                         ChrInfoTable chrInfo);

        /** This default constructor only makes sense, as long as the hg38 chromosome names are
            hard-coded. */
        Hg38ChrConverter();

        /** Number of chromosomes. */
        ChrIndex nChromosomes() const;

        /** Number of compressed mref chromosomes. */
        CompressedMrefIndex nChromosomesCompressedMref() const;

        /** Map an index position to a chromosome name. */
        ChrName indexToChrName(ChrIndex index) const;

        /** Map an index position to a chromosome name for compressed mref files. */
        ChrName indexToChrNameCompressedMref(CompressedMrefIndex index) const;

        // The following methods could also be implemented as isCategory(ChrIndex, ChrCategory),
        // but, for performance reason we provide them as separate methods.

        /** chr1-chr22 */
        bool isAutosome(ChrIndex index) const;

        /** chrX, chrY */
        bool isGonosome(ChrIndex index) const;

        /** phix index. */
        bool isTechnical(ChrIndex index) const;

        /** NC_007605, EBV. */
        bool isVirus(ChrIndex index) const;

        /** Mitochondrial chromosome index. */
        bool isExtrachromosomal(ChrIndex index) const;

        /** Decoy sequence index. */
        bool isDecoy(ChrIndex index) const;

        /** HLA chromosome index. */
        bool isHLA(ChrIndex index) const;

        /** ALT chromosome index. */
        bool isALT(ChrIndex index) const;

        /** Unplaced chromosome index. */
        bool isUnassigned(ChrIndex index) const;

        /** Whether the chromosome index is that of a compressed mref chromosome. */
        bool isCompressedMrefIndex(ChrIndex index) const;

        /** Map the compressed mref index to the uncompressed mref index. */
        ChrIndex compressedMrefIndexToIndex(CompressedMrefIndex index) const;

        /** Map compressed mref index to chromosome size. */
        ChrSize chrSizeCompressedMref(CompressedMrefIndex index) const;

        /** Map a chromosome name to an index position. */
        ChrIndex chrNameToIndex(ChrName chrName) const;

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
        static
        ChrName parseChr(std::string::const_iterator startIt,
                         std::string::const_iterator endIt,
                         char stopChar,
                         const std::string &stopCharExt = "");

        // The same as `parseChr`, but returns the index instead of the name.
        ChrIndex parseChrAndReturnIndex(std::string::const_iterator startIt,
                                        std::string::const_iterator endIt,
                                        char stopChar,
                                        const std::string &stopCharExt = nullptr) const;

    };

}

#endif /* HG38CHRCONVERTER_H_ */