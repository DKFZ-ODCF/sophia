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
#include <string>
#include <vector>
#include <optional>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered_set.hpp>


namespace sophia {

    class Hg38ChrConverter: public ChrConverter {

      public:

        using ChrNames = std::vector<ChrName>;
        using ChrSizes = std::vector<ChrSize>;
        using ChrIndexMap = boost::unordered::unordered_map<ChrName, ChrIndex>;
        using ChrCompressIndexMap = boost::unordered::unordered_map<ChrName, CompressedMrefIndex>;
        using ChrSizeMap = boost::unordered::unordered_map<ChrName, ChrSize>;


      protected:

        // ChrIndex -> ChrName
        const ChrNames allChromosomes;

        // ChrName -> ChrIndex
        const ChrIndexMap allChromosomeLookup;

        // ChrIndex -> ChrSize
        const ChrSizes allChromosomeLengths;

        // These are used for efficient membership checks of chromosomes in the different sets.
        const IndexRange primaryContigRange;
        const IndexRange extrachromosomalContigRange;
        const IndexRange virusContigRange;
        const IndexRange decoyContigRange;
        const IndexRange technicalContigRange;

        /** Initialize the hg38 chromosome converter with different types of contig/chromosome
          * names and the sizes of the corresponding chromosomes.
          *
          * @param primaryContigs           Primary contigs, e.g. chr1, chr2, ..., chr22, chrX, chrY
          *                                 Must not be empty.
          * @param extrachromosomalContigs  Extrachromosomal contigs, e.g. chrM, chrMT
          *                                 May be empty.
          * @param virusContigs             Virus contigs, e.g. NC_007605, EBV. This is for viruses
          *                                 that may insert into the nuclear genome.
          *                                 May be empty.
          * @param technicalContigs         Technical contigs, e.g. phiX
          *                                 May be empty.
          * @param decoyContigs             Decoy contigs.
          *                                 May be empty.
          *
          * These sets must be mutually non-overlapping, i.e. no chromosome name may be in more than
          * one set.
          *
          * The compressedMref chromosomes will be the union of
          *  - primaryContigs
          *  - extrachromosomalContigs
          *  - virusContigs
          **/
        Hg38ChrConverter(std::string assemblyName,
                         const ChrSizeMap &primaryContigs,
                         const ChrSizeMap &extrachromosomalContigs,
                         const ChrSizeMap &virusContigs,
                         const ChrSizeMap &technicalContigs,
                         const ChrSizeMap &decoyContigs);

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

        // ... for building an instance
        static
        ChrIndexMap
        buildChrToIndex(const std::vector<ChrName> &chromosomes);

        static
        ChrSizeMap
        mergeChrSizeMaps(const ChrSizeMap &primaryContigs,
                         const ChrSizeMap &extrachromosomalContigs,
                         const ChrSizeMap &virusContigs,
                         const ChrSizeMap &decoyContigs,
                         const ChrSizeMap &technicalContigs);

        static
        ChrNames
        collect(const ChrSizeMap &map,
                std::function<ChrName(std::pair<ChrName, ChrSize>)> fun);

        ChrSizes
        buildAllChromosomeLengths(const std::vector<ChrName> &allChromosomesIn,
                                  const ChrSizeMap &primaryContigs,
                                  const ChrSizeMap &extrachromosomalContigs,
                                  const ChrSizeMap &virusContigs,
                                  const ChrSizeMap &decoyContigs,
                                  const ChrSizeMap &technicalContigs) const;

      public:

        const std::string assemblyName;

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

        /** chr1-chr22, chrX, chrY. */
        bool isPrimary(ChrIndex index) const;

        /** phix index. */
        bool isTechnical(ChrIndex index) const;

        /** NC_007605, EBV. */
        bool isVirus(ChrIndex index) const;

        /** Mitochondrial chromosome index. */
        bool isExtrachromosal(ChrIndex index) const;

        /** Decoy sequence index. */
        bool isDecoy(ChrIndex index) const;

        /** Whether the chromosome index is that of a compressed mref chromosome. */
        bool isCompressedMrefIndex(ChrIndex index) const;

        /** Map the compressed mref index to the uncompressed mref index. */
        std::optional<ChrIndex> compressedMrefIndexToIndex(CompressedMrefIndex index) const;

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