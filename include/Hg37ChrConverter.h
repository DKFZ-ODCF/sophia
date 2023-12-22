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

#ifndef HG37CHRCONVERTER_H_
#define HG37CHRCONVERTER_H_

#include "ChrConverter.h"
#include "global.h"
#include <vector>
#include <string>
#include <optional>


namespace sophia {

    /** Hard-coded chromosome converter for hg37. This tries to encapsulate the implementation
        details of the original version. */
    class Hg37ChrConverter: public ChrConverter {
      protected:

        /** The constructor does additional checks of the dimensions of the input vectors. */
        Hg37ChrConverter(const std::vector<std::string>& indexToChr,
                         const std::vector<std::string>& indexToChrCompressedMref,
                         const std::vector<ChrSize>& chrSizesCompressedMref,
                         const std::vector<ChrIndex>& indexConverter);

        /** Mapping indices to chromosome names. */
        const std::vector<std::string> indexToChr;

        /** Mapping indices to chromosome names for compressed mref indices. */
        const std::vector<std::string> indexToChrCompressedMref;

        /** Chromosome sizes in base pairs, only for compressed mref chromosomes. */
        const std::vector<CompressedMrefIndex> chrSizesCompressedMref;

        /** Mapping compressed mref indices names to indices. */
        const std::vector<ChrIndex> indexConverter;

        ChrIndex parseChrAndReturnIndexImpl(std::string::const_iterator start,
                                            std::string::const_iterator end,
                                            char stopChar) const;

        static const ChrIndex phixChrIndex;

      public:

        static const std::string assemblyName;

        Hg37ChrConverter();

        /** Return the number of chromosomes. */
        ChrIndex nChromosomes() const;

        /** Number of compressed mref chromosomes. */
        CompressedMrefIndex nChromosomesCompressedMref() const;

        /** Map an index position to a chromosome name. */
        std::string indexToChrName(ChrIndex index) const;

        /** Map an index position to a chromosome name for compressed mref files. */
        std::string indexToChrNameCompressedMref(CompressedMrefIndex index) const;

        /** Whether the chromosome index is that of an ignored chromosome. Ignored chromosomes
          * are not the same as the ones that are not among the compressedMref chromosomes.
          * For hg37 this used to be phiX only. */
        bool isIgnoredChromosome(ChrIndex index) const;

        /** Map the compressed mref index to the uncompressed mref index. */
        std::optional<ChrIndex> compressedMrefIndexToIndex(CompressedMrefIndex index) const;

        /** Map compressed mref index to chromosome size. */
        ChrSize chrSizeCompressedMref(CompressedMrefIndex index) const;

        /** Map a chromosome name to an index position for compressed mref files. */
        CompressedMrefIndex chrNameToIndexCompressedMref(std::string chrName) const;

        /** Map a chromosome name to an index position. */
        ChrIndex chrNameToIndex(std::string chrName) const;

        /** Parse chromosome name given a iterator (start) and termination character.
            Validate against the pre-declared chromosome names. */
        ChrIndex parseChrAndReturnIndex(std::string::const_iterator startIt, char stopChar) const;

    };

}

#endif /* HG37CHRCONVERTER_H_ */