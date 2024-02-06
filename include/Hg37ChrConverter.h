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


namespace sophia {

    /** Hard-coded chromosome converter for hg37. This tries to encapsulate the implementation
        details of the original version. */
    class Hg37ChrConverter: public ChrConverter {
      protected:

        static std::vector<ChrIndex> buildCompressedMrefIndexToIndex(
            const std::vector<CompressedMrefIndex> &indexToCompressedMrefIndex);

        /** The constructor does additional checks of the dimensions of the input vectors. */
        Hg37ChrConverter(const std::vector<std::string>& indexToChrName,
                         const std::vector<std::string>& indexToChrCompressedMref,
                         const std::vector<ChrSize>& chrSizesCompressedMref,
                         const std::vector<CompressedMrefIndex>& indexToCompressedMrefIndex);

        /** Mapping indices to chromosome names. */
        const std::vector<std::string> _indexToChrName;

        /** Mapping indices to chromosome names for compressed mref indices. */
        const std::vector<std::string> _compressedMrefIndexToChrName;

        /** Chromosome sizes in base pairs, only for compressed mref chromosomes. */
        const std::vector<ChrSize> _chrSizesCompressedMref;

        /** Mapping compressed mref indices names to indices. */
        const std::vector<CompressedMrefIndex> _indexToCompressedMrefIndex;

        /* Mapping of compressed mref indices to indices. */
        const std::vector<ChrIndex> _compressedMrefIndexToIndex;

        inline static bool isValid(ChrIndex index);

        inline static void assertValid(ChrIndex index);

        inline static bool isValid(CompressedMrefIndex index);

        inline static void assertValid(CompressedMrefIndex index);

        // The following static methods are used for checks during construction, but also
        // to implement the public interface.

        /** chr1-chr22, GL00+ */
        inline static bool _isAutosome(ChrIndex index);

        /** chrX, Y, ... */
        inline static bool _isGonosome(ChrIndex index);

        /** phix index. */
        inline static bool _isTechnical(ChrIndex index);

        /** NC_007605. */
        inline static bool _isVirus(ChrIndex index);

        /** Mitochondrial chromosome index. */
        inline static bool _isExtrachromosomal(ChrIndex index);

        /** Decoy sequence index. */
        inline static bool _isDecoy(ChrIndex index);

        /** GL00.+ */
        inline static bool _isUnassigned(ChrIndex index);

        /** none */
        inline static bool _isALT(ChrIndex index);

        /** none */
        inline static bool _isHLA(ChrIndex index);

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
        std::string compressedMrefIndexToChrName(CompressedMrefIndex index) const;

        /** chr1-chr22, GL00+ */
        bool isAutosome(ChrIndex index) const;

        /** chrX, Y, ... */
        bool isGonosome(ChrIndex index) const;

        /** phix index. */
        bool isTechnical(ChrIndex index) const;

        /** NC_007605. */
        bool isVirus(ChrIndex index) const;

        /** Mitochondrial chromosome index. */
        bool isExtrachromosomal(ChrIndex index) const;

        /** Decoy sequence index. */
        bool isDecoy(ChrIndex index) const;

        /** GL00.+ */
        bool isUnassigned(ChrIndex index) const;

        /** none */
        bool isALT(ChrIndex index) const;

        /** none */
        bool isHLA(ChrIndex index) const;

        /** Whether the chromosome index is that of a compressed mref chromosome. */
        bool isCompressedMref(ChrIndex index) const;

        /** Map the compressed mref index to the uncompressed mref index. */
        ChrIndex compressedMrefIndexToIndex(CompressedMrefIndex index) const;

        /** Map an index from the global index-space to the compressed mref index-space. */
        CompressedMrefIndex indexToCompressedMrefIndex(ChrIndex index) const;

        /** Map compressed mref index to chromosome size. */
        ChrSize chrSizeCompressedMref(CompressedMrefIndex index) const;

        /** Map a chromosome name to an index position for compressed mref files. */
        CompressedMrefIndex chrNameToIndexCompressedMref(std::string chrName) const;

        /** Map a chromosome name to an index position. */
        ChrIndex chrNameToIndex(std::string chrName) const;

        bool isInBlockedRegion(ChrIndex chrIndex, ChrSize position) const;

        /* This is parsing code. It takes a position in a character stream, and translates the
           following character(s) into index positions (see ChrConverter::indexToChrName). It is
           slightly modified from the original implementation by Umut Toprak.

           If the first position is a digit, read up to the next stopChar.

             * (\d+)$ -> $1

           If the first position is *not* a digit return indices according to the following rules:

             * h -> 999
             * X -> 40
             * Y -> 41
             * MT -> 1001
             * G?(\d+)\. -> $1
             * N -> 1000
             * p -> 1002

           NOTE: Most of the matches are eager matches, which means the algorithm does not check for
                 whether the end iterator or the stopChar is actually reached! The actual stopChar is
                 not actually checked in these cases.

           All identifiers not matching any of these rules, with throw an exception (domain_error).

           IMPORTANT: The hg37 parser, here, ignores the `stopCharExt`, but instead keeps the
                      legacy behavior only using the `stopChar`
        */
        ChrIndex parseChrAndReturnIndex(std::string::const_iterator startIt,
                                        std::string::const_iterator endIt,
                                        char stopChar,
                                        const std::string &stopCharFirst = "") const;

    };

}

#endif /* HG37CHRCONVERTER_H_ */