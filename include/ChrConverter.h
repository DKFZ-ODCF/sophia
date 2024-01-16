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

#ifndef CHRCONVERTER_H_
#define CHRCONVERTER_H_

#include "global.h"
#include <stdexcept>
#include <string>
#include <vector>
#include <optional>


namespace sophia {

    /** ChrConverter manages information on chromosomes names, sizes, and index positions in
      * data arrays.
      *
      * This probably needs a redesign. The current situation is just an intermediate step away
      * from the former completely procedural implementation that heavily leaked implementation
      * details into the calling code and was highly tuned but very inflexible. */
    class ChrConverter {

      public:

        virtual ~ChrConverter();

        /** The name of the assembly. */
        static const std::string assemblyName;

        /** Number of chromosomes. */
        virtual ChrIndex nChromosomes() const = 0;

        /** Number of compressed mref chromosomes. */
        virtual CompressedMrefIndex nChromosomesCompressedMref() const = 0;

        /** Map an index position to a chromosome name. */
        virtual std::string indexToChrName(ChrIndex index) const = 0;

        /** Map an index position to a chromosome name for compressed mref files. */
        virtual std::string indexToChrNameCompressedMref(CompressedMrefIndex index) const = 0;
//
//        /** phix index. */
//        virtual bool isTechnical(ChrIndex index) const = 0;
//
//        /** NC_007605, EBV. */
//        virtual bool isVirus(ChrIndex index) const = 0;
//
//        /** Mitochondrial chromosome index. */
//        virtual bool isExtrachromosal(ChrIndex index) const = 0;
//
//          /** Decoy sequence index. */
//          virtual bool isDecoy(ChrIndex index) const = 0;


        /** Whether the chromosome index is that of a compressed mref chromosome. */
        virtual bool isCompressedMrefIndex(ChrIndex index) const = 0;

        /** Whether the chromosome index is that of an ignored chromosome. Ignored chromosomes
          * are not the same as the ones that are not among the compressedMref chromosomes.
          * This should be the index of the phiX chromosomes. */
        virtual bool isIgnoredChromosome(ChrIndex index) const = 0;

        /** Map the compressed mref index to the uncompressed mref index. */
        virtual std::optional<ChrIndex>
        compressedMrefIndexToIndex(CompressedMrefIndex index) const = 0;

        /** Map compressed mref index to chromosome size. */
        virtual ChrSize chrSizeCompressedMref(CompressedMrefIndex index) const = 0;

        /** Map a chromosome name to an index position. */
        virtual ChrIndex chrNameToIndex(std::string chrName) const = 0;

        /** Parse chromosome index.
          *
          * 1. Input is a plain chromosome string separated from the following string by '\t' when
          *    parsing a BED file (the chromosome identifier in the first column).
          * 2. ...
          *
          * If the `stopCharExt` parameter is an empty string, then it takes a position in a
          * character stream, and translates the following character(s) into index positions
          * (using ChrConverter::indexToChr). If the name cannot be parsed, throws a domain_error
          * exception.
          *
          * IMPORTANT: Implementations may or may not use the `stopCharExt` parameter. Therefore,
          *            the following behavior is optional. An implementation may not even actually
          *            validate that the stopChar or string-end terminates the parsed identifier!
          *
          * If the `stopCharExt` parameter is *not* empty, the method first parses up to the first
          * occurrence of the `stopCharExt`. Then within the identified start and end range, parses
          * up to the last occurrence of `stopChar`. This allows to parse a chromosome name
          * "HLA-DRB1*13:01:01" from a string "HLA-DRB1*13:01:01:2914|(4,0,0?/0)" by first
          * separating out the `|` separator (stopCharExt), and then finding the last `:`
          * separator (stopChar).
          *
          **/
        virtual ChrIndex
        parseChrAndReturnIndex(std::string::const_iterator startIt,
                               std::string::const_iterator endIt,
                               char stopChar,
                               const std::string &stopCharExt = "") const = 0;

    };

}

#endif /* CHRCONVERTER_H_ */