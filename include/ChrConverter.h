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


namespace sophia {

    /** ChrConverter manages information on chromosomes names, sizes, and index positions in
      * data arrays.
      *
      * The concept of compressed master ref (mref) chromosomes is used to separate out a
      * subset of particularly important chromosomes for which SOPHIA will calculate its
      * statistics. Note that sophiaMref will allocate a huge vector of the size of the sum
      * of the lengths of the compressed mref chromosomes (see `MasterRefProcessor` constructor.
      *
      * This class also is managing the mapping between all chromosome indices and indices of
      * the compressed master ref chromosomes (called "all index space" and "compressed master-ref
      * index space").
      **/
    class ChrConverter {

      public:

        virtual ~ChrConverter();

        /** The name of the assembly. */
        static const std::string assemblyName;

        /** Number of chromosomes. */
        virtual ChrIndex nChromosomes() const = 0;

        /** Map an index position to a chromosome name. */
        virtual ChrName indexToChrName(ChrIndex index) const = 0;

        /** Map a chromosome name to an index position. */
        virtual ChrIndex chrNameToIndex(ChrName chrName) const = 0;

        /** chr1-chr22 */
        virtual bool isAutosome(ChrIndex index) const = 0;

        /** chrX, Y, ...*/
        virtual bool isGonosome(ChrIndex index) const = 0;

        /** phix index. */
        virtual bool isTechnical(ChrIndex index) const = 0;

        /** NC_007605, EBV. */
        virtual bool isVirus(ChrIndex index) const = 0;

        /** Mitochondrial chromosome index. */
        virtual bool isExtrachromosomal(ChrIndex index) const = 0;

        /** Decoy sequence index. */
        virtual bool isDecoy(ChrIndex index) const = 0;

        /** Chromosomes that are not assigned to a specific position in a chromosome. This
          * includes unplaced, unlocalized, and random contigs, such as GL000192.1. */
        virtual bool isUnassigned(ChrIndex index) const = 0;

        /** HLA contigs */
        virtual bool isHLA(ChrIndex index) const = 0;

        /** ALT contigs */
        virtual bool isALT(ChrIndex index) const = 0;

        // Methods for working with the subset of compressed master-ref chromosomes.

        /** Number of compressed mref chromosomes. */
        virtual CompressedMrefIndex nChromosomesCompressedMref() const = 0;

        /** Map an index position to a chromosome name for compressed mref files. */
        virtual ChrName compressedMrefIndexToChrName(CompressedMrefIndex index) const = 0;

        /** Map an index from the global index-space to the compressed mref index-space. */
        virtual CompressedMrefIndex indexToCompressedMrefIndex(ChrIndex index) const = 0;

        /** Whether the chromosome index is that of a compressed mref chromosome. */
        virtual bool isCompressedMref(ChrIndex index) const = 0;

        /** Map from compressed mref index space to all chromosome index space. */
        virtual ChrIndex compressedMrefIndexToIndex(CompressedMrefIndex index) const = 0;

        /** Map compressed mref index to chromosome size. */
        virtual ChrSize chrSizeCompressedMref(CompressedMrefIndex index) const = 0;

        /** Returns true, if the region of the read is aligned to is blocked. */
        virtual bool isInBlockedRegion(ChrIndex chrIndex, ChrSize position) const;

        /** Parse chromosome index.
          *
          * 1. Input is a plain chromosome string separated from the following string by '\t' when
          *    parsing a BED file (the chromosome identifier in the first column).
          * 2. ...
          *
          * If the `stopCharExt` parameter is an empty string, then it takes a position in a
          * character stream, and translates the following character(s) into index positions
          * (using ChrConverter::indexToChrName). If the name cannot be parsed, throws a domain_error
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
          * If no chromosome name can be parsed, throws a std::domain_error enriched with
          * boost::exception information.
          **/
        virtual ChrIndex
        parseChrAndReturnIndex(std::string::const_iterator startIt,
                               std::string::const_iterator endIt,
                               char stopChar,
                               const std::string &stopCharExt = "") const = 0;

    };

}

#endif /* CHRCONVERTER_H_ */