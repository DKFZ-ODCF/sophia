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

    // These two are only to make the code clearer, but are not type checked. There are no opaque
    // or strongly type-checked typedefs in C++17.
    typedef size_t ChrIndex;
    typedef size_t CompressedMrefIndex;
    typedef long unsigned int ChrSize;

    /** ChrConverter manages information on chromosomes names, sizes, and index positions in
        data arrays.

        This probably needs a redesign. The current situation is just an intermediate step away
        from the former completely procedural implementation that heavily leaked implementation
        details into the calling code and was highly tuned but very inflexible. */
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

        /** Map the compressed mref index to the uncompressed mref index. */
        virtual ChrIndex compressedMrefIndexToIndex(CompressedMrefIndex index) const = 0;

        /** Map compressed mref index to chromosome size. */
        virtual ChrSize chrSizeCompressedMref(CompressedMrefIndex index) const = 0;

        /** Map a chromosome name to an index position. */
        virtual ChrIndex chrNameToIndex(std::string chrName) const = 0;

        /** Parse chromosome index. It takes a position in a character stream, and translates the
            following character(s) into index positions (using ChrConverter::indexToChr).
            If the name cannot be parsed, throws a domain_error exception. */
        virtual ChrIndex parseChrAndReturnIndex(std::string::const_iterator startIt,
                                                char stopChar) const = 0;

    };

}

#endif /* _CHRCONVERTER_H_ */