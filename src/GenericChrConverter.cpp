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

#include "GenericChrConverter.h"
#include "global.h"
#include "ChrInfo.h"
#include "ChrInfoTable.h"
#include "ChrCategory.h"

#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string/join.hpp>


namespace sophia {

    /* This makes a search structure to map chromosome names to index positions in an array.
       This is useful mostly during the parsing of chromosome names from the input. In most
       downstream code just the indices are used, but the parsing has to be done for every input
       line.

       unordered_map is quite fast and O(1), but we might consider making a tuned structure here,
       that takes into account the hit-probability of reads against contigs, based on the
       assumption of a uniform hit-probability proportional to chromosome size. */
    GenericChrConverter::ChrToIndexMap
    GenericChrConverter::buildAllChromosomeLookup(const ChrInfoTable::ChrNames &chromosomes) {
        ChrToIndexMap mapping;
        mapping.reserve(chromosomes.size());
        for (ChrIndex i = 0; i < (ChrIndex) chromosomes.size(); ++i) {
            mapping[chromosomes[static_cast<long>(i)]] = i;
        }
        return mapping;
    }

    std::vector<ChrIndex>
    GenericChrConverter::buildCompressedMrefToAllMapping(ChrInfoTable chrInfoIn) {
        std::vector<ChrIndex> mapping;
        mapping.reserve((size_t) chrInfoIn.nChromosomes());
        for (ChrIndex idx = 0; idx < (ChrIndex) chrInfoIn.nChromosomes(); ++idx) {
            if (chrInfoIn.getChrInfos()[static_cast<long>(idx)].isCompressedMref()) {
                mapping.emplace_back(idx);
            }
        }
        return mapping;
    }

    std::vector<std::optional<CompressedMrefIndex>>
    GenericChrConverter::buildAllToCompressedMrefMapping(ChrInfoTable chrInfoIn) {
        std::vector<std::optional<CompressedMrefIndex>> mapping;
        mapping.reserve((size_t) chrInfoIn.nChromosomes());
        CompressedMrefIndex compressedMrefIndex = 0;
        for (ChrIndex idx = 0; idx < chrInfoIn.nChromosomes(); ++idx) {
            std::optional<CompressedMrefIndex> compressedMrefIndexO = std::nullopt;
            if (chrInfoIn.getChrInfos()[static_cast<long>(idx)].isCompressedMref()) {
                compressedMrefIndexO = std::optional<CompressedMrefIndex>(compressedMrefIndex);
                ++compressedMrefIndex;
            }
            mapping.emplace_back(compressedMrefIndexO);
        }
        return mapping;
    }

    GenericChrConverter::GenericChrConverter(
        std::string assemblyNameIn,
        ChrInfoTable chrInfoTableIn)
            : ChrConverter(assemblyNameIn),
              chrInfoTable { chrInfoTableIn },
              allChromosomeLookup { buildAllChromosomeLookup(chrInfoTableIn.getNames()) },
              compressedToAllMapping { buildCompressedMrefToAllMapping(chrInfoTableIn) },
              allToCompressedMapping { buildAllToCompressedMrefMapping(chrInfoTableIn) } {}

    /** Number of all chromosomes. */
    ChrIndex GenericChrConverter::nChromosomes() const {
        return chrInfoTable.nChromosomes();
    }

    /** Map an index position to a chromosome name. */
    ChrName GenericChrConverter::indexToChrName(ChrIndex index) const {
        return chrInfoTable.getChrInfos().at(static_cast<long>(index)).getName();
    }

    /** Map a chromosome name to an index position. */
    ChrIndex
    GenericChrConverter::chrNameToIndex(ChrName chrName) const {
        ChrIndex result;
        try {
            result = allChromosomeLookup.at(chrName);
        } catch (std::out_of_range &e) {
            throw_with_trace(DomainError("Chromosome name not found: '" + chrName + "'"));
        }
        return result;
    }


    /** chr1-chr22 */
    bool GenericChrConverter::isAutosome(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::AUTOSOME;
    }

    /** chrX */
    bool GenericChrConverter::isX(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::X;
    }

    /** chrY */
    bool GenericChrConverter::isY(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::Y;
    }

    /** chrX, chrY */
    bool GenericChrConverter::isGonosome(ChrIndex index) const {
        return isX(index) || isY(index);
    }

    /** phix index. */
    bool GenericChrConverter::isTechnical(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::TECHNICAL;
    }

    /** NC_007605, EBV. */
    bool GenericChrConverter::isVirus(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::VIRUS;
    }

    /** Mitochondrial chromosome index. */
    bool GenericChrConverter::isExtrachromosomal(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::EXTRACHROMOSOMAL;
    }

    /** Decoy sequence index. */
    bool GenericChrConverter::isDecoy(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::DECOY;
    }

    /** ALT sequence index. */
    bool GenericChrConverter::isALT(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::ALT;
    }

    /** HLA sequence index. */
    bool GenericChrConverter::isHLA(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::HLA;
    }

    /** Unassigned (unplaced, random, unlocalized) sequence index. */
    bool GenericChrConverter::isUnassigned(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].getCategory() == ChrCategory::UNASSIGNED;
    }

    bool GenericChrConverter::isCompressedMref(ChrIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(index)].isCompressedMref();
    }

    /** Number of compressedMref chromosomes. */
    CompressedMrefIndex GenericChrConverter::nChromosomesCompressedMref() const {
        return CompressedMrefIndex(compressedToAllMapping.size());
    }

    /** Map the compressed mref index to the uncompressed mref index. */
    ChrIndex GenericChrConverter::compressedMrefIndexToIndex(CompressedMrefIndex compressedMrefIndex) const {
        if (compressedMrefIndex >= nChromosomesCompressedMref()) {
            throw_with_trace(std::logic_error("Compressed mref index out of range."));
        }
        ChrIndex result = compressedToAllMapping[static_cast<unsigned int>(compressedMrefIndex)];
        // The following is just a crude logic test. It will fail, if there is something wrong
        // with the index space mapping and it is rather a guard against programming errors.
        // If global and compressed mref indices are properly type-checked (instead of
        // using/typedef declarations, which are not type-checked!), then this should be removed.
        // TODO Remove when switching to typed ChrIndex and CompressedMrefIndex.
        if (!chrInfoTable.getChrInfos()[static_cast<long>(result)].isCompressedMref())
            throw_with_trace(DomainError(
                "Compressed mref index does not map back to a compressed mref chromosome."));
        return result;
    }

    /** Map an index from the global index-space to the compressed mref index-space. */
   CompressedMrefIndex
   GenericChrConverter::indexToCompressedMrefIndex(ChrIndex index) const {
        if (allToCompressedMapping.at(static_cast<long>(index)) == std::nullopt) {
            throw_with_trace(std::logic_error(
                "Index does not map to a compressed mref chromosome."));
        }
        return allToCompressedMapping.at(static_cast<long>(index)).value();
   }

    /** Map compressed mref index to chromosome size. */
    ChrSize GenericChrConverter::chrSizeCompressedMref(CompressedMrefIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(compressedMrefIndexToIndex(index))].getSize();
    }

    /** Map an compressed mref index to a chromosome name. */
    ChrName GenericChrConverter::compressedMrefIndexToChrName(CompressedMrefIndex index) const {
        return chrInfoTable.getChrInfos()[static_cast<long>(compressedMrefIndexToIndex(index))].getName();
    }


    /** Parse chromosome index. It takes a position in a character stream, and translates the
        following character(s) into index positions (using ChrConverter::indexToChrName).
        If the name cannot be parsed, throws a domain_error exception.

        This method parses up to the first occurrence of the `stopChar1`. Then within the identified
        start and end positions, parses up to the last occurrence of `stopChar2`. This allows to
        parse a chromosome name "HLA-DRB1*13:01:01" from a string
        "HLA-DRB1*13:01:01:2914|(4,0,0?/0)" by first separating out the `|` separator, and then
        finding the last `:` separator before position.
        */
    ChrName GenericChrConverter::parseChrBreakPoint(std::string::const_iterator startIt,
                                                 std::string::const_iterator endIt,
                                                 char stopChar,
                                                 const std::string &stopCharsExt) {
        if (stopCharsExt.empty()) {
            throw_with_trace(std::invalid_argument("stopCharsExt must not be empty."));
        }

        // First find the outer separator (one of `stopCharsExt`).
        auto isStopCharExt = [stopCharsExt](char c) {
            return stopCharsExt.find(c) != std::string::npos;
        };
        auto endMatchIt = std::find_if(startIt, endIt, isStopCharExt);

        // Then find the inner separator (`stopChar`) by searching backwards from the outer
        // separator.
        auto isStopChar = [stopChar](char c) { return c == stopChar; };
        auto reverseStartIt = std::reverse_iterator(endMatchIt);
        auto reverseEndMatchIt = std::find_if(reverseStartIt,
                                              std::reverse_iterator(startIt),
                                              isStopChar);
        // The reverseEndMatchIt now points onto the stopChar. We need to reverse it again (base())
        // which will *include* the stopChar in the result, which we don't want. Therefore, we
        // increment the reverseEndMatchIt once.
        ++reverseEndMatchIt;

        // Finally, prepare and return the result.
        ChrName chrName;
        chrName.reserve(50);   // Should be sufficient for most chromosome names.
        std::copy(startIt,
                  reverseEndMatchIt.base(),  // back-convert reverse_iterator to normal iterator.
                  std::back_inserter(chrName));
        return chrName;
    }

    /** Parse the chromosome index just by finding the `stopChar`. Everything between the `startIt`,
           and the first occurrence of the `stopChar` is returned as chromosome name. */
    ChrName GenericChrConverter::parseChrSimple(std::string::const_iterator startIt,
                                             std::string::const_iterator endIt,
                                             char stopChar) {
        auto isStopChar = [stopChar](char c) { return c == stopChar; };
        auto endMatchIt = std::find_if(startIt, endIt, isStopChar);

        // Prepare and return the result.
        std::string chrName;
        chrName.reserve(50);   // Should be sufficient for most chromosome names.
        std::copy(startIt,
                  endMatchIt,
                  std::back_inserter(chrName));
        return chrName;
    }

    ChrName GenericChrConverter::parseChr(std::string::const_iterator startIt,
                                       std::string::const_iterator endIt,
                                       char stopChar,
                                       const std::string &stopCharsExt) {
        if (stopCharsExt.empty()) {
            return parseChrSimple(startIt, endIt, stopChar);
        } else {
            return parseChrBreakPoint(startIt, endIt, stopChar, stopCharsExt);
        }
    }

    ChrIndex GenericChrConverter::parseChrAndReturnIndex(std::string::const_iterator startIt,
                                                         std::string::const_iterator endIt,
                                                         char stopChar,
                                                         const std::string &stopCharsExt) const {
        ChrName chrName = parseChr(startIt, endIt, stopChar, stopCharsExt);

        // Map to ChrIndex and return it, of if the chromosome is not registered, give a helpful
        // error message, that shows the parsed name and from what input it was parsed.
        try {
            return allChromosomeLookup.at(chrName);
        } catch (std::out_of_range& e) {
            throw_with_trace(DomainError(
                "Chromosome name '" + chrName + "' not found for assembly '" +
                getAssemblyName() + "'."));
        }
        // Just to get rid of a warning.
        return std::numeric_limits<ChrIndex>::max();
    }

} /* namespace sophia */
