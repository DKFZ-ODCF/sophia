/*
 * MasterRefProcessor.cpp
 *
 *  Created on: 27 Apr 2016
 *      Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical
 * Bioinformatics, Bioinformatics and Omics Data Analytics and currently
 * Neuroblastoma Genomics) Copyright (C) 2018 Umut H. Toprak, Matthias
 * Schlesner, Roland Eils and DKFZ Heidelberg
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
 *      LICENSE: GPL
 */

#include "GlobalAppConfig.h"
#include "DeFuzzier.h"
#include "HelperFunctions.h"
#include "strtk-wrap.h"
#include <MasterRefProcessor.h>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/exception/all.hpp>
#include <chrono>
#include <cmath>
#include <iostream>

namespace sophia {

    /**
      * This constructor has a side-effect. It reads from the filesIn and write breakpoint
      * information to
      *
      *    outputRootName + "_" + NUM_PIDS + "_mergedBpCounts.bed"
      *
      * @param filesIn            vector if input gzFile names.
      * @param outputRootName     base name/path for the output files
      * @param version            the version is matched in the gzFile name to find the realPidName.
      * @param defaultreadlength  Value for the default read length used for the DeFuzzier.
      */
    MasterRefProcessor::MasterRefProcessor(const vector<string> &filesIn,
                                           const string &outputRootName,
                                           const string &version,
                                           const ChrSize defaultReadLengthIn)
        : NUM_PIDS { static_cast<int>(filesIn.size()) },
          DEFAULT_READ_LENGTH{ defaultReadLengthIn },
          mrefDb {} {

        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();

        // Preallocate the full memory in one go. Otherwise, the vector will be repeatedly
        // copied and reallocated, which is slow. Also the finally reserved size will be
        // larger than necessary.
        // NOTE: This will allocate a lot of memory as the total size of the vectors is the
        //       genome size (3.7 giga-bases for hg19).
        std::vector<std::vector<sophia::MrefEntry>>::size_type totalSize = 0;
        for (CompressedMrefIndex i = 0; i < chrConverter.nChromosomesCompressedMref(); ++i) {
            totalSize += static_cast<std::vector<std::vector<sophia::MrefEntry>>::size_type>(
                chrConverter.chrSizeCompressedMref(i) + 1);
        }
        cerr << "Allocating " << (totalSize / 1024 / 1024) << " MB for mrefDb ..." << endl;
        mrefDb.reserve(totalSize);

        // Initialize the mrefDb with default values.
        for (CompressedMrefIndex i = 0; i < chrConverter.nChromosomesCompressedMref(); ++i) {
            // It is unclear, why here +1 is added to the chromosomes sizes, in particular, as
            // the original chromosome size data already had sized incremented by 1 summing up
            // to total genome size here of N*2 additional positions, with N being the number of
            // compressed master-ref chromosomes. This is just kept for all changes :|
            mrefDb.emplace_back(chrConverter.chrSizeCompressedMref(i) + 1, MrefEntry{});
        }

        // Construct the output file header. This collects the `realPidName`s from the gzFile
        // and appends them to the header. The `version` is matched in the gzFile name.
        vector<string> header {"#chr", "start", "end"};
        for (const auto &gzFile : filesIn) {
            signed int posOnVersion = version.size() - 1;
            bool counting { false };
            string realPidName;
            for (auto rit = gzFile.crbegin(); rit != gzFile.crend(); ++rit) {
                if (!counting) {
                    // Match the version in the gzFile name. Note that we traverse the gzFile name
                    // in reverse order (from end), and therefore the matching algorithm is
                    // formulated in reverse order.
                    if (*rit != version[static_cast<unsigned long>(posOnVersion)]) {
                        // No match, means continue searching.
                        posOnVersion = version.size() - 1;
                    } else {
                        // Match, means continue matching.
                        --posOnVersion;
                        if (posOnVersion == -1) {
                            // We have a match, therefore continue in the other branch that
                            // collects the letters for the realPidName.
                            ++rit;
                            counting = true;
                        }
                    }
                } else {
                    if (*rit == '/') {
                        // We matched a '/' character, i.e. a path separator. Therefore, we are done
                        // with the realPidName. Stop collecting letters.
                        break;
                    } else {
                        // Everything we see is part of the realPidName.
                        realPidName.push_back(*rit);
                    }
                }
            }
            if (realPidName.size() == 0) {
                throw_with_trace(runtime_error(
                                    "Could not match realPidName in gzFile '" + gzFile + "'. "
                                    "The version value '" + version + "' has to be contained "
                                    "in the gzFile name. Rename the gzFile to match the pattern "
                                    "'.*/$realPidName?$version.+'."));
            }

            reverse(realPidName.begin(), realPidName.end());
            cerr << "Matched realPidName '" << realPidName
                 << "' in gzFile '" << gzFile << "'" << endl;
            header.push_back(realPidName);
        }

        // This reopens the gzFiles from filesIn and processes them. It logs information to the
        // standard error output.
        short fileIndex{0};
        for (const auto &gzFile : filesIn) {
            chrono::time_point<chrono::steady_clock> start =
                chrono::steady_clock::now();
            // newBreakpoints contains only information for chromosomes from the compressedMref set.
            auto newBreakpoints = processFile(gzFile, fileIndex);
            chrono::time_point<chrono::steady_clock> end = chrono::steady_clock::now();
            chrono::seconds diff = chrono::duration_cast<chrono::seconds>(end - start);
            ++fileIndex;
            cerr << gzFile << "\t" << diff.count() << "\t" << newBreakpoints << "\t"
                 << fileIndex << "\t" << 100 * (fileIndex + 0.0) / NUM_PIDS << "%\n";
        }

        // Finally, open the output file, and write the header and the breakpoint information.
        mergedBpsOutput = make_unique<ofstream>(
            outputRootName + "_" + strtk::type_to_string<int>(NUM_PIDS) +
            "_mergedBpCounts.bed");
        auto defuzzier = DeFuzzier {DEFAULT_READ_LENGTH * 3, true};
        CompressedMrefIndex compressedMrefChrIndex = chrConverter.nChromosomesCompressedMref() - 1;
        while (!mrefDb.empty()) {
            // Remove all invalid breakpoints.
            mrefDb.back().erase(
                remove_if(mrefDb.back().begin(), mrefDb.back().end(),
                          [](const MrefEntry &bp) { return !bp.isValid(); }),
                mrefDb.back().end());

            // Run the DeFuzzier.
            defuzzier.deFuzzyDb(mrefDb.back());

            // Again, remove all invalid breakpoints.
            mrefDb.back().erase(
                remove_if(mrefDb.back().begin(), mrefDb.back().end(),
                          [](const MrefEntry &bp) { return !bp.isValid(); }),
                mrefDb.back().end());

            // Write the breakpoint information.
            std::string chromosome =
                chrConverter.compressedMrefIndexToChrName(compressedMrefChrIndex);
            --compressedMrefChrIndex;
            for (auto &bp : mrefDb.back()) {
                if (bp.isValid()) {
                    // cout << bp.printArtifactRatios(chromosome);
                    *mergedBpsOutput << bp.printBpInfo(chromosome);
                }
            }
            mrefDb.pop_back();
        }
    }

    /**
      * Process the file at gzPath. Chromosomes not in the compressedMref set are ignored.
      * The file format the one produced by the `sophia` tool.
      */
    unsigned long long
    MasterRefProcessor::processFile(const string &gzPath, short fileIndex) {
        cerr << "Processing file '" << gzPath << "'" << endl;
        unsigned long long newBreakpoints{0};
        ifstream refHandle(gzPath, ios_base::in | ios_base::binary);
        boost::iostreams::filtering_istream gzStream{};
        gzStream.push(boost::iostreams::gzip_decompressor());
        gzStream.push(refHandle);
        string sophiaLine{};

        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        CompressedMrefIndex vectorSize = chrConverter.nChromosomesCompressedMref();

        vector<vector<BreakpointReduced>> fileBps = vector<vector<BreakpointReduced>>
            { static_cast<unsigned int>(vectorSize), vector<BreakpointReduced>{}};
        auto lineIndex = 0;

        while (error_terminating_getline(gzStream, sophiaLine)) {
            // Ignore comment lines.
            if (sophiaLine[0] != '#') {
                // Parse the chromosome name in the first column of the gzip file.
                ChrIndex globalIndex;
                try {
                    globalIndex = chrConverter.parseChrAndReturnIndex(
                        sophiaLine.cbegin(), sophiaLine.cend(), '\t');
                } catch (const DomainError &e) {
                    e <<
                        error_info_string("file = " + gzPath + ", line = " + sophiaLine);
                    throw e;
                }

                // Ignore chromosomes not in the compressedMref set.
                if (chrConverter.isCompressedMref(globalIndex)) {
                    CompressedMrefIndex compressedMrefIndex =
                        chrConverter.indexToCompressedMrefIndex(globalIndex);
                    Breakpoint tmpBp = Breakpoint::parse(sophiaLine, true);
                    fileBps[static_cast<unsigned int>(compressedMrefIndex)].emplace_back(
                        tmpBp, lineIndex++,
                        (sophiaLine.back() != '.' && sophiaLine.back() != '#'));
                }
            }
        }

        ChrIndex chrIndex = 0;
        for (auto &chromosome : fileBps) {
            DeFuzzier deFuzzierControl {DEFAULT_READ_LENGTH * 6, false};
            deFuzzierControl.deFuzzyDb(chromosome);
            for (auto &bp : chromosome) {
                if (processBp(bp, chrIndex, fileIndex)) {
                    ++newBreakpoints;
                }
            }
            ++chrIndex;
        }
        return newBreakpoints;
    }

    bool
    MasterRefProcessor::processBp(BreakpointReduced &bp,
                                  ChrIndex chrIndex,
                                  short fileIndex) {
        MrefEntry tmpMrefEntry{};
        tmpMrefEntry.addEntry(bp, fileIndex);

        unsigned long pos = static_cast<unsigned long>(tmpMrefEntry.getPos());
        unsigned int idx = static_cast<unsigned int>(chrIndex);

        auto validityInit = mrefDb[idx][pos].getValidityScore();
        mrefDb[idx][pos].mergeMrefEntries(tmpMrefEntry);
        auto validityFinal = mrefDb[idx][pos].getValidityScore();

        return validityFinal > validityInit;
    }

}   // namespace sophia
