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
#include <chrono>
#include <cmath>
#include <iostream>

namespace sophia {

    /**
      * This constructor has a sideeffect. It reads from the filesIn and write breakpoint
      * information to
      *
      *    outputRootName + "_" + NUMPIDS + "_mergedBpCounts.bed"
      *
      * @param filesIn            vector if input file names.
      * @param outputRootName     base name/path for the output files
      * @param version            ??
      * @param defaultReadLength  Value for the default read length used for the DeFuzzier.
      */
    MasterRefProcessor::MasterRefProcessor(const vector<string> &filesIn,
                                           const string &outputRootName,
                                           const string &version,
                                           const int defaultReadLengthIn)
        : NUMPIDS{static_cast<int>(filesIn.size())},
          DEFAULTREADLENGTH{defaultReadLengthIn},
          mrefDb{} {

        // Initialize the mrefDb with default values. Only for compressed Mref indices.
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        for (std::vector<int>::size_type i = 0;
             i < chrConverter.nChromosomesCompressedMref();
             ++i) {
            mrefDb.emplace_back(chrConverter.chrSizeCompressedMref(i) + 1, MrefEntry{});
        }

        // Construct the output file header. This collects the `realPidName`s from the gzFile
        // and appends them to the header.
        vector<string> header{"#chr", "start", "end"};
        for (const auto &gzFile : filesIn) {
            int posOnVersion = version.size() - 1;
            bool counting{false};
            string realPidName;
            for (auto rit = gzFile.crbegin(); rit != gzFile.crend(); ++rit) {
                if (!counting) {
                    if (*rit != version[posOnVersion]) {
                        posOnVersion = version.size() - 1;
                    } else {
                        --posOnVersion;
                        if (posOnVersion == -1) {
                            ++rit;
                            counting = true;
                        }
                    }
                } else {
                    if (*rit == '/') {
                        break;
                    } else {
                        realPidName.push_back(*rit);
                    }
                }
            }
            reverse(realPidName.begin(), realPidName.end());
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
            chrono::time_point<chrono::steady_clock> end =
                chrono::steady_clock::now();
            chrono::seconds diff =
                chrono::duration_cast<chrono::seconds>(end - start);
            ++fileIndex;
            cerr << gzFile << "\t" << diff.count() << "\t" << newBreakpoints << "\t"
                 << fileIndex << "\t" << 100 * (fileIndex + 0.0) / NUMPIDS << "%\n";
        }

        // Finally, open the output file, and write the header and the breakpoint information.
        mergedBpsOutput = make_unique<ofstream>(
            outputRootName + "_" + strtk::type_to_string<int>(NUMPIDS) +
            "_mergedBpCounts.bed");
        auto defuzzier = DeFuzzier{DEFAULTREADLENGTH * 3, true};
        ChrIndex compressedMrefChrIndex = chrConverter.nChromosomesCompressedMref() - 1;
        while (!mrefDb.empty()) {
            // Remove all breakpoints with position -1.
            mrefDb.back().erase(
                remove_if(mrefDb.back().begin(), mrefDb.back().end(),
                          [](const MrefEntry &bp) { return bp.getPos() == -1; }),
                mrefDb.back().end());

            // Run the DeFuzzier.
            defuzzier.deFuzzyDb(mrefDb.back());

            // Again, remove all breakpoints with position -1.
            mrefDb.back().erase(
                remove_if(mrefDb.back().begin(), mrefDb.back().end(),
                          [](const MrefEntry &bp) { return bp.getPos() == -1; }),
                mrefDb.back().end());

            // Write the breakpoint information.
            std::string chromosome =
                chrConverter.indexToChrNameCompressedMref(compressedMrefChrIndex);
            --compressedMrefChrIndex;
            for (auto &bp : mrefDb.back()) {
                if (bp.getPos() != -1 && bp.getValidityScore() != -1) {
                    // cout << bp.printArtifactRatios(chromosome);
                    *mergedBpsOutput << bp.printBpInfo(chromosome);
                }
            }
            mrefDb.pop_back();
        }
    }

    /**
      * Process the file at gzPath. Chromosomes not in the compressedMref set are ignored.
      */
    unsigned long long
    MasterRefProcessor::processFile(const string &gzPath, short fileIndex) {
        unsigned long long newBreakpoints{0};
        ifstream refHandle(gzPath, ios_base::in | ios_base::binary);
        boost::iostreams::filtering_istream gzStream{};
        gzStream.push(boost::iostreams::gzip_decompressor());
        gzStream.push(refHandle);
        string sophiaLine{};

        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();

        vector<vector<BreakpointReduced>> fileBps = vector<vector<BreakpointReduced>>
            {chrConverter.nChromosomesCompressedMref(), vector<BreakpointReduced>{}};
        auto lineIndex = 0;

        while (error_terminating_getline(gzStream, sophiaLine)) {
            // Ignore comment lines.
            if (sophiaLine[0] != '#') {
                auto chrIndexO =
                    chrConverter.compressedMrefIndexToIndex(
                        chrConverter.parseChrAndReturnIndex(
                            sophiaLine.cbegin(),
                            sophiaLine.cend(),
                            '\t'));
                // Ignore chromosomes not in the compressedMref set.
                if (chrIndexO.has_value()) {
                    ChrIndex chrIndex = chrIndexO.value();
                    Breakpoint tmpBp{sophiaLine, true};
                    fileBps[chrIndex].emplace_back(
                        tmpBp, lineIndex++,
                        (sophiaLine.back() != '.' && sophiaLine.back() != '#'));
                }
            }
        }

        ChrIndex chrIndex = 0;
        for (auto &chromosome : fileBps) {
            DeFuzzier deFuzzierControl{DEFAULTREADLENGTH * 6, false};
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
        auto validityInit =
            mrefDb[chrIndex][tmpMrefEntry.getPos()].getValidityScore();
        mrefDb[chrIndex][tmpMrefEntry.getPos()].mergeMrefEntries(tmpMrefEntry);
        auto validityFinal =
            mrefDb[chrIndex][tmpMrefEntry.getPos()].getValidityScore();
        return validityFinal > validityInit;
    }

}   // namespace sophia
