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

#include "ChrConverter.h"
#include "GlobalAppConfig.h"
#include "DeFuzzier.h"
#include "HelperFunctions.h"
#include "strtk.hpp"
#include <MasterRefProcessor.h>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <chrono>
#include <cmath>
#include <iostream>

namespace sophia {

    MasterRefProcessor::MasterRefProcessor(const vector<string> &filesIn,
                                           const string &outputRootName,
                                           const string &version,
                                           const int defaultReadLengthIn)
        : NUMPIDS{static_cast<int>(filesIn.size())},
          DEFAULTREADLENGTH{defaultReadLengthIn},
          mrefDb{} {

        const vector<int> &chrSizes =
            GlobalAppConfig::getInstance().getChrConverter().chrSizesCompressedMref;
        for (std::vector<int>::size_type i = 0; i < chrSizes.size(); ++i) {
            //		mrefDbPtrs.emplace_back(chrSizes[i] + 1, nullptr);
            mrefDb.emplace_back(chrSizes[i] + 1, MrefEntry{});
        }
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
        mergedBpsOutput = make_unique<ofstream>(
            outputRootName + "_" + strtk::type_to_string<int>(NUMPIDS) +
            "_mergedBpCounts.bed");
        short fileIndex{0};
        for (const auto &gzFile : filesIn) {
            chrono::time_point<chrono::steady_clock> start =
                chrono::steady_clock::now();
            auto newBreakpoints = processFile(gzFile, fileIndex);
            chrono::time_point<chrono::steady_clock> end =
                chrono::steady_clock::now();
            chrono::seconds diff =
                chrono::duration_cast<chrono::seconds>(end - start);
            ++fileIndex;
            cerr << gzFile << "\t" << diff.count() << "\t" << newBreakpoints << "\t"
                 << fileIndex << "\t" << 100 * (fileIndex + 0.0) / NUMPIDS << "%\n";
        }
        auto defuzzier = DeFuzzier{DEFAULTREADLENGTH * 3, true};
        auto i = 84;
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        while (!mrefDb.empty()) {
            mrefDb.back().erase(
                remove_if(mrefDb.back().begin(), mrefDb.back().end(),
                          [](const MrefEntry &bp) { return bp.getPos() == -1; }),
                mrefDb.back().end());
            defuzzier.deFuzzyDb(mrefDb.back());
            mrefDb.back().erase(
                remove_if(mrefDb.back().begin(), mrefDb.back().end(),
                          [](const MrefEntry &bp) { return bp.getPos() == -1; }),
                mrefDb.back().end());
            auto chromosome = chrConverter.indexToChrCompressedMref[i];
            --i;
            for (auto &bp : mrefDb.back()) {
                if (bp.getPos() != -1 && bp.getValidityScore() != -1) {
                    //				cout <<
                    //bp.printArtifactRatios(chromosome);
                    *mergedBpsOutput << bp.printBpInfo(chromosome);
                }
            }
            mrefDb.pop_back();
        }
    }

    unsigned long long
    MasterRefProcessor::processFile(const string &gzPath, short fileIndex) {
        unsigned long long newBreakpoints{0};
        ifstream refHandle(gzPath, ios_base::in | ios_base::binary);
        boost::iostreams::filtering_istream gzStream{};
        gzStream.push(boost::iostreams::gzip_decompressor());
        gzStream.push(refHandle);
        string sophiaLine{};
        vector<vector<BreakpointReduced>> fileBps{85, vector<BreakpointReduced>{}};
        auto lineIndex = 0;
        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
        while (error_terminating_getline(gzStream, sophiaLine)) {
            if (sophiaLine[0] != '#') {
                auto chrIndex =
                    chrConverter.indexConverter[chrConverter.readChromosomeIndex(
                            sophiaLine.cbegin(), '\t')];
                if (chrIndex < 0) {
                    continue;
                }
                Breakpoint tmpBp{sophiaLine, true};
                fileBps[chrIndex].emplace_back(
                    tmpBp, lineIndex++,
                    (sophiaLine.back() != '.' && sophiaLine.back() != '#'));
            }
        }
        auto chrIndex = 0;
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
    MasterRefProcessor::processBp(BreakpointReduced &bp, int chrIndex,
                                  short fileIndex) {
        MrefEntry tmpMrefEntry{};
        tmpMrefEntry.addEntry(bp, fileIndex);
        auto validitiyInit =
            mrefDb[chrIndex][tmpMrefEntry.getPos()].getValidityScore();
        mrefDb[chrIndex][tmpMrefEntry.getPos()].mergeMrefEntries(tmpMrefEntry);
        auto validitiyFinal =
            mrefDb[chrIndex][tmpMrefEntry.getPos()].getValidityScore();
        if (validitiyFinal > validitiyInit) {
            return true;
        }
        return false;
    }

}   // namespace sophia
