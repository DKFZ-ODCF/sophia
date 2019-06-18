/*
 * MasterRefProcessor.cpp
 *
 *  Created on: 27 Apr 2016
 *      Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical Bioinformatics, Bioinformatics and Omics Data Analytics and currently Neuroblastoma Genomics)
 *      Copyright (C) 2018 Umut H. Toprak, Matthias Schlesner, Roland Eils and DKFZ Heidelberg
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

#include <MasterRefProcessor.h>
#include "strtk.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include "DeFuzzier.h"
#include "ChrConverter.h"
#include <boost/algorithm/string/join.hpp>
#include <chrono>
#include "HelperFunctions.h"

namespace sophia {
MasterRefProcessor::MasterRefProcessor(const std::vector<std::string>& filesIn, const std::string &outputRootName, const std::string &version, const int defaultReadLengthIn) :
				NUMPIDS { static_cast<int>(filesIn.size()) },
				DEFAULTREADLENGTH { defaultReadLengthIn },
				mrefDb { } {
	const std::vector<int> CHRSIZES { 249250622, 243199374, 198022431, 191154277, 180915261, 171115068, 159138664, 146364023, 141213432, 135534748, 135006517, 133851896, 115169879, 107349541, 102531393, 90354754, 81195211, 78077249, 59128984, 63025521, 48129896, 51304567, 155270561, 59373567, 106434, 547497, 189790, 191470, 182897, 38915, 37176, 90086, 169875, 187036, 36149, 40104, 37499, 81311,
			174589, 41002, 4263, 92690, 159170, 27683, 166567, 186859, 164240, 137719, 172546, 172295, 172150, 161148, 179199, 161803, 155398, 186862, 180456, 179694, 211174, 15009, 128375, 129121, 19914, 43692, 27387, 40653, 45942, 40532, 34475, 41935, 45868, 39940, 33825, 41934, 42153, 43524, 43342, 39930, 36652, 38155, 36423, 39787, 38503, 35477944, 171824 };
	for (auto i = 0; i < 85; ++i) {
//		mrefDbPtrs.emplace_back(CHRSIZES[i] + 1, nullptr);
		mrefDb.emplace_back(CHRSIZES[i] + 1, MrefEntry { });
	}
	std::vector<std::string> header { "#chr", "start", "end" };
	for (const auto& gzFile : filesIn) {
		int posOnVersion = version.size() - 1;
		bool counting { false };
		std::string realPidName;
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
		std::reverse(realPidName.begin(), realPidName.end());
		header.push_back(realPidName);
	}
	mergedBpsOutput = std::make_unique<std::ofstream>(outputRootName + "_" + strtk::type_to_string<int>(NUMPIDS) + "_mergedBpCounts.bed");
	short fileIndex { 0 };
	for (const auto& gzFile : filesIn) {
		std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
		auto newBreakpoints = processFile(gzFile, fileIndex);
		std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
		std::chrono::seconds diff = std::chrono::duration_cast<std::chrono::seconds>(end - start);
		++fileIndex;
		std::cerr << gzFile << "\t" << diff.count() << "\t" << newBreakpoints << "\t" << fileIndex << "\t" << 100 * (fileIndex + 0.0) / NUMPIDS << "%\n";
	}
	auto defuzzier = DeFuzzier { DEFAULTREADLENGTH * 3, true };
	auto i = 84;
	while (!mrefDb.empty()) {
		mrefDb.back().erase(std::remove_if(mrefDb.back().begin(), mrefDb.back().end(), [](const MrefEntry& bp) {return bp.getPos()==-1;}), mrefDb.back().end());
		defuzzier.deFuzzyDb(mrefDb.back());
		mrefDb.back().erase(std::remove_if(mrefDb.back().begin(), mrefDb.back().end(), [](const MrefEntry& bp) {return bp.getPos()==-1;}), mrefDb.back().end());
		auto chromosome = ChrConverter::indexToChrCompressedMref[i];
		--i;
		for (auto &bp : mrefDb.back()) {
			if (bp.getPos() != -1 && bp.getValidityScore() != -1) {
				//				std::cout << bp.printArtifactRatios(chromosome);
				*mergedBpsOutput << bp.printBpInfo(chromosome);
			}
		}
		mrefDb.pop_back();
	}
}
unsigned long long MasterRefProcessor::processFile(const std::string& gzPath, short fileIndex) {
	unsigned long long newBreakpoints { 0 };
	std::ifstream refHandle(gzPath, std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_istream gzStream { };
	gzStream.push(boost::iostreams::gzip_decompressor());
	gzStream.push(refHandle);
	std::string sophiaLine { };
	std::vector<std::vector<BreakpointReduced>> fileBps { 85, std::vector<BreakpointReduced> { } };
	auto lineIndex = 0;
	while (error_terminating_getline(gzStream, sophiaLine)) {
		if (sophiaLine[0] != '#') {
			auto chrIndex = ChrConverter::indexConverter[ChrConverter::readChromosomeIndex(sophiaLine.cbegin(), '\t')];
			if (chrIndex < 0) {
				continue;
			}
			Breakpoint tmpBp { sophiaLine, true };
			fileBps[chrIndex].emplace_back(tmpBp, lineIndex++, (sophiaLine.back() != '.' && sophiaLine.back() != '#'));
		}
	}
	auto chrIndex = 0;
	for (auto &chromosome : fileBps) {
		DeFuzzier deFuzzierControl { DEFAULTREADLENGTH * 6, false };
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

bool MasterRefProcessor::processBp(BreakpointReduced &bp, int chrIndex, short fileIndex) {
	MrefEntry tmpMrefEntry { };
	tmpMrefEntry.addEntry(bp, fileIndex);
	auto validitiyInit = mrefDb[chrIndex][tmpMrefEntry.getPos()].getValidityScore();
	mrefDb[chrIndex][tmpMrefEntry.getPos()].mergeMrefEntries(tmpMrefEntry);
	auto validitiyFinal = mrefDb[chrIndex][tmpMrefEntry.getPos()].getValidityScore();
	if (validitiyFinal > validitiyInit) {
		return true;
	}
	return false;
}

} /* namespace sophiaMref */

