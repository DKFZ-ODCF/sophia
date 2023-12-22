/*
 * SuppAlignment.cpp
 *
 *  Created on: 16 Apr 2016
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

#include "SuppAlignment.h"
#include <vector>
#include <algorithm>
#include "strtk-wrap.h"
#include "GlobalAppConfig.h"
//#include <iostream>

namespace sophia {

using namespace std;

double SuppAlignment::ISIZEMAX { };

int SuppAlignment::DEFAULTREADLENGTH { };

SuppAlignment::SuppAlignment(string::const_iterator saCbegin,
                             string::const_iterator saCend,
                             bool primaryIn,
                             bool lowMapqSourceIn,
                             bool nullMapqSourceIn,
                             bool alignmentOnForwardStrand,
                             bool bpEncounteredM,
                             int originIndexIn,
                             int bpChrIndex,
                             int bpPos) :
				matchFuzziness { 5 * DEFAULTREADLENGTH },
				chrIndex { 0 },
				pos { 0 },
				extendedPos { 0 },
				mapq { 0 },
				supportingIndices { },
				supportingIndicesSecondary { },
				distinctReads { 1 },
				support { 0 },
				secondarySupport { 0 },
				mateSupport { 0 },
				expectedDiscordants { 0 },
				encounteredM { bpEncounteredM },
				toRemove { false },
				inverted { false },
				fuzzy { false },
				strictFuzzy { false },
				distant { false },
				lowMapqSource { lowMapqSourceIn },
				nullMapqSource { nullMapqSourceIn },
				suspicious { false },
				semiSuspicious { false },
				properPairErrorProne { false },
				primary { primaryIn } {
	if (primary) {
		supportingIndices.push_back(originIndexIn);
	} else {
		supportingIndicesSecondary.push_back(originIndexIn);
	}
//	//"SA:Z:10,24753146,+,68S33M,48,1;X,135742083,-,47S22M32S,0,0;8,72637925,-,29S19M53S,0,0;"
//	//"10,24753146,+,68S33M,48,1"
//	for (auto cigarString_cit = saCbegin; cigarString_cit != saCend; ++cigarString_cit) {
//		cerr << *cigarString_cit;
//	}
//	cerr << endl;

	vector<string::const_iterator> fieldBegins = { saCbegin };
	vector<string::const_iterator> fieldEnds;
	for (auto it = saCbegin; it != saCend; ++it) {
		if (*it == ',') {
			fieldEnds.push_back(it);
			fieldBegins.push_back(it + 1);
		}
	}
	fieldEnds.push_back(saCend);

    const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
	chrIndex = chrConverter.parseChrAndReturnIndex(fieldBegins[0], ',');
	if (chrConverter.isIgnoredChromosome(chrIndex)) {
		return;
	}
	for (auto it = fieldBegins[1]; it != fieldEnds[1]; ++it) {
		pos = 10 * pos + (*it - '0');
	}

    //cerr << "guessSupplementOffset" << endl;
	vector<CigarChunk> cigarChunks;
	auto cigarEncounteredM = false;
	auto cumulativeNucleotideCount = 0, currentNucleotideCount = 0, chunkIndex = 0, bestChunkIndex = 0, indelAdjustment = 0;
	auto largestClip = 0;
	auto leftClipAdjustment = 0;
	for (auto cigarString_cit = fieldBegins[3]; cigarString_cit != fieldEnds[3]; ++cigarString_cit) {
		if (isdigit(*cigarString_cit)) {
			currentNucleotideCount = currentNucleotideCount * 10 + (*cigarString_cit - '0');
		} else {
			switch (*cigarString_cit) {
			case 'M':
				cigarEncounteredM = true;
				cumulativeNucleotideCount += currentNucleotideCount;
				break;
			case 'S':
				if (!cigarEncounteredM) {
					leftClipAdjustment = currentNucleotideCount;
				}
				cigarChunks.emplace_back(*cigarString_cit,
				                         cigarEncounteredM,
				                         cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
				                         currentNucleotideCount);
				if (largestClip < currentNucleotideCount) {
					largestClip = currentNucleotideCount;
					bestChunkIndex = chunkIndex;
				}
				++chunkIndex;
				cumulativeNucleotideCount += currentNucleotideCount;
				break;
			case 'H':
				cigarChunks.emplace_back(*cigarString_cit,
				                         cigarEncounteredM,
				                         cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
				                         currentNucleotideCount);
				if (largestClip < currentNucleotideCount) {
					largestClip = currentNucleotideCount;
					bestChunkIndex = chunkIndex;
				}
				++chunkIndex;
				break;
			case 'I':
				indelAdjustment -= currentNucleotideCount;
				cumulativeNucleotideCount += currentNucleotideCount;
				break;
			case 'D':
				indelAdjustment += currentNucleotideCount;
				break;
			default:
				break;
			}
			currentNucleotideCount = 0;
		}
	}
	if (cigarChunks[bestChunkIndex].encounteredM) {
		pos += cigarChunks[bestChunkIndex].startPosOnRead;
	}
	extendedPos = pos;
	//cerr << "done" << endl;
	for (auto it = fieldBegins[4]; it != fieldEnds[4]; ++it) {
		mapq = 10 * mapq + (*it - '0');
	}
	if (alignmentOnForwardStrand) {
		inverted = ('+' != *fieldBegins[2]);
	} else {
		inverted = ('-' != *fieldBegins[2]);
	}
	distant = (bpChrIndex != chrIndex || (abs(bpPos - pos) > ISIZEMAX));
	if (bpChrIndex == chrIndex) {
		matchFuzziness = min(abs(bpPos - pos), matchFuzziness);
	}
	strictFuzzy = fuzzy || (support + secondarySupport) < 3;
}

void SuppAlignment::finalizeSupportingIndices() {
	sort(supportingIndices.begin(), supportingIndices.end());
	sort(supportingIndicesSecondary.begin(), supportingIndicesSecondary.end());
	supportingIndices.erase(unique(supportingIndices.begin(), supportingIndices.end()), supportingIndices.end());
	supportingIndicesSecondary.erase(unique(supportingIndicesSecondary.begin(),
	                                        supportingIndicesSecondary.end()),
	                                        supportingIndicesSecondary.end());
	support = static_cast<int>(supportingIndices.size());
	secondarySupport = static_cast<int>(supportingIndicesSecondary.size());
}

SuppAlignment::SuppAlignment(int chrIndexIn,
                             int posIn,
                             int mateSupportIn,
                             int expectedDiscordantsIn,
                             bool encounteredMIn,
                             bool invertedIn,
                             int extendedPosIn,
                             bool primaryIn,
                             bool lowMapqSourceIn,
                             bool nullMapqSourceIn,
                             int originIndexIn) :
				matchFuzziness { 5 * DEFAULTREADLENGTH },
				chrIndex { chrIndexIn },
				pos { posIn },
				extendedPos { extendedPosIn },
				mapq { 0 },
				supportingIndices { },
				supportingIndicesSecondary { },
				distinctReads { 1 },
				support { 0 },
				secondarySupport { 0 },
				mateSupport { mateSupportIn },
				expectedDiscordants { expectedDiscordantsIn },
				encounteredM { encounteredMIn },
				toRemove { false },
				inverted { invertedIn },
				fuzzy { extendedPosIn != posIn },
				strictFuzzy { false },
				distant { true },
				lowMapqSource { lowMapqSourceIn },
				nullMapqSource { nullMapqSourceIn },
				suspicious { false },
				semiSuspicious { false },
				properPairErrorProne { false },
				primary { primaryIn } {
	if (originIndexIn != -1) {
		if (primary) {
			supportingIndices.push_back(originIndexIn);
		} else {
			supportingIndicesSecondary.push_back(originIndexIn);
		}
	} else {
		distinctReads = 0;
	}
	strictFuzzy = fuzzy || (support + secondarySupport) < 3;
}

string SuppAlignment::print() const {
	string outStr;
	outStr.reserve(36);
	string invStr { };
	if (inverted) {
		invStr.append("_INV");
	}
	if (encounteredM) {
		outStr.append("|");
	} else {
		invStr.append("|");
	}
	const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
	if (!fuzzy) {
		outStr.
		    append(chrConverter.indexToChrName(chrIndex)).
		    append(":").
		    append(strtk::type_to_string<int>(pos));
	} else {
		outStr.
		    append(chrConverter.indexToChrName(chrIndex)).
		    append(":").
		    append(strtk::type_to_string<int>(pos)).
		    append("-").
		    append(strtk::type_to_string<int>(extendedPos));
	}
	outStr.
	    append(invStr).
	    append("(").
	    append(strtk::type_to_string<int>(support)).
	    append(",").
	    append(strtk::type_to_string<int>(secondarySupport)).
	    append(",");
	if (!suspicious) {
		outStr.append(strtk::type_to_string<int>(mateSupport));
		if (semiSuspicious || nullMapqSource) {
			outStr.append("?");
		}
	} else {
		outStr.append("!");
	}
	outStr.append("/").append(strtk::type_to_string<int>(expectedDiscordants)).append(")");
	if (properPairErrorProne) {
		outStr.append("#");
	}
	return outStr;
}

SuppAlignment::SuppAlignment(const string& saIn) :
				matchFuzziness { 5 * DEFAULTREADLENGTH },
				chrIndex { 0 },
				pos { 0 },
				extendedPos { 0 },
				mapq { 0 },
				distinctReads { 1 },
				support { 0 },
				secondarySupport { 0 },
				mateSupport { 0 },
				expectedDiscordants { 0 },
				encounteredM { saIn[0] == '|' },
				toRemove { false },
				inverted { false },
				fuzzy { false },
				strictFuzzy { false },
				distant { false },
				lowMapqSource { false },
				nullMapqSource { false },
				suspicious { false },
				semiSuspicious { false },
				properPairErrorProne { saIn.back() == '#' },
				primary { true } {
	auto index = 0;
	if (encounteredM) {
		++index;
	}
	const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
	chrIndex = chrConverter.parseChrAndReturnIndex(next(saIn.cbegin(), index), ':');
	if (chrConverter.isIgnoredChromosome(chrIndex)) {
		return;
	}
	while (saIn[index] != ':') {
		++index;
	}
	++index;
	for (; saIn[index] != '('; ++index) {
		if (saIn[index] == '-') {
			fuzzy = true;
		} else if (saIn[index] == '_') {
			inverted = true;
			while (saIn[index] != '(') {
				++index;
			}
			break;
		} else if (saIn[index] != '|') {
			if (!fuzzy) {
				pos = 10 * pos + (saIn[index] - '0');
			} else {
				extendedPos = 10 * extendedPos + (saIn[index] - '0');
			}
		}
	}
	if (!fuzzy) {
		extendedPos = pos;
	}
	++index;
	for (; saIn[index] != ','; ++index) {
		support = 10 * support + (saIn[index] - '0');
	}
	++index;
	for (; saIn[index] != ','; ++index) {
		secondarySupport = 10 * secondarySupport + (saIn[index] - '0');
	}
	++index;
	if (saIn[index] == '!') {
		suspicious = true;
		index += 2;
	} else {
		for (; saIn[index] != '/'; ++index) {
			if (saIn[index] == '?') {
				semiSuspicious = true;
			} else {
				mateSupport = 10 * mateSupport + (saIn[index] - '0');
			}
		}
		++index;
	}
	for (; saIn[index] != ')'; ++index) {
		expectedDiscordants = 10 * expectedDiscordants + (saIn[index] - '0');
	}
	distant = expectedDiscordants > 0 || suspicious;
	strictFuzzy = fuzzy || (support + secondarySupport) < 3;
}

bool SuppAlignment::saCloseness(const SuppAlignment& rhs, int fuzziness) const {
	if (inverted == rhs.isInverted() && chrIndex == rhs.getChrIndex() && encounteredM == rhs.isEncounteredM()) {
		if (strictFuzzy || rhs.isStrictFuzzy()) {
			fuzziness = 2.5 * DEFAULTREADLENGTH;
			return (rhs.getPos() - fuzziness) <= (extendedPos + fuzziness) &&
			            (pos - fuzziness) <= (rhs.getExtendedPos() + fuzziness);
		} else {
			return abs(pos - rhs.getPos()) <= fuzziness;
		}
	} else {
		return false;
	}
}


bool SuppAlignment::saDistHomologyRescueCloseness(const SuppAlignment& rhs, int fuzziness) const {
	if (!distant || !rhs.isDistant()) {
		return false;
	}
	if (chrIndex == rhs.getChrIndex() && encounteredM == rhs.isEncounteredM()) {
		if (strictFuzzy || rhs.isStrictFuzzy()) {
			return (rhs.getPos() - fuzziness) <= (extendedPos + fuzziness) &&
			            (pos - fuzziness) <= (rhs.getExtendedPos() + fuzziness);
		} else {
			return abs(pos - rhs.getPos()) <= fuzziness;
		}
	} else {
		return false;
	}
}

}/* namespace sophia */

