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


namespace sophia {

using namespace std;

double SuppAlignment::ISIZEMAX { };

int SuppAlignment::DEFAULTREADLENGTH { };


// Default constructor.
SuppAlignment::SuppAlignment() :
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
    encounteredM { false },
    toRemove { false },
    inverted { false },
    fuzzy { false },
    strictFuzzy { false },
    distant { false },
    lowMapqSource { false },
    nullMapqSource { false },
    suspicious { false },
    semiSuspicious { false },
    properPairErrorProne { false },
    primary { true } {}


SuppAlignment::SuppAlignment(ChrIndex chrIndexIn,
                             int posIn,
                             int mateSupportIn,
                             int expectedDiscordantsIn,
                             bool encounteredMIn,
                             bool invertedIn,
                             int extendedPosIn,
                             bool primaryIn,
                             bool lowMapqSourceIn,
                             bool nullMapqSourceIn,
                             int originIndexIn) : SuppAlignment() {
    chrIndex = chrIndexIn;
    pos = posIn;
    extendedPos = extendedPosIn;
    mateSupport = mateSupportIn;
    expectedDiscordants = expectedDiscordantsIn;
    encounteredM = encounteredMIn;
    inverted = invertedIn;
    fuzzy = extendedPosIn != posIn;
    distant = true;
    lowMapqSource = lowMapqSourceIn;
    nullMapqSource = nullMapqSourceIn;
    primary = primaryIn;

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

SuppAlignment SuppAlignment::parse(string::const_iterator saCbegin,
                                   string::const_iterator saCend,
                                   bool primaryIn,
                                   bool lowMapqSourceIn,
                                   bool nullMapqSourceIn,
                                   bool alignmentOnForwardStrand,
                                   bool bpEncounteredM,
                                   int originIndexIn,
                                   ChrIndex bpChrIndex,
                                   int bpPos) {

    SuppAlignment result = SuppAlignment();
    result.encounteredM = bpEncounteredM;
    result.lowMapqSource = lowMapqSourceIn;
    result.nullMapqSource = nullMapqSourceIn;
    result.primary = primaryIn;

	if (result.primary) {
		result.supportingIndices.push_back(originIndexIn);
	} else {
		result.supportingIndicesSecondary.push_back(originIndexIn);
	}
//	//"SA:Z:10,24753146,+,68S33M,48,1;X,135742083,-,47S22M32S,0,0;8,72637925,-,29S19M53S,0,0;"
//	//"10,24753146,+,68S33M,48,1"
//	for (auto cigarString_cit = saCbegin; cigarString_cit != saCend; ++cigarString_cit) {
//		cerr << *cigarString_cit;
//	}
//	cerr << endl;

    // Split the SA tag into fields. From the SAM specification:
    //   SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ Other canonical alignments in a chimeric alignment, for-
    //   matted as a semicolon-delimited list. Each element in the list represents a part of the chimeric align-
    //   ment. Conventionally, at a supplementary line, the first element points to the primary line. Strand is
    //   either ‘+’ or ‘-’, indicating forward/reverse strand, corresponding to FLAG bit 0x10. Pos is a 1-based
    //   coordinate.
    //
    // NOTE: This parser does *NOT* cover the case with multiple semicolon-separated alignments.
    static const unsigned int
        RNAME = 0,
        POS = 1,
        STRAND = 2,
        CIGAR = 3,
        MAPQ = 4,
        NM = 5;

	vector<string::const_iterator> fieldBegins = { saCbegin };
	vector<string::const_iterator> fieldEnds;
	for (auto it = saCbegin; it != saCend; ++it) {
		if (*it == ',') {
			fieldEnds.push_back(it);
			fieldBegins.push_back(it + 1);
		}
	}
	fieldEnds.push_back(saCend);

    // Update `chrIndex` field.
    const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
	result.chrIndex = chrConverter.parseChrAndReturnIndex(
	    fieldBegins[RNAME],
	    fieldEnds[RNAME],
	    ',');

    // If the chromosome is to be ignored, don't update any of the other fields.
	if (chrConverter.isIgnoredChromosome(result.chrIndex)) {
		return result;
	}
	// else

    // Update `pos` field.
	for (auto it = fieldBegins[POS]; it != fieldEnds[POS]; ++it) {
		result.pos = 10 * result.pos + (*it - '0');
	}

    // Update `mapq` field.
    for (auto it = fieldBegins[MAPQ]; it != fieldEnds[MAPQ]; ++it) {
		result.mapq = 10 * result.mapq + (*it - '0');
	}

	// Update `inverted` field
	if (alignmentOnForwardStrand) {
		result.inverted = ('+' != *fieldBegins[STRAND]);
	} else {
		result.inverted = ('-' != *fieldBegins[STRAND]);
	}

	// Update `strictFuzzy` field.
	result.strictFuzzy = result.fuzzy || (result.support + result.secondarySupport) < 3;

    // Now, parse the CIGAR string and identify soft/hard-clipped segments.
	vector<CigarChunk> cigarChunks;
	auto cigarEncounteredM = false;
	auto cumulativeNucleotideCount = 0,
	     currentNucleotideCount = 0,
	     chunkIndex = 0,
	     largestClipIndex = 0,
	     indelAdjustment = 0;
	auto largestClipSize = 0;
	auto leftClipAdjustment = 0;
	for (auto cigarString_cit = fieldBegins[CIGAR];
	     cigarString_cit != fieldEnds[CIGAR];
	     ++cigarString_cit) {
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
				cigarChunks.emplace_back(
				    *cigarString_cit,
                    cigarEncounteredM,
                    cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
                    currentNucleotideCount);
				if (largestClipSize < currentNucleotideCount) {
					largestClipSize = currentNucleotideCount;
					largestClipIndex = chunkIndex;
				}
				++chunkIndex;
				cumulativeNucleotideCount += currentNucleotideCount;
				break;
			case 'H':
				cigarChunks.emplace_back(
				    *cigarString_cit,
                    cigarEncounteredM,
                    cumulativeNucleotideCount + indelAdjustment - leftClipAdjustment,
                    currentNucleotideCount);
				if (largestClipSize < currentNucleotideCount) {
					largestClipSize = currentNucleotideCount;
					largestClipIndex = chunkIndex;
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

    // NOTE: If there are no supplementary alignments with no soft/hard-clipped segments, then
    //       `cigarChunks` will be empty. In this case, we still return the SuppAlignment, even
    //       if it is located on, e.g. a decoy chromosome.
	if (cigarChunks.size() != 0) {
        // We found soft/hard-clipped segments. Update `pos` and `extendedPos` fields.
        if (cigarChunks[largestClipIndex].encounteredM) {
            result.pos += cigarChunks[largestClipIndex].startPosOnRead;
        }
        result.extendedPos = result.pos;

        result.distant = (bpChrIndex != result.chrIndex || (abs(bpPos - result.pos) > ISIZEMAX));
        if (bpChrIndex == result.chrIndex) {
            result.matchFuzziness = min(abs(bpPos - result.pos), result.matchFuzziness);
        }
    }

	return result;
}

static const string STOP_CHARS = "|(,!/?)\t";
inline bool isStopChar(char c) {
    return STOP_CHARS.find(c) != std::string::npos;
};

SuppAlignment SuppAlignment::parse(const string& saIn) {
    const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();
    SuppAlignment result = SuppAlignment();
    result.properPairErrorProne = saIn.back() == '#';
    result.encounteredM = saIn[0] == '|';

	auto index = 0;

	// If the current index is the field separator '|', skip it.
	if (result.encounteredM) {
		++index;
	}

	// Parse the first column with the chromosome name (BED format)
	result.chrIndex = chrConverter.parseChrAndReturnIndex(
	    next(saIn.cbegin(), index),
	    saIn.cend(),
	    ':');

	if (chrConverter.isIgnoredChromosome(result.chrIndex)) {
		return result;
	}

	// else, skip forward to the first colon ':' character.
	while (saIn[index] != ':') {
		++index;
	}
	++index;

    // ... and continue to parse the breakpoint specification, which gives information about the
    // position, and whether the breakpoint is fuzzy or inverted.
	for (; saIn[index] != '('; ++index) {
		if (saIn[index] == '-') {
			result.fuzzy = true;
		} else if (saIn[index] == '_') {
			result.inverted = true;
			while (saIn[index] != '(') {
				++index;
			}
			break;
		} else if (saIn[index] != '|') {
			if (!result.fuzzy) {
				result.pos = 10 * result.pos + (saIn[index] - '0');
			} else {
				result.extendedPos = 10 * result.extendedPos + (saIn[index] - '0');
			}
		}
	}

	if (!result.fuzzy) {
		result.extendedPos = result.pos;
	}
	++index;

	for (; saIn[index] != ','; ++index) {
		result.support = 10 * result.support + (saIn[index] - '0');
	}
	++index;

	for (; saIn[index] != ','; ++index) {
		result.secondarySupport = 10 * result.secondarySupport + (saIn[index] - '0');
	}

	++index;
	if (saIn[index] == '!') {
		result.suspicious = true;
		index += 2;
	} else {
		for (; saIn[index] != '/'; ++index) {
			if (saIn[index] == '?') {
				result.semiSuspicious = true;
			} else {
				result.mateSupport = 10 * result.mateSupport + (saIn[index] - '0');
			}
		}
		++index;
	}
	for (; saIn[index] != ')'; ++index) {
		result.expectedDiscordants = 10 * result.expectedDiscordants + (saIn[index] - '0');
	}

	result.distant = result.expectedDiscordants > 0 || result.suspicious;

	result.strictFuzzy = result.fuzzy || (result.support + result.secondarySupport) < 3;

	return result;
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

