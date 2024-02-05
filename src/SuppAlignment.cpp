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

#include "global.h"
#include "strtk-wrap.h"
#include "SuppAlignment.h"
#include <boost/exception/all.hpp>
#include <vector>
#include <algorithm>


namespace sophia {

    using namespace std;

    double SuppAlignment::ISIZEMAX { };

    ChrSize SuppAlignment::DEFAULT_READ_LENGTH { };


    // Default constructor.
    SuppAlignment::SuppAlignment() :
        matchFuzziness { 5 * DEFAULT_READ_LENGTH },
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


    SuppAlignment
    SuppAlignment::create(ChrIndex chrIndexIn,
                          ChrSize posIn,
                          int mateSupportIn,
                          int expectedDiscordantsIn,
                          bool encounteredMIn,
                          bool invertedIn,
                          ChrSize extendedPosIn,
                          bool primaryIn,
                          bool lowMapqSourceIn,
                          bool nullMapqSourceIn,
                          int originIndexIn) {
        SuppAlignment result = SuppAlignment();

        result.chrIndex = chrIndexIn;
        result.pos = posIn;
        result.extendedPos = extendedPosIn;
        result.mateSupport = mateSupportIn;
        result.expectedDiscordants = expectedDiscordantsIn;
        result.encounteredM = encounteredMIn;
        result.inverted = invertedIn;
        result.fuzzy = extendedPosIn != posIn;
        result.distant = true;
        result.lowMapqSource = lowMapqSourceIn;
        result.nullMapqSource = nullMapqSourceIn;
        result.primary = primaryIn;

        if (originIndexIn != -1) {
            if (result.primary) {
                result.supportingIndices.push_back(originIndexIn);
            } else {
                result.supportingIndicesSecondary.push_back(originIndexIn);
            }
        } else {
            result.distinctReads = 0;
        }
        result.strictFuzzy = result.fuzzy || (result.support + result.secondarySupport) < 3;

        return result;
    }

    /** Parse the supplementary alignment information from an SA:Z: tag according to the a
      * SAM specification, such as
      *
      * SA:Z:10,24753146,+,68S33M,48,1;X,135742083,-,47S22M32S,0,0;8,72637925,-,29S19M53S,0,0"
      */
    SuppAlignment SuppAlignment::parseSamSaTag(string::const_iterator saCbegin,
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

        // Split the SA tag (from a SAM file) into fields. From the SAM specification:
        //
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
            NM [[gnu::unused]] = 5;   // not used in this parser and only provided for documentation

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
        try {
            result.chrIndex = chrConverter.parseChrAndReturnIndex(
                fieldBegins[RNAME],
                fieldEnds[RNAME],
                ',');
        } catch (DomainError& e) {
            e <<
                error_info_string("from = " + string(fieldBegins[RNAME], fieldEnds[RNAME]));
            throw e;
        }

        // If the chromosome is to be ignored, don't update any of the other fields.
        if (chrConverter.isTechnical(result.chrIndex)) {
            return result;
        }
        // else

        // Update `pos` field.
        for (auto it = fieldBegins[POS]; it != fieldEnds[POS]; ++it) {
            result.pos = 10 * result.pos + ChrPosition(*it - '0');
        }

        // Update `mapq` field.
        for (auto it = fieldBegins[MAPQ]; it != fieldEnds[MAPQ]; ++it) {
            result.mapq = 10 * result.mapq + (int) (*it - '0');
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
        unsigned int chunkIndex = 0,
                     largestClipIndex = 0;
        auto cumulativeNucleotideCount = 0,
             currentNucleotideCount = 0,
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

            result.distant = (bpChrIndex != result.chrIndex || (abs((long) bpPos - (long) result.pos) > ISIZEMAX));
            if (bpChrIndex == result.chrIndex) {
                result.matchFuzziness = min((ChrSize) abs((long) bpPos - (long) result.pos), result.matchFuzziness);
            }
        }

        return result;
    }

    static const string STOP_CHARS = "|(\t";
    inline bool isStopChar(char c) {
        return STOP_CHARS.find(c) != std::string::npos;
    };

    /** The syntax parsed here is the same as the one generated by `SuppAlignment::print()`.
      *  TODO FIX: This parser parses the same syntax as `SuppAlignmentAnno::parseSaSupport(const string&)`.
      */
    SuppAlignment SuppAlignment::parseSaSupport(const string& saIn) {
        try {
            const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();

            SuppAlignment result = SuppAlignment();

            // If the last character is a `#` then properPairErrorProne is true.
            result.properPairErrorProne = saIn.back() == '#';

            unsigned int index = 0;

            // If the string starts with a `|` then encounteredM is true.
            result.encounteredM = saIn.at(0) == '|';
            if (result.encounteredM) {
                ++index;
            }

            // Parse chromosome name. The chromosome name will be separated from the position information
            // by a colon ':' character, but as the chromosome name itself may also contain colons, we
            // need to anchor first character after the position which is either a `|` or a `(`, and then
            // track back to the *last* colon.
            try {
                result.chrIndex = chrConverter.parseChrAndReturnIndex(
                    next(saIn.cbegin(), index),
                    saIn.cend(),
                    ':',
                    STOP_CHARS);
            } catch (DomainError& e) {
                e <<
                    error_info_string("from = " + string(next(saIn.cbegin(), index), saIn.cend()));
                throw e;
            }

            // If this is an ignored chromosome, don't bother parsing the rest.
            if (chrConverter.isTechnical(result.chrIndex)) {
                return result;
            }

            // else, skip forward to the first colon ':' character. This ':' will be in column 6 or 7,
            // dependent on the support information there.
            while (saIn.at(index) != ':') {
                ++index;
            }
            ++index;

            // ... and continue to parse the breakpoint specification, which gives information about the
            // position, and whether the breakpoint is fuzzy or inverted.
            for (; saIn.at(index) != '('; ++index) {
                if (saIn.at(index) == '-') {
                    result.fuzzy = true;
                } else if (saIn.at(index) == '_') {
                    result.inverted = true;
                    while (saIn.at(index) != '(') {
                        ++index;
                    }
                    break;
                } else if (saIn.at(index) != '|') {
                    if (!result.fuzzy) {
                        result.pos = 10 * result.pos + (unsigned int) (saIn.at(index) - '0');
                    } else {
                        result.extendedPos = 10 * result.extendedPos + (unsigned int) (saIn.at(index) - '0');
                    }
                }
            }

            if (!result.fuzzy) {
                result.extendedPos = result.pos;
            }
            ++index;

            for (; saIn.at(index) != ','; ++index) {
                result.support = 10 * result.support + (saIn.at(index) - '0');
            }
            ++index;

            for (; saIn.at(index) != ','; ++index) {
                result.secondarySupport = 10 * result.secondarySupport + (saIn.at(index) - '0');
            }

            ++index;
            if (saIn.at(index) == '!') {
                result.suspicious = true;
                index += 2;
            } else {
                for (; saIn.at(index) != '/'; ++index) {
                    if (saIn.at(index) == '?') {
                        result.semiSuspicious = true;
                    } else {
                        result.mateSupport = 10 * result.mateSupport + (saIn.at(index) - '0');
                    }
                }
                ++index;
            }
            for (; saIn.at(index) != ')'; ++index) {
                result.expectedDiscordants = 10 * result.expectedDiscordants + (saIn.at(index) - '0');
            }

            result.distant = result.expectedDiscordants > 0 || result.suspicious;
            bool strictFuzzyCandidate = (result.support + result.secondarySupport) < 3;
            result.strictFuzzy = result.fuzzy || strictFuzzyCandidate;

            return result;
        } catch (out_of_range& e) {
            throw boost::enable_error_info(e) <<
                error_info_string("from = " + saIn);
        }
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


    /** The syntax is as follows:
     *
     *  spec ::= encounteredM position inverted notEncounteredM
     *             '(' support ',' secondarySupport ',' mateInfo '/' expectedDiscordants ')'
     *             properPairErrorProne
     *  encounteredM ::= <epsilon> | ('|'  iff sa.encounteredM == true)
     *  position ::= chrName ':' position2
     *  position2 ::= pos | pos '-' extendedPos
     *  pos := [0-9]+
     *  extendedPos := [0-9]+
     *  inverted ::= <epsilon> | ('_INV' iff sa.inverted == true)
     *  notEncounteredM ::= <epsilon> | ('|' iff sa.encounteredM == false)
     *  support ::= [0-9]+
     *  secondarySupport ::= [0-9]+
     *  mateInfo ::= (`!` iff sa.suspicious == true) | mateSupport
     *  mateSupport ::= [0-9]+ mateAddition
     *  mateAddition ::= <epsilon> | (`?` iff sa.semiSuspicious == true or sa.nullMapqSource == true)
     *  properPairErrorProne ::= <epsilon> | (`#` iff sa.properPairErrorProne == true)
     *
     * <epsilon> here refers to the empty string ""
     **/
    string SuppAlignment::print() const {
        string outStr;
        outStr.reserve(36);

        string invStr { };
        if (inverted) {
            invStr.append("_INV");
        }

        // If the encounteredM attribute is set then the string starts with a `|`, otherwise, the
        // inv-string starts with a `|`.
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
                fuzziness = 2.5 * DEFAULT_READ_LENGTH;
                return ((long) rhs.getPos() - (long) fuzziness) <= ((long) extendedPos + (long) fuzziness) &&
                            ((long) pos - (long) fuzziness) <= ((long) rhs.getExtendedPos() + (long) fuzziness);
            } else {
                return abs((long) pos - (long) rhs.getPos()) <= (long) fuzziness;
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
                return ((long) rhs.getPos() - (long) fuzziness) <= ((long) extendedPos + (long) fuzziness) &&
                            ((long) pos - (long) fuzziness) <= ((long) rhs.getExtendedPos() + (long) fuzziness);
            } else {
                return abs((long) pos - (long) rhs.getPos()) <= (long) fuzziness;
            }
        } else {
            return false;
        }
    }

}/* namespace sophia */

