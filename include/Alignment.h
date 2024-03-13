/*
 * Alignment.h
 *
 *  Created on: 16 Apr 2016
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

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_
#include "ChosenBp.h"
#include "CigarChunk.h"
#include "CoverageAtBase.h"
#include "SuppAlignment.h"
#include "global.h"
#include <OverhangRange.h>
#include <algorithm>
#include <deque>
#include <memory>
#include <string>
#include <vector>

namespace sophia {

    using namespace std;

    class Alignment {

      public:
        Alignment();

        void continueConstruction();

        static ChrSize LOW_QUAL_CLIP_THRESHOLD;

        static int
            BASE_QUALITY_THRESHOLD,
            BASE_QUALITY_THRESHOLD_LOW;

        static ChrSize CLIPPED_NUCLEOTIDE_COUNT_THRESHOLD,
                       INDEL_NUCLEOTIDE_COUNT_THRESHOLD;

        static double ISIZEMAX;

        ChrSize getStartPos() const { return startPos; }

        ChrSize getEndPos() const { return endPos; }

        int getReadType() const { return readType; }

        const vector<ChrSize> &getReadBreakpoints() const { return readBreakpoints; }

        bool isValidLine() const { return validLine; }

        const string &getSamLine() const { return samLine; }

        const vector<unsigned int> &getSamChunkPositions() const {
            return samTabPositions;
        }

        bool assessOutlierMateDistance();

        ChrIndex getMateChrIndex() const { return mateChrIndex; }

        ChrSize getMatePos() const { return matePos; }

        const vector<char> &getReadBreakpointTypes() const {
            return readBreakpointTypes;
        }

        void setChosenBp(ChrSize chosenBpLoc, int alignmentIndex);

        bool isOverhangEncounteredM() const { return chosenBp->bpEncounteredM; }

        ChrSize getOverhangLength() const { return ChrSize(chosenBp->overhangLength); }

        ChrSize getOverhangStartIndex() const { return ChrSize(chosenBp->overhangStartIndex); }

        vector<SuppAlignment> generateSuppAlignments(ChrIndex bpChrIndex, int bpPos);

        const vector<SuppAlignment> &getSupplementaryAlignments() const {
            return chosenBp->supplementaryAlignments;
        }

        ChrIndex getChrIndex() const { return chrIndex; }

        /** This returns a signed integer, because break-point sizes can be negative. **/
        const vector<signed int> &getReadBreakpointsSizes() const {
            return readBreakpointSizes;
        }

        /** true if mapq < 13 */
        bool isLowMapq() const { return lowMapq; }

        /** mapq 0 is treated as a special case, where number of SAs and
          * base qualities will be the sole determinants of read quality */
        bool isNullMapq() const { return nullMapq; }

        bool isSupplementary() const { return supplementary; }

        void addChildNode(int indexIn) { chosenBp->addChildNode(indexIn); }

        void addSupplementaryAlignments(const vector<SuppAlignment> &suppAlignments) {
            chosenBp->addSupplementaryAlignments(suppAlignments);
        }

        const vector<int> &getChildrenNodes() const {
            return chosenBp->childrenNodes;
        }

        int getOriginIndex() const { return chosenBp->selfNodeIndex; }

        string printOverhang() const;

        double overhangComplexityMaskRatio() const;

        bool isInvertedMate() const { return invertedMate; }

        bool isDistantMate() const { return distantMate == 1; }


      private:

        void mappingQualityCheck();

        /** The `Alignment::isEventCandidate` is true, if the last CIGAR code indicates a match,
         *  or if the CIGAR indicates a soft-clip, hard-clip, insertion, or deletion.
         */
        bool isEventCandidate() const;

        void createCigarChunks();

        void assignBreakpointsAndOverhangs();

        void qualityCheckCascade();

        bool clipCountCheck();

        bool uniqueSuppCheck();

        double overhangMedianQuality(const CigarChunk &cigarChunk) const;

        template <typename Iterator>
        void fullMedianQuality(Iterator qualBegin, Iterator qualEnd,
                               vector<int> &overhangPerBaseQuality) const;

        template <typename Iterator>
        double getMedian(Iterator begin, Iterator end) const;

        void assessReadType();

        bool lowMapq;

        bool nullMapq;

        int distantMate;

        unique_ptr<ChosenBp> chosenBp;

        ChrIndex chrIndex;

        int readType;

        ChrSize startPos, endPos;

        ChrIndex mateChrIndex;

        ChrSize matePos;

        string samLine;

        bool validLine;

        vector<unsigned int> samTabPositions;

        string::const_iterator saCbegin, saCend;

        bool hasSa;

        bool supplementary;

        bool fwdStrand;

        bool invertedMate;

        bool qualChecked;

        vector<CigarChunk> cigarChunks;

        vector<ChrSize> readBreakpoints;

        vector<char> readBreakpointTypes;

        vector<signed int> readBreakpointSizes;

        vector<double> readBreakpointComplexityMaskRatios;

        deque<bool> readBreakpointsEncounteredM;

        vector<OverhangRange> readOverhangCoords;

    };

} /* namespace sophia */

#endif /* ALIGNMENT_H_ */
