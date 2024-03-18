#include "MateInfo.h"

namespace sophia {


        bool MateInfo::operator<(const MateInfo &rhs) const {
            if (mateChrIndex < rhs.mateChrIndex)
                return true;
            if (mateChrIndex > rhs.mateChrIndex)
                return false;
            if (mateStartPos < rhs.mateStartPos)
                return true;
            return false;
        }

        bool MateInfo::suppAlignmentFuzzyMatch(const SuppAlignment &sa) const {
            if (mateChrIndex != sa.getChrIndex()) {
                return false;
            } else {
                if (!sa.isFuzzy()) {
                    return static_cast<long>(sa.getPos()) >= (static_cast<long>(mateStartPos) - static_cast<long>(sa.getMatchFuzziness())) &&
                           static_cast<long>(sa.getPos()) <= (static_cast<long>(mateEndPos) + static_cast<long>(sa.getMatchFuzziness()));
                } else {
                    return (static_cast<long>(mateStartPos) - static_cast<long>(sa.getMatchFuzziness())) <= static_cast<long>(sa.getExtendedPos()) &&
                           static_cast<long>(sa.getPos()) <= (static_cast<long>(mateEndPos) + static_cast<long>(sa.getMatchFuzziness()));
                }
            }
        }

        MateInfo::MateInfo(ChrSize readStartPosIn,
                           ChrSize readEndPosIn,
                           ChrIndex mateChrIndexIn,
                           ChrSize mateStartPosIn,
                           int sourceType,
                           bool invertedIn)
            : readStartPos{readStartPosIn},
              readEndPos{readEndPosIn},
              mateChrIndex{mateChrIndexIn},
              mateStartPos{mateStartPosIn},
              mateEndPos{mateStartPosIn},
              inverted{invertedIn},
              source{sourceType},
              evidenceLevel{sourceType == 2 ? 3 : 1},
              matePower{1},
              inversionSupport{invertedIn},
              straightSupport{!invertedIn},
              bpLocs{},
              saSupporter{false},
              toRemove{false} {}

        MateInfo::MateInfo(ChrSize readStartPosIn,
                           ChrSize readEndPosIn,
                           ChrIndex mateChrIndexIn,
                           ChrSize mateStartPosIn,
                           int sourceType,
                           bool invertedIn,
                           const std::vector<ChrSize> &bpLocsIn)
            : readStartPos{readStartPosIn},
              readEndPos{readEndPosIn},
              mateChrIndex{mateChrIndexIn},
              mateStartPos{mateStartPosIn},
              mateEndPos{mateStartPosIn},
              inverted{invertedIn},
              source{sourceType},
              evidenceLevel{sourceType == 2 ? 3 : 1},
              matePower{1},
              inversionSupport{invertedIn},
              straightSupport{!invertedIn},
              bpLocs{bpLocsIn},
              saSupporter{false},
              toRemove{false} {}

        bool MateInfo::isToRemove() const { return toRemove; }


} /* namespace sophia */