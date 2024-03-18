#include "IndexRange.h"
#include "global.h"

namespace sophia {

    IndexRange::IndexRange(ChrIndex start, ChrIndex end)
            : start_(start), end_(end) {
        if (start > end) {
            throw_with_trace(std::invalid_argument("IndexRange: start must be <= end"));
        }
    }

    ChrIndex IndexRange::start() const {
        return start_;
    }
    ChrIndex IndexRange::end() const {
        return end_;
    }

    ChrSize IndexRange::width() const {
        return static_cast<ChrSize>(end_ - start_);
    }

    bool IndexRange::contains(const ChrIndex &index) const {
        // 0-based, left-inclusive, right-exclusive range.
        return index >= start_ && index < end_;
    }

} // namespace sophia
