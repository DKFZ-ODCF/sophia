#ifndef INDEXRANGE_H_
#define INDEXRANGE_H_

#include "global.h"

namespace sophia {

    /** A simple class to represent a 0-based range of indices.
      * Chromosome ranges need to be 0-based and right-exclusive. For instance, in the range
      * [0, 10), the first contained index is 0 and the last contained index is 9.
      *
      * The condition start <= end has to hold.
      *
      * If you set start == end, then you obtain a range of `width()` == 0.
      * you can use this to define a range that does not contain any values, i.e. where
      * `contains()` always returns false.
      **/
    class IndexRange {

    private:

        ChrIndex start_;
        ChrIndex end_;

    public:

        /* Throws invalid_argument exception, if start > end */
         IndexRange(ChrIndex start, ChrIndex end);
        ~IndexRange() {}

        ChrIndex start() const;
        ChrIndex end() const;

        ChrSize width() const;

        bool contains(const ChrIndex &index) const;
    };

}

#endif // INDEXRANGE_H_