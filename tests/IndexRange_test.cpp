#include <gtest/gtest.h>
#include <memory>

#include "IndexRange.h"

namespace sophia {

    TEST(IndexRangeTest, IndexRangeTest_constructor) {
        IndexRange range(1, 10);
        EXPECT_EQ(range.start(), 1);
        EXPECT_EQ(range.end(), 10);

        EXPECT_THROW(IndexRange(1, 0), std::invalid_argument);
    }

    TEST(IndexRangeTest, IndexRangeTest_contains) {
        IndexRange range(1, 11);
        EXPECT_FALSE(range.contains(0));
        EXPECT_TRUE(range.contains(1));
        EXPECT_TRUE(range.contains(5));
        EXPECT_TRUE(range.contains(10));
        EXPECT_FALSE(range.contains(11));
    }

    TEST(IndexRangeTest, IndexRangeTest_width) {
        IndexRange range(1, 11);
        EXPECT_EQ(range.width(), 10);
    }

}  // namespace sophia