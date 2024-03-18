#include <gtest/gtest.h>
#include <string>
#include <stdexcept>
#include <vector>
#include "global.h"
#include "Hg37ChrConverter.h"

namespace sophia {

    TEST(Hg37ChrConverterTest, Hg37ChrConverterTest_chrNameToIndex) {
        std::vector<CompressedMrefIndex> mrefIndex {1003, 0, 1003, 1};
        std::vector<ChrIndex> result {1, 3};
        EXPECT_EQ(Hg37ChrConverter::_buildCompressedMrefIndexToIndex(2, mrefIndex),
                  result);
    }
} // namespace sophia