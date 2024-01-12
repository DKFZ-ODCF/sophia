#include <gtest/gtest.h>
#include "Hg38ChrConverter.h"

namespace sophia {

    Hg38ChrConverter converter = Hg38ChrConverter();

    TEST(Hg38ChrConverterTest_ParseSimpleStrings, BasicAssertions) {
        const std::string test1 = "chr1\tsomething\telse\n";
        EXPECT_EQ(converter.parseChr(test1.begin(), test1.end(), '\t'),
                  "chr1");
    }

    TEST(Hg38ChrConverterTest_ParseBreakPointStrings, BasicAssertions) {
        const std::string stopChars = "|(,!/?";

        const std::string test1 = "HLA-DRB1*13:01:01:2914|(4,0,0?/0)";
        EXPECT_EQ(converter.parseChr(test1.begin(), test1.end(), ':', stopChars),
                  "HLA-DRB1*13:01:01");

        const std::string test2 = "chrUn_KI270749v1:13653-13654(1,0,3?/4)";
        EXPECT_EQ(converter.parseChr(test2.begin(), test2.end(), ':', stopChars),
                  "chrUn_KI270749v1");

        const std::string test3 = test1 + ";" + test2;
        EXPECT_EQ(converter.parseChr(test3.begin(), test3.end(), ':', stopChars),
                  "HLA-DRB1*13:01:01");
    }

} // namespace sophia