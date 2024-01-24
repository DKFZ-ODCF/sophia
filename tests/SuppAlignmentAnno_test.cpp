#include <gtest/gtest.h>
#include <memory>

#include "GenericChrConverter.h"
#include "GlobalAppConfig.h"
#include "SuppAlignmentAnno.h"
#include "Fixtures.h"

namespace sophia {

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor1) {
        const std::string test = "chr1:1041693|(8,0,!/0)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr1"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor2) {
        const std::string test = "chr16:1041693_INV|(8,0,!/0)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor3) {
        const std::string test = "chr16:32084614-32084615_INV|(0,1,3?/6)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor4) {
        const std::string test = "|chr4:49107794-49107795(0,1,2?/89)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr4"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor5) {
        const std::string test = "HLA-DRB1*13:01:01:2914|(4,0,0?/0)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("HLA-DRB1*13:01:01"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor6) {
        const std::string test = "|chr17:64795390_INV(8,0,6?/12)";
        EXPECT_TRUE(SuppAlignmentAnno(test).isInverted());
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr17"));
    }

} // namespace sophia