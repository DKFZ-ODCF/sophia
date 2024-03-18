#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <memory>

#include "GenericChrConverter.h"
#include "GlobalAppConfig.h"
#include "SuppAlignmentAnno.h"
#include "Fixtures.h"

namespace sophia {

    using namespace testing;

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor1) {
        const std::string test = "chr1:1041693|(8,0,!/0)";
        SuppAlignmentAnno anno = SuppAlignmentAnno(test);

        EXPECT_EQ(anno.getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr1"));

        EXPECT_EQ(anno.isEncounteredM(), false);
        EXPECT_EQ(anno.isInverted(), false);

        EXPECT_EQ(anno.getSupport(), 8);
        anno.setSupport(6);
        EXPECT_EQ(anno.getSupport(), 6);

        EXPECT_EQ(anno.getMateSupport(), 0);
        anno.incrementMateSupport();
        EXPECT_EQ(anno.getMateSupport(), 1);
        anno.setMateSupport(4);
        EXPECT_EQ(anno.getMateSupport(), 4);

        EXPECT_EQ(anno.getExtendedPos(), 1041693);
        EXPECT_EQ(anno.getPos(), 1041693);

        EXPECT_EQ(anno.getSecondarySupport(), 0);
        anno.setSecondarySupport(5);
        EXPECT_EQ(anno.getSecondarySupport(), 5);

        EXPECT_EQ(anno.isToRemove(), false);
        anno.setToRemove(true);
        EXPECT_EQ(anno.isToRemove(), true);

        EXPECT_EQ(anno.isDistant(), true);

        EXPECT_EQ(anno.isSuspicious(), true);
        anno.setSuspicious(false);
        EXPECT_EQ(anno.isSuspicious(), false);

        EXPECT_EQ(anno.isSemiSuspicious(), false);
        anno.setSemiSuspicious(true);
        EXPECT_EQ(anno.isSemiSuspicious(), true);

        EXPECT_EQ(anno.getExpectedDiscordants(), 0);
        anno.setExpectedDiscordants(10);
        EXPECT_EQ(anno.getExpectedDiscordants(), 10);

        EXPECT_EQ(anno.isFuzzy(), false);
        anno.setFuzzy(true);
        EXPECT_EQ(anno.isFuzzy(), true);

        EXPECT_EQ(anno.isStrictFuzzy(), false);

        EXPECT_EQ(anno.isStrictFuzzyCandidate(), false);

        EXPECT_EQ(anno.isProperPairErrorProne(), false);

        EXPECT_EQ(anno.getSupportingIndices().size(), 0);
        anno.addSupportingIndices(vector<int>({1, 2, 3}));
        EXPECT_THAT(anno.getSupportingIndices(), ElementsAre(1, 2, 3));

        // !suspicious, semisuspicious, 10 expected discordants, 6 supporting, 5 secondary supporting, 4 mate supporting
        EXPECT_EQ(anno.print(), "chr1:1041693|(6,5,4?/10)");
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor2) {
        const std::string test = "chr16:1041693_INV|(8,0,!/0)";

        SuppAlignmentAnno anno = SuppAlignmentAnno(test);

        EXPECT_EQ(anno.getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));

        EXPECT_EQ(anno.isEncounteredM(), false);
        EXPECT_EQ(anno.isInverted(), true);
        EXPECT_EQ(anno.isSuspicious(), true);
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor3) {
        const std::string test = "chr16:32084614-32084615_INV|(0,1,3?/6)";

        SuppAlignmentAnno anno = SuppAlignmentAnno(test);

        EXPECT_EQ(anno.getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));

        EXPECT_EQ(anno.isEncounteredM(), false);
        EXPECT_EQ(anno.isInverted(), true);
        EXPECT_EQ(anno.isSuspicious(), false);
        EXPECT_EQ(anno.isSemiSuspicious(), true);
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor4) {
        const std::string test = "|chr4:49107794-49107795(0,1,2?/89)";
        SuppAlignmentAnno anno = SuppAlignmentAnno(test);
        EXPECT_EQ(anno.getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr4"));

        EXPECT_EQ(anno.isEncounteredM(), true);
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor5) {
        const std::string test = "HLA-DRB1*13:01:01:2914|(4,0,0?/0)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("HLA-DRB1*13:01:01"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor6) {
        const std::string test = "|chr17:64795390_INV(8,0,6?/12)";

        SuppAlignmentAnno anno = SuppAlignmentAnno(test);

        EXPECT_TRUE(anno.isInverted());
        EXPECT_EQ(anno.isEncounteredM(), true);

        EXPECT_EQ(anno.getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr17"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor7) {
        const std::string test = "chr5:10019-11689|(11,0,13?/13)";

        SuppAlignmentAnno anno = SuppAlignmentAnno(test);

        EXPECT_EQ(anno.getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr5"));

        EXPECT_EQ(anno.isEncounteredM(), false);
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor8) {
        const std::string test = "chr5:10006-10487|(11,0,18?/18)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr5"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor9) {
        const std::string test = "chr5:10041-11708|(21,0,90?/90)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr5"));
    }

    TEST_F(GenericChrConverterFixture, SuppAlignmentAnnoTest_ParsingConstructor10) {
        const std::string test = "chr18:10007-10586|(8,0,9/9)";
        EXPECT_EQ(SuppAlignmentAnno(test).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr18"));
    }

} // namespace sophia