#include <gtest/gtest.h>
#include <memory>

#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"
#include "SuppAlignment.h"

namespace sophia {

    TEST(SuppAlignmentTest_ParseSaSupport, BasicAssertions) {
        GlobalAppConfig::init(move(make_unique<Hg38ChrConverter>()));

        const std::string test1 = "chr16:1041693_INV|(8,0,!/0)";
        EXPECT_EQ(SuppAlignment::parseSaSupport(test1).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));

        const std::string test2 = "chr16:32084614-32084615_INV|(0,1,3?/6)";
        EXPECT_EQ(SuppAlignment::parseSaSupport(test2).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));

        const std::string test3 = "|chr4:49107794-49107795(0,1,2?/89)";
        EXPECT_EQ(SuppAlignment::parseSaSupport(test3).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr4"));
    }

} // namespace sophia