#include <gtest/gtest.h>
#include <memory>

#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"
#include "Breakpoint.h"

namespace sophia {

    TEST(BreakpointTest_Parse, BasicAssertions) {
        GlobalAppConfig::init(move(make_unique<Hg38ChrConverter>()));

        const std::string test1 = "";
        EXPECT_EQ(Breakpoint::parse(test1, true).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));

        const std::string test2 = "";
        EXPECT_EQ(Breakpoint::parse(test2, true).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));

        const std::string test3 = "";
        EXPECT_EQ(Breakpoint::parse(test3, true).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr4"));
    }

} // namespace sophia