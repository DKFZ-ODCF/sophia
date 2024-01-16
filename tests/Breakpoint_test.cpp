#include <gtest/gtest.h>
#include <memory>

#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"
#include "Breakpoint.h"

namespace sophia {

    TEST(BreakpointTest_Parse, BasicAssertions) {
        try {
            // TODO Find a way to deal with the singleton, without having to refactor the whole code
            GlobalAppConfig::init(move(make_unique<Hg38ChrConverter>()));
        } catch (const std::exception& e) {
        }

        const std::string test1 = "chr22\t10525762\t10525763\t0,0,7,0,0,0,4,0,0,9,0,0\t4,4\t|chrUn_KI270749v1:73502-73503(0,2,4?/8);|chr22_KI270735v1_random:17508-17509_INV(0,1,3?/6)\t.\t.";
        EXPECT_EQ(Breakpoint::parse(test1, true).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr22"));
//
//        const std::string test2 = "";
//        EXPECT_EQ(Breakpoint::parse(test2, true).getChrIndex(),
//                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr16"));
//
//        const std::string test3 = "";
//        EXPECT_EQ(Breakpoint::parse(test3, true).getChrIndex(),
//                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr4"));
    }

} // namespace sophia