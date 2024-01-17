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

        // Breakpoint contains `:`
        const std::string test2 = "chr22\t37934897\t37934898\t0,0,4,4,0,0,32,2,0,0,0,0\t32,36\t.\tHLA-DRB1*13:01:01:2914|(4,0,0?/0);chr1:28514292-28514293_INV|(4,0,1?/5);HLA-DRB1*13:01:01:2919|(4,0,0?/0);chr13:20527706-20527707_INV|(4,0,1?/5);HLA-DRB1*13:01:01:2922|(4,0,0?/0);chr1:120948617-120948618_INV|(4,0,1?/5);HLA-DRB1*13:01:01:2923|(4,0,0?/0);chr13:92964918-92964919_INV|(4,0,1?/5)\t>28974_1:ATGTCCACGGTAAAAAATTTGAATTTTATTT|(4)";
        EXPECT_EQ(Breakpoint::parse(test2, true).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr22"));

    }

} // namespace sophia