#include <gtest/gtest.h>
#include <memory>

#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"
#include "SuppAlignment.h"

namespace sophia {

    TEST(SuppAlignmentTest_Parse, BasicAssertions) {
        GlobalAppConfig::init(move(make_unique<Hg38ChrConverter>()));

        const std::string test1 = "chr2\t10525827\t10525828\t0,0,3,0,0,0,12,1,3,7,0,0\t12,12\tchr16:32084614-32084615_INV|(0,1,3?/6)\t.\t.";
        EXPECT_EQ(SuppAlignment::parse(test1).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr2"));

        const std::string test2 = "chr22\t10722811\t10722812\t0,0,25,0,0,0,306,32,2,10,1,0\t304,306\t|chr4:49107794-49107795(0,1,2?/89);|chr17_KI270729v1_random:22495-22496_INV(0,1,2?/89);|chr20:31075717-31075718_INV(0,1,3?/91);|chr17_KI270729v1_random:997-998(0,1,2?/89);|chr5:49601606-49601607_INV(0,1,2?/89);|chr20:31069449-31069450_INV(0,1,2?/89);|chr17_KI270729v1_random:21047-21048_INV(0,1,2?/89);|chr20:62657-62658(0,1,3?/91);|chr17_KI270729v1_random:24683-24684_INV(0,1,2?/89);|chr5:49666980-49666981(0,1,3?/91)\t|chr22:10723582-10723872(0,0,37/85)\t.";
        EXPECT_EQ(SuppAlignment::parse(test2).getChrIndex(),
                  GlobalAppConfig::getInstance().getChrConverter().chrNameToIndex("chr22"));
    }

} // namespace sophia