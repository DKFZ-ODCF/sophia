#include <gtest/gtest.h>
#include <string>
#include <stdexcept>
#include "Fixtures.h"
#include "GlobalAppConfig.h"
#include "Hg38ChrConverter.h"

namespace sophia {

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverterTest_chrNameToIndex) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.chrNameToIndex("chr1"), 0);
        EXPECT_EQ(converter.chrNameToIndex("chr22"), 21);
        EXPECT_EQ(converter.chrNameToIndex("chrX"), 22);
        EXPECT_EQ(converter.chrNameToIndex("chrY"), 23);
        EXPECT_EQ(converter.chrNameToIndex("chrM"), 24);
        EXPECT_EQ(converter.chrNameToIndex("chrUn_KI270423v1"), 88);
        EXPECT_EQ(converter.chrNameToIndex("phix"), 3366);
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverter_assemblyName) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.assemblyName, "hg38");
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverter_nChromosomes) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.nChromosomes(), 3367);
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverterTest_indexToChrName) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.indexToChrName(0), "chr1");
        EXPECT_EQ(converter.indexToChrName(21), "chr22");
        EXPECT_EQ(converter.indexToChrName(22), "chrX");
        EXPECT_EQ(converter.indexToChrName(23), "chrY");
        EXPECT_EQ(converter.indexToChrName(24), "chrM");
        EXPECT_EQ(converter.indexToChrName(88), "chrUn_KI270423v1");
        EXPECT_EQ(converter.indexToChrName(3366), "phix");
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverter_chrNameToIndex) {
        const Hg38ChrConverter &converter =
            dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.chrNameToIndex("chr1"), 0);
        EXPECT_EQ(converter.chrNameToIndex("chr22"), 21);
        EXPECT_EQ(converter.chrNameToIndex("chrX"), 22);
        EXPECT_EQ(converter.chrNameToIndex("chrY"), 23);
        EXPECT_EQ(converter.chrNameToIndex("chrM"), 24);
        EXPECT_EQ(converter.chrNameToIndex("chrUn_KI270423v1"), 88);
        EXPECT_EQ(converter.chrNameToIndex("phix"), 3366);
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverterTest_nChromosomesCompressedMref) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.nChromosomesCompressedMref(), 3365);
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverterTest_is_category) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_TRUE(converter.isAutosome(0));  // chr1
        EXPECT_TRUE(converter.isAutosome(21)); // chr22
        EXPECT_TRUE(converter.isGonosome(22)); // chrX
        EXPECT_TRUE(converter.isGonosome(23)); // chrY
        EXPECT_TRUE(converter.isExtrachromosomal(24)); // chrM
        EXPECT_TRUE(converter.isUnassigned(88)); // chrUn_KI270423v1
        EXPECT_TRUE(converter.isTechnical(3366)); // phix
        EXPECT_TRUE(converter.isVirus(455)); // chrEBV
        EXPECT_TRUE(converter.isDecoy(456)); // chrUn_KN707606v1_decoy
        EXPECT_TRUE(converter.isALT(194)); // chr1_KI270762v1_alt
        EXPECT_TRUE(converter.isHLA(2841)); // HLA-A*01:01:01:01
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverter_isCompressedMrefIndex) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_TRUE(converter.isCompressedMrefIndex(0));  // chr1
        EXPECT_TRUE(converter.isCompressedMrefIndex(21)); // chr22
        EXPECT_TRUE(converter.isCompressedMrefIndex(22)); // chrX
        EXPECT_TRUE(converter.isCompressedMrefIndex(23)); // chrY
        EXPECT_TRUE(! converter.isCompressedMrefIndex(24)); // chrM
        EXPECT_TRUE(converter.isCompressedMrefIndex(88)); // chrUn_KI270423v1
        EXPECT_TRUE(! converter.isCompressedMrefIndex(3366)); // phix
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverter_compressedMrefIndexToIndex) {
        const Hg38ChrConverter &converter =
                dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.compressedMrefIndexToIndex(0), 0);  // chr1
        EXPECT_EQ(converter.compressedMrefIndexToIndex(23), 23); // chrY
        // chrM = global::24 is missing from compressed mrefs,
        // therefore, compressed::24 = global::25
        EXPECT_EQ(converter.compressedMrefIndexToIndex(24), 25); // chr1_KI270706v1_random
        EXPECT_THROW(converter.compressedMrefIndexToIndex(3366), std::logic_error); // phix
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverterTest_ParseSimpleStrings) {
        const std::string test1 = "chr1\tsomething\telse\n";
        const Hg38ChrConverter &converter =
            dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.parseChr(test1.begin(), test1.end(), '\t'),
                  "chr1");
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverter_chrSizeCompressedMref) {
        const Hg38ChrConverter &converter =
            dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.chrSizeCompressedMref(0), 248956422);  // chr1
        EXPECT_EQ(converter.chrSizeCompressedMref(23), 57227415);  // chrY
        EXPECT_EQ(converter.chrSizeCompressedMref(454), 171823);  // chrEBV; 454 is the index in compressed mrefs!
        // chrM is not in compressed mrefs
        EXPECT_EQ(converter.chrSizeCompressedMref(24), 175055); // chr1_KI270706v1_random
    }

    TEST_F(Hg38ChrConverterFixture, Hg38ChrConverterTest_ParseBreakPointStrings) {
        const std::string stopChars = "|(,!/?;";
        const Hg38ChrConverter &converter =
            dynamic_cast<const Hg38ChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());


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