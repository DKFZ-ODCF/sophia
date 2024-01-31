#include <gtest/gtest.h>
#include <string>
#include <stdexcept>
#include "Fixtures.h"
#include "GlobalAppConfig.h"
#include "GenericChrConverter.h"

namespace sophia {

    TEST_F(GenericChrConverterFixture, GenericChrConverterTest_chrNameToIndex) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.chrNameToIndex("chr1"), 0);
        EXPECT_EQ(converter.chrNameToIndex("chr22"), 21);
        EXPECT_EQ(converter.chrNameToIndex("chrX"), 22);
        EXPECT_EQ(converter.chrNameToIndex("chrY"), 23);
        EXPECT_EQ(converter.chrNameToIndex("chrM"), 24);
        EXPECT_EQ(converter.chrNameToIndex("chrUn_KI270423v1"), 88);
        EXPECT_EQ(converter.chrNameToIndex("phix"), 3366);
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverter_assemblyName) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.assemblyName, "hg38");
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverter_nChromosomes) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.nChromosomes(), 3367);
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverterTest_indexToChrName) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.indexToChrName(0), "chr1");
        EXPECT_EQ(converter.indexToChrName(21), "chr22");
        EXPECT_EQ(converter.indexToChrName(22), "chrX");
        EXPECT_EQ(converter.indexToChrName(23), "chrY");
        EXPECT_EQ(converter.indexToChrName(24), "chrM");
        EXPECT_EQ(converter.indexToChrName(88), "chrUn_KI270423v1");
        EXPECT_EQ(converter.indexToChrName(3366), "phix");
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverter_chrNameToIndex) {
        const GenericChrConverter &converter =
            dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.chrNameToIndex("chr1"), 0);
        EXPECT_EQ(converter.chrNameToIndex("chr22"), 21);
        EXPECT_EQ(converter.chrNameToIndex("chrX"), 22);
        EXPECT_EQ(converter.chrNameToIndex("chrY"), 23);
        EXPECT_EQ(converter.chrNameToIndex("chrM"), 24);
        EXPECT_EQ(converter.chrNameToIndex("chrUn_KI270423v1"), 88);
        EXPECT_EQ(converter.chrNameToIndex("phix"), 3366);
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverterTest_nChromosomesCompressedMref) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.nChromosomesCompressedMref(), 3365);
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverterTest_is_category) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
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

    TEST_F(GenericChrConverterFixture, GenericChrConverter_isCompressedMrefIndex) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_TRUE(converter.isCompressedMref(0));  // chr1
        EXPECT_TRUE(converter.isCompressedMref(21)); // chr22
        EXPECT_TRUE(converter.isCompressedMref(22)); // chrX
        EXPECT_TRUE(converter.isCompressedMref(23)); // chrY
        EXPECT_TRUE(! converter.isCompressedMref(24)); // chrM
        EXPECT_TRUE(converter.isCompressedMref(88)); // chrUn_KI270423v1
        EXPECT_TRUE(! converter.isCompressedMref(3366)); // phix
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverter_compressedMrefIndexToIndex) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.compressedMrefIndexToIndex(0), 0);  // chr1
        EXPECT_EQ(converter.compressedMrefIndexToIndex(23), 23); // chrY
        // chrM = global::24 is missing from compressed mrefs,
        // therefore, compressed::24 = global::25
        EXPECT_EQ(converter.compressedMrefIndexToIndex(24), 25); // chr1_KI270706v1_random
        EXPECT_THROW(converter.compressedMrefIndexToIndex(3366), std::logic_error); // phix
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverter_indexToCompressedMrefIndex) {
        const GenericChrConverter &converter =
                dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.indexToCompressedMrefIndex(0), 0);   // chr1
        EXPECT_EQ(converter.indexToCompressedMrefIndex(23), 23); // chrY
        // chrM = global::24 is missing from compressed mrefs, ...
        // [this is an interesting case, because it demonstrates that the GenericChrConverter
        //  is not constrained to support only compressed Mref chromosomes strictly separated from
        //  the rest of the chromosomes, but that they isCompressedMref() flag can freely be used.]
        EXPECT_THROW(converter.indexToCompressedMrefIndex(24), std::logic_error); // phix
        // ... therefore, compressed::24 = global::25
        EXPECT_EQ(converter.indexToCompressedMrefIndex(25), 24); // chr1_KI270706v1_random
        EXPECT_THROW(converter.indexToCompressedMrefIndex(3366), std::logic_error); // phix
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverterTest_ParseSimpleStrings) {
        const std::string test1 = "chr1\tsomething\telse\n";
        const GenericChrConverter &converter =
            dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.parseChr(test1.begin(), test1.end(), '\t'),
                  "chr1");
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverter_chrSizeCompressedMref) {
        const GenericChrConverter &converter =
            dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());
        EXPECT_EQ(converter.chrSizeCompressedMref(0), 248956422);  // chr1
        EXPECT_EQ(converter.chrSizeCompressedMref(23), 57227415);  // chrY
        EXPECT_EQ(converter.chrSizeCompressedMref(454), 171823);  // chrEBV; 454 is the index in compressed mrefs!
        // chrM is not in compressed mrefs
        EXPECT_EQ(converter.chrSizeCompressedMref(24), 175055); // chr1_KI270706v1_random
    }

    TEST_F(GenericChrConverterFixture, GenericChrConverterTest_ParseBreakPointStrings) {
        const std::string stopChars = "|(,!/?;";
        const GenericChrConverter &converter =
            dynamic_cast<const GenericChrConverter&>(GlobalAppConfig::getInstance().getChrConverter());


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