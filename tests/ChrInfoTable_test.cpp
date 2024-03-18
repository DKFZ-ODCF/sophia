#include <gtest/gtest.h>
#include <string>
#include <stdexcept>
#include "Fixtures.h"
#include "GlobalAppConfig.h"
#include "GenericChrConverter.h"

namespace sophia {

    TEST_F(ChrInfoTableFixture, ChrInfoTableTest_nChromosomes) {
        EXPECT_EQ(chr_info_table->nChromosomes(), 3367);
    }

    TEST_F(ChrInfoTableFixture, ChrInfoTableTest_getNames) {
        EXPECT_EQ(chr_info_table->getNames()[0], "chr1");
        EXPECT_EQ(chr_info_table->getNames()[21], "chr22");
        EXPECT_EQ(chr_info_table->getNames()[22], "chrX");
        EXPECT_EQ(chr_info_table->getNames()[23], "chrY");
        EXPECT_EQ(chr_info_table->getNames()[24], "chrM");
        EXPECT_EQ(chr_info_table->getNames()[3366], "phix");
        EXPECT_EQ(chr_info_table->getNames()[88], "chrUn_KI270423v1");
    }

    TEST_F(ChrInfoTableFixture, ChrInfoTableTest_getNamesByCategory) {
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::AUTOSOME).size(), 22);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::AUTOSOME)[0], "chr1");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::X).size(), 1);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::X)[0], "chrX");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::Y).size(), 1);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::Y)[0], "chrY");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::EXTRACHROMOSOMAL).size(), 1);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::EXTRACHROMOSOMAL)[0], "chrM");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::UNASSIGNED).size(), 169);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::UNASSIGNED)[0], "chr1_KI270706v1_random");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::TECHNICAL).size(), 1);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::TECHNICAL)[0], "phix");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::VIRUS).size(), 1);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::VIRUS)[0], "chrEBV");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::DECOY).size(), 2385);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::DECOY)[0], "chrUn_KN707606v1_decoy");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::ALT).size(), 261);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::ALT)[0], "chr1_KI270762v1_alt");

        EXPECT_EQ(chr_info_table->getNames(ChrCategory::HLA).size(), 525);
        EXPECT_EQ(chr_info_table->getNames(ChrCategory::HLA)[0], "HLA-A*01:01:01:01");
    }

} // namespace sophia