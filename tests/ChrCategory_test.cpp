#include <gtest/gtest.h>
#include <string>
#include <stdexcept>
#include "ChrCategory.h"

namespace sophia {

    TEST(ChrCategoryTest, ChrCategoryTest_from_string) {
        ASSERT_EQ(ChrCategory::from_string("AUTOSOME"), ChrCategory::AUTOSOME);
        ASSERT_EQ(ChrCategory::from_string("gonosome"), ChrCategory::GONOSOME);
        ASSERT_EQ(ChrCategory::from_string("GONOSOME"), ChrCategory::GONOSOME);
        ASSERT_EQ(ChrCategory::from_string("extraChromosomal"), ChrCategory::EXTRACHROMOSOMAL);
        ASSERT_EQ(ChrCategory::from_string("uNasSIGNED"), ChrCategory::UNASSIGNED);
        ASSERT_EQ(ChrCategory::from_string("alt"), ChrCategory::ALT);
        ASSERT_EQ(ChrCategory::from_string("hLa"), ChrCategory::HLA);
        ASSERT_EQ(ChrCategory::from_string("decoy"), ChrCategory::DECOY);
        ASSERT_EQ(ChrCategory::from_string("virus"), ChrCategory::VIRUS);
        ASSERT_EQ(ChrCategory::from_string("technical"), ChrCategory::TECHNICAL);
    }

    TEST(ChrCategoryTest, ChrCategoryTest_from_string_invalid) {
        ASSERT_THROW(ChrCategory::from_string("invalid"), std::invalid_argument);
    }

    TEST(ChrCategoryTest, ChrCategoryTest_getName) {
        ASSERT_EQ(ChrCategory::AUTOSOME.getName(), "AUTOSOME");
        ASSERT_EQ(ChrCategory::GONOSOME.getName(), "GONOSOME");
        ASSERT_EQ(ChrCategory::EXTRACHROMOSOMAL.getName(), "EXTRACHROMOSOMAL");
        ASSERT_EQ(ChrCategory::UNASSIGNED.getName(), "UNASSIGNED");
        ASSERT_EQ(ChrCategory::ALT.getName(), "ALT");
        ASSERT_EQ(ChrCategory::HLA.getName(), "HLA");
        ASSERT_EQ(ChrCategory::DECOY.getName(), "DECOY");
        ASSERT_EQ(ChrCategory::VIRUS.getName(), "VIRUS");
        ASSERT_EQ(ChrCategory::TECHNICAL.getName(), "TECHNICAL");
    }

    TEST(ChrCategoryTest, ChrCategoryTest_numCategories) {
        ASSERT_EQ(ChrCategory::numCategories(), 9);
    }

    TEST(ChrCategoryTest, ChrCategoryTest_getCategories) {
        std::vector<ChrCategory> categories = ChrCategory::getCategories();
        ASSERT_EQ(categories.size(), 9);
        ASSERT_EQ(categories[0], ChrCategory::AUTOSOME);
        ASSERT_EQ(categories[1], ChrCategory::GONOSOME);
        ASSERT_EQ(categories[2], ChrCategory::EXTRACHROMOSOMAL);
        ASSERT_EQ(categories[3], ChrCategory::UNASSIGNED);
        ASSERT_EQ(categories[4], ChrCategory::ALT);
        ASSERT_EQ(categories[5], ChrCategory::HLA);
        ASSERT_EQ(categories[6], ChrCategory::VIRUS);
        ASSERT_EQ(categories[7], ChrCategory::DECOY);
        ASSERT_EQ(categories[8], ChrCategory::TECHNICAL);
    }

} // namespace sophia