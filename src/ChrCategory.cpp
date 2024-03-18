#include "global.h"
#include "ChrCategory.h"
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/unordered/unordered_map.hpp>

namespace sophia {

    const boost::unordered::unordered_map<std::string, const ChrCategory> ChrCategory::categories = {
            {"AUTOSOME", ChrCategory("AUTOSOME", 0)},
            {"X", ChrCategory("X", 1)},
            {"Y", ChrCategory("Y", 2)},
            {"EXTRACHROMOSOMAL", ChrCategory("EXTRACHROMOSOMAL", 3)},
            {"UNASSIGNED", ChrCategory("UNASSIGNED", 4)},
            {"ALT", ChrCategory("ALT", 5)},
            {"HLA", ChrCategory("HLA", 6)},
            {"VIRUS", ChrCategory("VIRUS", 7)},
            {"DECOY", ChrCategory("DECOY", 8)},
            {"TECHNICAL", ChrCategory("TECHNICAL", 9)}
        };

    const ChrCategory& ChrCategory::AUTOSOME = ChrCategory::categories.at("AUTOSOME");
    const ChrCategory& ChrCategory::X = ChrCategory::categories.at("X");
    const ChrCategory& ChrCategory::Y = ChrCategory::categories.at("Y");
    const ChrCategory& ChrCategory::EXTRACHROMOSOMAL = ChrCategory::categories.at("EXTRACHROMOSOMAL");
    const ChrCategory& ChrCategory::UNASSIGNED = ChrCategory::categories.at("UNASSIGNED");
    const ChrCategory& ChrCategory::ALT = ChrCategory::categories.at("ALT");
    const ChrCategory& ChrCategory::HLA = ChrCategory::categories.at("HLA");
    const ChrCategory& ChrCategory::VIRUS = ChrCategory::categories.at("VIRUS");
    const ChrCategory& ChrCategory::DECOY = ChrCategory::categories.at("DECOY");
    const ChrCategory& ChrCategory::TECHNICAL = ChrCategory::categories.at("TECHNICAL");


    const std::vector<ChrCategory> ChrCategory::sorted_categories = {
        ChrCategory::AUTOSOME,
        ChrCategory::X,
        ChrCategory::Y,
        ChrCategory::EXTRACHROMOSOMAL,
        ChrCategory::UNASSIGNED,
        ChrCategory::ALT,
        ChrCategory::HLA,
        ChrCategory::VIRUS,
        ChrCategory::DECOY,
        ChrCategory::TECHNICAL
    };

    ChrCategory::ChrCategory(const std::string &in, size_type index)
        : category_name { in },
          category_index { index } {}

    ChrCategory::~ChrCategory() {}

    const ChrCategory& ChrCategory::from_string(const std::string &in) {
        std::string normalizedIn = boost::algorithm::to_upper_copy(in);
        if (categories.find(normalizedIn) == categories.end()) {
            throw_with_trace(std::invalid_argument("Unknown chromosome category: '" + in + "'"));
        }
        return categories.at(normalizedIn);
    }

    std::string ChrCategory::getName() const {
        return category_name;
    }

    ChrCategory::size_type ChrCategory::numCategories() {
        return categories.size();
    }

    const std::vector<ChrCategory>& ChrCategory::getCategories() {
        return sorted_categories;
    }

    bool ChrCategory::operator==(const ChrCategory &other) const {
        return category_index == other.category_index;
    }

    bool ChrCategory::operator!=(const ChrCategory &other) const {
        return category_index != other.category_index;
    }

    bool ChrCategory::operator<(const ChrCategory &other) const {
        return category_index < other.category_index;
    }

    bool ChrCategory::operator>(const ChrCategory &other) const {
        return category_index > other.category_index;
    }

} // namespace sophia
