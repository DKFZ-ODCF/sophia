#include "global.h"
#include "ChrCategory.h"
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/unordered/unordered_set.hpp>

namespace sophia {

    const boost::unordered::unordered_map<std::string, const ChrCategory> ChrCategory::categories = {
            {"AUTOSOME", ChrCategory("AUTOSOME", 0)},
            {"GONOSOME", ChrCategory("GONOSOME", 1)},
            {"EXTRACHROMOSOMAL", ChrCategory("EXTRACHROMOSOMAL", 2)},
            {"UNASSIGNED", ChrCategory("UNASSIGNED", 3)},
            {"ALT", ChrCategory("ALT", 4)},
            {"HLA", ChrCategory("HLA", 5)},
            {"VIRUS", ChrCategory("VIRUS", 6)},
            {"DECOY", ChrCategory("DECOY", 7)},
            {"TECHNICAL", ChrCategory("TECHNICAL", 8)}
        };

    const ChrCategory& ChrCategory::AUTOSOME = ChrCategory::categories.at("AUTOSOME");
    const ChrCategory& ChrCategory::GONOSOME = ChrCategory::categories.at("GONOSOME");
    const ChrCategory& ChrCategory::EXTRACHROMOSOMAL = ChrCategory::categories.at("EXTRACHROMOSOMAL");
    const ChrCategory& ChrCategory::UNASSIGNED = ChrCategory::categories.at("UNASSIGNED");
    const ChrCategory& ChrCategory::ALT = ChrCategory::categories.at("ALT");
    const ChrCategory& ChrCategory::HLA = ChrCategory::categories.at("HLA");
    const ChrCategory& ChrCategory::VIRUS = ChrCategory::categories.at("VIRUS");
    const ChrCategory& ChrCategory::DECOY = ChrCategory::categories.at("DECOY");
    const ChrCategory& ChrCategory::TECHNICAL = ChrCategory::categories.at("TECHNICAL");

    ChrCategory::ChrCategory(const std::string &in, size_type index)
        : category_name { in },
          category_index { index } {}

    ChrCategory::~ChrCategory() {}

    const ChrCategory& ChrCategory::from_string(const std::string &in) {
        std::string normalizedIn = boost::algorithm::to_upper_copy(in);
        if (categories.find(normalizedIn) == categories.end()) {
            throw_with_trace(std::runtime_error("Unknown chromosome category: '" + in + "'"));
        }
        return categories.at(normalizedIn);
    }

    std::string ChrCategory::getName() const {
        return category_name;
    }

    ChrCategory::size_type ChrCategory::numCategories() {
        return categories.size();
    }

    std::vector<ChrCategory> ChrCategory::getCategories() {
        std::vector<ChrCategory> result;
        result.reserve(numCategories());
        std::transform(
            categories.begin(),
            categories.end(),
            std::back_inserter(result),
            [](const std::pair<std::string, ChrCategory>& kv) { return kv.second; });
        return result;
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
