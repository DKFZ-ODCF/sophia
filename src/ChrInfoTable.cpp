#include "ChrInfoTable.h"
#include <algorithm>
#include <stdexcept>
#include <boost/unordered/unordered_map.hpp>


namespace sophia {

    /** Just a helper for the constructor. */
    boost::unordered::unordered_map<ChrCategory, std::vector<ChrInfo>>
    ChrInfoTable::buildChrInfosByCategory(const std::vector<ChrInfo> &chr_info) {
        // There must be a vector of `ChrInfo` for each existing category.
        boost::unordered::unordered_map<ChrCategory, std::vector<ChrInfo>> result;
        result.reserve(ChrCategory::numCategories());
        for (auto category : ChrCategory::getCategories()) {
            result.try_emplace(category, std::vector<ChrInfo>());
        }
        // Elements are added one by one to their corresponding category vector.
        for (const auto &info : chr_info) {
            result.at(info.getCategory()).push_back(info);
        }
        return result;
    }

    /** Create the map of `ChrName` to `ChrInfo` and while doing that check for duplicate
      * chromosome names. */
    boost::unordered::unordered_map<ChrName, ChrInfo>
    ChrInfoTable::buildChrInfosByName(const std::vector<ChrInfo> &chr_info) {
        boost::unordered::unordered_map<ChrName, ChrInfo> result;
        result.reserve(chr_info.size());
        for (const auto &info : chr_info) {
            if (result.contains(info.getName())) {
                throw_with_trace(std::invalid_argument("Duplicate chromosome name '" +
                                                       info.getName() + "'"));
            }
            result.try_emplace(info.getName(), info);
        }
        return result;
    }

    ChrInfoTable::ChrInfoTable(const std::vector<ChrInfo> &chr_infos)
            : chrInfos { chr_infos },
              chrInfosByCategory { buildChrInfosByCategory(chr_infos) },
              chrInfosByName { buildChrInfosByName(chr_infos) } {

    }

    std::vector<ChrInfo>::size_type ChrInfoTable::nChromosomes() const {
        return chrInfos.size();
    }

    const std::vector<ChrInfo>& ChrInfoTable::getChrInfos() const {
        return chrInfos;
    }

    const std::vector<ChrInfo> &ChrInfoTable::getChrInfos(ChrCategory category) const {
        return chrInfosByCategory.at(category);
    }

    ChrInfoTable::ChrNames ChrInfoTable::getNames() const {
        ChrNames result;
        std::transform(chrInfos.begin(),
                       chrInfos.end(),
                       std::back_inserter(result),
                       [](const ChrInfo &info) { return info.getName(); });
        return result;
    }

    ChrInfoTable::ChrNames ChrInfoTable::getNames(ChrCategory category) const {
        ChrNames result;
        std::transform(chrInfosByCategory.at(category).begin(),
                       chrInfosByCategory.at(category).end(),
                       std::back_inserter(result),
                       [](const ChrInfo &info) { return info.getName(); });
        return result;
    }

    ChrInfoTable::ChrSizes ChrInfoTable::getSizes() const {
        ChrSizes result;
        std::transform(chrInfos.begin(),
                       chrInfos.end(),
                       std::back_inserter(result),
                       [](const ChrInfo &info) { return info.getSize(); });
        return result;
    }

    ChrInfoTable::ChrSizes ChrInfoTable::getSizes(ChrCategory category) const {
        ChrSizes result;
        std::transform(chrInfosByCategory.at(category).begin(),
                       chrInfosByCategory.at(category).end(),
                       std::back_inserter(result),
                       [](const ChrInfo &info) { return info.getSize(); });
        return result;
    }

} // namespace sophia
