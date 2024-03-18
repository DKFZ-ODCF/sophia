#include "global.h"
#include "ChrInfo.h"
#include "rapidcsv.h"
#include <vector>
#include <stdexcept>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
#include <boost/algorithm/string.hpp>


namespace sophia {

    ChrInfo::ChrInfo(ChrName _name,
                     ChrSize _size,
                     bool _compressedMref,
                     ChrCategory _category)
            : name {_name},
              size {_size},
              compressedMref {_compressedMref},
              category {_category} {
        if (name.size() == 0) {
            throw std::invalid_argument("ChrInfo: name cannot be empty");
        }
        if (size <= 0) {
            throw std::invalid_argument("ChrInfo: size must be larger than zero for '" + name + "'");
        }
    }

    ChrName ChrInfo::getName() const {
        return name;
    }

    ChrSize ChrInfo::getSize() const {
        return size;
    }

    bool ChrInfo::isCompressedMref() const {
        return compressedMref;
    }

    ChrCategory ChrInfo::getCategory() const {
        return category;
    }


    bool to_boolean (const std::string& str) {
        static const boost::unordered::unordered_set<std::string> trues {
            "true",
            "t",
            "1",
            "no",
            "n"
        };
        static const boost::unordered::unordered_set<std::string> falses {
            "false",
            "f",
            "0",
            "yes",
            "y"
        };
        std::string lowercaseStr = str;
        boost::algorithm::to_lower(lowercaseStr);
        bool result;
        if (trues.contains(lowercaseStr)) {
            result = true;
        } else if (falses.contains(lowercaseStr)) {
            result = false;
        } else {
            throw_with_trace(std::invalid_argument("Could not parse boolean from '" + str + "'"));
        }
        return result;
    }

    std::vector<ChrInfo> read_chr_info(std::istream &in) {
        rapidcsv::Document doc(in,
                               rapidcsv::LabelParams(0, -1),
                               rapidcsv::SeparatorParams('\t'));

        // Used for checking uniqueness of chromosome names
        boost::unordered::unordered_set<ChrName> names;
        names.reserve(doc.GetRowCount());

        // Convert the file into a vector of ChrInfo.
        std::vector<ChrInfo> chr_info;
        chr_info.reserve(doc.GetRowCount());
        for (size_t i = 0; i < doc.GetRowCount(); ++i) {
            ChrName name = doc.GetCell<std::string>("chromosome", i);
            if (names.find(name) != names.end()) {
                throw_with_trace(std::invalid_argument(
                    "Chromosome name '" + name + "' is not unique."));
            }
            ChrSize size = doc.GetCell<ChrSize>("size", i);
            std::string category_string = doc.GetCell<std::string>("category", i);
            ChrCategory category = ChrCategory::from_string(category_string);
            bool compressedMref = to_boolean(doc.GetCell<std::string>("compressedMref", i));
            chr_info.emplace_back(ChrInfo(name, size, compressedMref, category));
        }

        return chr_info;
    }


    std::vector<ChrInfo> read_chr_info(const std::string &filename) {
        std::ifstream in(filename);
        if (!in) {
            throw_with_trace(std::invalid_argument("Cannot open file '" + filename + "'"));
        }
        return read_chr_info(in);
    }

} // namespace sophia
