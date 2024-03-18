#ifndef CHRINFO_H_
#define CHRINFO_H_

#include "global.h"
#include "ChrCategory.h"
#include <vector>
#include <string>
#include <iostream>
#include <boost/unordered/unordered_map.hpp>

namespace sophia {

    class ChrInfo {

    private:
        ChrName name;
        ChrSize size;
        bool compressedMref;
        ChrCategory category;

    public:
        ChrInfo(ChrName _name,
                ChrSize _size,
                bool _compressedMref,
                ChrCategory _category);

        ChrName getName() const;
        ChrSize getSize() const;
        bool isCompressedMref() const;
        ChrCategory getCategory() const;

    };



    /** Read the chromosome information from a TSV table with a header (chromosome, size, class),
        where class is defined for each chromosome ID as "primary", "extrachromosomal", "random",
        "unplaced", "alt", "hla", "decoy", "technical". All classes but "primary" are allowed to
        be empty. */
    bool to_boolean(const std::string& str);
    std::vector<ChrInfo> read_chr_info(std::istream &in);
    std::vector<ChrInfo> read_chr_info(const std::string &filename);

    /** Convert a sequence of ChrInfo objects into a map from category to ChrInfo. */
    boost::unordered::unordered_map<ChrCategory, std::vector<ChrInfo>>
    to_chr_info_map(const std::vector<ChrInfo> &chr_info);


} // namespace sophia

#endif /* CHRINFO_H_ */
