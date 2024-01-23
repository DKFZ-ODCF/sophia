#ifndef CHRINFOTABLE_H_
#define CHRINFOTABLE_H_

#include "global.h"
#include "ChrInfo.h"
#include "ChrInfoTable.h"
#include <boost/unordered/unordered_map.hpp>

namespace sophia {

    /** Provide access to chromosome information.
      * Note that the order of the input vector is preserved. This means, if the vector is not
      * sorted by the category, then also the results will not be sorted.
      */
    class ChrInfoTable {

      public:
        using ChrNames = std::vector<ChrName>;
        using ChrSizes = std::vector<ChrSize>;

      private:
        const std::vector<ChrInfo> chrInfos;

        const
        boost::unordered::unordered_map<ChrCategory, std::vector<ChrInfo>>
        chrInfosByCategory;

        /** Helper for the constructor. */
        static
        boost::unordered::unordered_map<ChrCategory, std::vector<ChrInfo>>
        buildChrInfosByCategory(const std::vector<ChrInfo> &chr_info);

        const boost::unordered::unordered_map<ChrName, ChrInfo> chrInfosByName;

        /** Helper for the constructor. */
        static
        boost::unordered::unordered_map<ChrName, ChrInfo>
        buildChrInfosByName(const std::vector<ChrInfo> &chr_info);

      public:

        ChrInfoTable(const std::vector<ChrInfo> &chr_info);

        std::vector<ChrInfo>::size_type nChromosomes() const;

        /** Return the `ChrInfo` in the same order they were provided to the constructor. */
        const std::vector<ChrInfo> &getChrInfos() const;
        const std::vector<ChrInfo> &getChrInfos(ChrCategory category) const;

        ChrNames getNames() const;

        /** Get the names of all chromosomes of the given category. The order of values is exactly
          * the same as they were (maybe interrupted by other chromosomes) in the vector provided
          * to the constructor.
          */
        ChrNames getNames(ChrCategory category) const;

        /** The the chromosome sizes, again in the same order as provided to the constructor. */
        ChrSizes getSizes() const;

        /** Get the lengths of all chromosomes of the given category. The order of values is exactly
          * the same as they were (maybe interrupted by other chromosomes) in the vector provided
          * to the constructor.
          */
        ChrSizes getSizes(ChrCategory category) const;

        ChrInfo getChrInfo(ChrName name) const;

    };


} /* namespace sophia */

#endif /* CHRINFOTABLE_H_ */
