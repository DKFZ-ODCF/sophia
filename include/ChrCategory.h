#ifndef CHRCATEGORY_H_
#define CHRCATEGORY_H_


#include <string>
#include <iterator>
#include <boost/unordered/unordered_map.hpp>

namespace sophia {

    /** C++ enums suck (also enum class). This is just a manually implemented rich enum. */
    class ChrCategory {

      public:

        using size_type =
            boost::unordered::unordered_map<std::string, const ChrCategory>::size_type;

      private:

        std::string category_name;

        // This is mostly for comparison with `operator<`.
        std::size_t category_index;

        // Used for initialization.
        static const boost::unordered::unordered_map<std::string, const ChrCategory> categories;

        static const std::vector<ChrCategory> sorted_categories;

      public:
        // Only used to define categories.
        ChrCategory(const std::string &s, std::size_t index);

        // Predefined instances.


        /** Autosomal contigs, e.g. chr1, chr2, ..., chr22 */
        static const ChrCategory& AUTOSOME;
        /** X chromosome */
        static const ChrCategory& X;
        /** Y chromosome */
        static const ChrCategory& Y;
        /** extrachromosomalContigs  Extrachromosomal contigs, e.g. chrM, chrMT */
        static const ChrCategory& EXTRACHROMOSOMAL;
        /** Joined category for unlocalized, unplaced, or random placed contigs
          * Compare https://www.ncbi.nlm.nih.gov/grc/help/definitions
          */
        static const ChrCategory& UNASSIGNED;
        /** _alt contigs */
        static const ChrCategory& ALT;
        /** HLA contigs */
        static const ChrCategory& HLA;
        /** Virus contigs, e.g. NC_007605, EBV.
          * This is for viruses that may insert into the nuclear genome. */
        static const ChrCategory& VIRUS;
        /** Decoy contigs. */
        static const ChrCategory& DECOY;
        /** Technical contigs, e.g. phiX or lambda */
        static const ChrCategory& TECHNICAL;

        // Parser
        static const ChrCategory& from_string(const std::string &s);

        ~ChrCategory();

        static size_type numCategories();

        static const std::vector<ChrCategory>& getCategories();

        std::string getName() const;

        bool operator==(const ChrCategory &other) const;

        bool operator!=(const ChrCategory &other) const;

        bool operator<(const ChrCategory &other) const;

        bool operator>(const ChrCategory &other) const;

    };


} // namespace sophia

/** The following defines hash and equal_to functions such that ChrCategory can be used as a key in
 *  unordered containers without explicitly setting these two functions. Thus we can continue to
 *  use `boost::unordered::unordered_set<ChrCategory>`. */

namespace boost {

    template<>
    struct hash<sophia::ChrCategory> {
        std::size_t operator()(const sophia::ChrCategory& chrCategory) const {
            return std::hash<std::string>()(chrCategory.getName());
        }
    };

} // namespace boost

namespace std {

    template<>
    struct equal_to<sophia::ChrCategory> {
        bool operator()(const sophia::ChrCategory& lhs, const sophia::ChrCategory& rhs) const {
            return lhs.operator==(rhs);
        }
    };

    template<>
    struct less<sophia::ChrCategory> {
        bool operator()(const sophia::ChrCategory& lhs, const sophia::ChrCategory& rhs) const {
            return lhs.operator<(rhs);
        }
    };

} // namespace std

#endif /* CHRCATEGORY_H_ */
