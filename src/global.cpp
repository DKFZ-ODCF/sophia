#include "global.h"
#include <exception>
#include <sstream>
#include <iostream>
#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>

namespace sophia {

    std::string get_trace(const boost::exception &e) {
        std::stringstream ss;
        const boost::stacktrace::stacktrace *st = boost::get_error_info<traced>(e);
        if (st != nullptr) {
            ss << *st << std::endl;
        }
        return ss.str();
    }

} // namespace sophia