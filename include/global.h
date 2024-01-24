/*
 * GlobalAppConfig.h
 *
 * Author: Philip R. Kensche Copyright (C) 2023 DKFZ Heidelberg
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *     LICENSE: GPL
 */

#ifndef GLOBALAPPCONFIG_H_
#define GLOBALAPPCONFIG_H_

#include <cstddef>
#include <iostream>
#include <optional>
#include <string>
#include <boost/exception/all.hpp>
#include <boost/stacktrace.hpp>


namespace sophia {

    // Enrich boost::exceptions with additional std::string information
    using error_info_string = boost::error_info<struct tag_string, std::string>;

    class DomainError :
            public virtual boost::exception,
            public virtual std::domain_error {
      public:
        DomainError(const std::string &msg) : std::domain_error(msg) {}
    };

    // These two are only to make the code clearer, but are not type checked. There are no opaque
    // or strongly type-checked type-"aliases" in C++17. A type-safe solution would use classes.
    using ChrName = std::string;
    using ChrIndex = size_t;
    using CompressedMrefIndex = size_t ;
    using ChrSize = long unsigned int;

    std::string get_trace(const boost::exception &e);

    typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> traced;

    template <class E>
    void throw_with_trace(const E &e) {
        throw boost::enable_error_info(e) <<
              traced(boost::stacktrace::stacktrace());
    }

    void setApplicationConfig(std::optional<std::string> assemblyname);
}


#endif /* GLOBALAPPCONFIG_H_ */