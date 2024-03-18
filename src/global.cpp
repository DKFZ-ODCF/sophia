#include <exception>
#include <sstream>
#include <iostream>
#include <optional>
#include <string>
#include <memory>
#include <boost/stacktrace.hpp>
#include <boost/exception/all.hpp>
#include "global.h"
#include "ChrConverter.h"
#include "ChrInfo.h"
#include "ChrInfoTable.h"
#include "Hg37ChrConverter.h"
#include "GenericChrConverter.h"
#include "GlobalAppConfig.h"


namespace sophia {

    std::string get_trace(const boost::exception &e) {
        std::stringstream ss;
        const boost::stacktrace::stacktrace *st = boost::get_error_info<traced>(e);
        if (st != nullptr) {
            ss << *st << std::endl;
        }
        return ss.str();
    }

    void setApplicationConfig(std::optional<std::string> assembly_name) {
        std::unique_ptr<ChrConverter> converter;

        if (!assembly_name.has_value() || assembly_name.value() == "classic_hg37") {
            converter = std::unique_ptr<ChrConverter>(new Hg37ChrConverter());

        } else if (assembly_name.value().size() > 0) {
            // Compose the chromosome set name from the assembly name.
            std::string chromosome_file = "resources/" + assembly_name.value() + ".tsv";
            std::vector<ChrInfo> chr_info = read_chr_info(chromosome_file);
            ChrInfoTable chr_info_table { chr_info };
            converter = std::unique_ptr<ChrConverter>(
                new GenericChrConverter(assembly_name.value(), chr_info_table));
        } else {
            throw std::invalid_argument("Empty assembly name.");
        }

        // Initialize the global application configuration.
        GlobalAppConfig::init(std::move(converter));

    }

} // namespace sophia