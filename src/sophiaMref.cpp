/*
 * sophiaMref.cpp
 *
 *  Created on: 27 Apr 2016
 *      Author: umuttoprak
 */

#include <string>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <vector>
#include "MasterRefProcessor.h"
#include "MrefEntry.h"
#include "HelperFunctions.h"
#include "GlobalAppConfig.h"
#include "Hg37ChrConverter.h"
#include "GenericChrConverter.h"
#include "ChrInfo.h"
#include "ChrInfoTable.h"
#include "ChrCategory.h"


int main(int argc, char** argv) {

    using namespace sophia;

    ChrSize defaultReadLength { 0 };
    std::string assemblyName = "classic_hg37";

    try {
        boost::program_options::options_description desc("Allowed options for sophiaMref");
            desc.add_options()
            ("help",
                "produce help message")
            ("gzins",
                boost::program_options::value<std::string>(),
                "A file containing the the paths of the of all gzipped control beds, line-by-line")
            ("outputrootname",
                boost::program_options::value<std::string>(),
                "base name/path for the output files")
            ("version",
                boost::program_options::value<std::string>(),
                "version string used to match the PID in the BED files with the pattern\n  `.*/$pidName.{1}$version.+`")
            ("assemblyname",
                boost::program_options::value<std::string>(&assemblyName)->default_value("classic_hg37"),
                "assembly name (classic_hg37, hg38, ...)")
            ("defaultreadlength",
                boost::program_options::value<ChrSize>(&defaultReadLength),
                "Default read length for the technology used in sequencing, e.g. 101 or 151.")
        ;

        boost::program_options::variables_map inputVariables { };
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc),
                                      inputVariables);
        boost::program_options::notify(inputVariables);

        if (inputVariables.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }

        std::optional<std::string> assemblyNameOpt { };
        if (inputVariables.count("assemblyname")) {
            assemblyNameOpt = inputVariables["assemblyname"].as<std::string>();
        }
        setApplicationConfig(assemblyNameOpt);

        std::string gzInFilesList;
        if (inputVariables.count("gzins")) {
            gzInFilesList = inputVariables["gzins"].as<std::string>();
        } else {
            std::cerr << "No gzipped control bed list file given, exiting" << std::endl;
            return 1;
        }

        std::ifstream gzInFilesHandle { gzInFilesList };
        std::vector<std::string> gzListIn;
        for (std::string line; error_terminating_getline(gzInFilesHandle, line);) {
            gzListIn.push_back(line);
        }

        std::string version { };
        if (inputVariables.count("version")) {
            version = inputVariables["version"].as<std::string>();
        } else {
            std::cerr << "No input version given, exiting" << std::endl;
            return 1;
        }

        if (inputVariables.count("defaultreadlength")) {
            defaultReadLength = inputVariables["defaultreadlength"].as<ChrSize>();
        } else {
            std::cerr << "Default read length not given, exiting" << std::endl;
            return 1;
        }
        if (defaultReadLength < 1) {
            std::cerr << "Default read length " << std::to_string(defaultReadLength)
                 << " is invalid." << std::endl;
            return 1;
        }


        std::string outputRoot { };
        if (inputVariables.count("outputrootname")) {
            outputRoot = inputVariables["outputrootname"].as<std::string>();
        } else {
            std::cerr << "No output file root name given, exiting" << std::endl;
            return 1;
        }

        SuppAlignment::DEFAULT_READ_LENGTH = defaultReadLength;
        SuppAlignmentAnno::DEFAULT_READ_LENGTH = defaultReadLength;
        MrefEntry::NUM_PIDS = gzListIn.size();

        std::cerr << "Running sophiaMref on " << MrefEntry::NUM_PIDS << " PIDs ..." << std::endl;
        MasterRefProcessor mRefProcessor { gzListIn, outputRoot, version, defaultReadLength };

        return 0;
    } catch (boost::exception &e) {
        std::cerr << "Error: " << boost::diagnostic_information(e) << std::endl;
        return 1;
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
