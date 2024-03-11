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
#include "HelperFunctions.h"
#include "ChrInfo.h"
#include "ChrInfoTable.h"
#include "ChrCategory.h"


int main(int argc, char** argv) {
    using namespace std;
    using namespace sophia;

    ChrSize defaultReadLength { 0 };
    string assemblyName = "classic_hg37";

    try {
        boost::program_options::options_description desc("Allowed options for sophiaMref");
            desc.add_options()
            ("help",
                "produce help message")
            ("gzins",
                boost::program_options::value<string>(),
                "A file containing the the paths of the of all gzipped control beds, line-by-line")
            ("outputrootname",
                boost::program_options::value<string>(),
                "base name/path for the output files")
            ("version",
                boost::program_options::value<string>(),
                "version string used to match the PID in the BED files with the pattern\n  `.*/$pidName.{1}$version.+`")
            ("assemblyname",
                boost::program_options::value<string>(&assemblyName)->default_value("classic_hg37"),
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
            cout << desc << endl;
            return 0;
        }

        std::optional<std::string> assemblyNameOpt { };
        if (inputVariables.count("assemblyname")) {
            assemblyNameOpt = inputVariables["assemblyname"].as<string>();
        }
        setApplicationConfig(assemblyNameOpt);

        string gzInFilesList;
        if (inputVariables.count("gzins")) {
            gzInFilesList = inputVariables["gzins"].as<string>();
        } else {
            cerr << "No gzipped control bed list file given, exiting" << endl;
            return 1;
        }

        ifstream gzInFilesHandle { gzInFilesList };
        vector<string> gzListIn;
        for (string line; error_terminating_getline(gzInFilesHandle, line);) {
            gzListIn.push_back(line);
        }

        string version { };
        if (inputVariables.count("version")) {
            version = inputVariables["version"].as<string>();
        } else {
            cerr << "No input version given, exiting" << endl;
            return 1;
        }

        if (inputVariables.count("defaultreadlength")) {
            defaultReadLength = inputVariables["defaultreadlength"].as<ChrSize>();
        } else {
            cerr << "Default read length not given, exiting" << endl;
            return 1;
        }
        if (defaultReadLength < 1) {
            cerr << "Default read length " << std::to_string(defaultReadLength)
                 << " is invalid." << endl;
            return 1;
        }


        string outputRoot { };
        if (inputVariables.count("outputrootname")) {
            outputRoot = inputVariables["outputrootname"].as<string>();
        } else {
            cerr << "No output file root name given, exiting" << endl;
            return 1;
        }

        SuppAlignment::DEFAULT_READ_LENGTH = defaultReadLength;
        SuppAlignmentAnno::DEFAULT_READ_LENGTH = defaultReadLength;
        MrefEntry::NUM_PIDS = gzListIn.size();

        cerr << "Running sophiaMref on " << MrefEntry::NUM_PIDS << " PIDs ..." << endl;
        MasterRefProcessor mRefProcessor { gzListIn, outputRoot, version, defaultReadLength };

        return 0;
    } catch (boost::exception &e) {
        cerr << "Error: " << boost::diagnostic_information(e) << endl;
        return 1;
    } catch (std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}
