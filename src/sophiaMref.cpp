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
#include "Hg38ChrConverter.h"
#include "HelperFunctions.h"


int main(int argc, char** argv) {
    using namespace std;
    using namespace sophia;

    int defaultReadLength { 0 };
    string assemblyName = "hg37";

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
                boost::program_options::value<string>(&assemblyName)->default_value("hg37"),
                "assembly name (hg37, hg38)")
            ("defaultreadlength",
                boost::program_options::value<int>(&defaultReadLength)->default_value(defaultReadLength),
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
            defaultReadLength = inputVariables["defaultreadlength"].as<int>();
        } else {
            cerr << "Default read Length not given, exiting" << endl;
            return 1;
        }

        string outputRoot { };
        if (inputVariables.count("outputrootname")) {
            outputRoot = inputVariables["outputrootname"].as<string>();
        } else {
            cerr << "No output file root name given, exiting" << endl;
            return 1;
        }

        unique_ptr<ChrConverter> chrConverter;
        if (!inputVariables.count("assemblyname") ||
              inputVariables["assemblyname"].as<string>() == Hg37ChrConverter::assemblyName) {
            chrConverter = unique_ptr<ChrConverter>(new Hg37ChrConverter());
        } else if (inputVariables["assemblyname"].as<string>() == Hg38ChrConverter::assemblyName) {
            chrConverter = unique_ptr<ChrConverter>(new Hg38ChrConverter());
        } else {
            cerr << "Unknown assembly name " << inputVariables["assemblyname"].as<string>() << ". I know "
                 << Hg37ChrConverter::assemblyName << " and "
                 << Hg38ChrConverter::assemblyName << endl;
                return 1;
        }

        // Initialize the global application configuration.
        GlobalAppConfig::init(move(chrConverter));

        SuppAlignment::DEFAULTREADLENGTH = defaultReadLength;
        SuppAlignmentAnno::DEFAULTREADLENGTH = defaultReadLength;
        MrefEntry::NUMPIDS = gzListIn.size();
        MasterRefProcessor mRefProcessor { gzListIn, outputRoot, version, defaultReadLength };

        return 0;
    } catch (boost::exception &e) {
        cerr << get_trace(e) << endl;
        return 1;
    } catch (std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}
