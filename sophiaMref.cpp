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


int main(int argc, char** argv) {
    using namespace std;

	boost::program_options::options_description desc("Allowed options");
	desc.add_options() //
	("help", "produce help message") //
	("gzins", boost::program_options::value<string>(), "list of all gzipped control beds") //
	("version", boost::program_options::value<string>(), "version") //
	("defaultreadlength", boost::program_options::value<int>(), "Default read length for the technology used in sequencing 101,151 etc.") //
	("outputrootname", boost::program_options::value<string>(), "outputrootname");
	boost::program_options::variables_map inputVariables { };
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), inputVariables);
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
	int defaultReadLength { 0 };
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
	sophia::SuppAlignment::DEFAULTREADLENGTH = defaultReadLength;
	sophia::SuppAlignmentAnno::DEFAULTREADLENGTH = defaultReadLength;
	sophia::MrefEntry::NUMPIDS = gzListIn.size();
	sophia::MasterRefProcessor mRefProcessor { gzListIn, outputRoot, version, defaultReadLength };
}
