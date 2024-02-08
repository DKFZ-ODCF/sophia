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
	boost::program_options::options_description desc("Allowed options");
	desc.add_options() //
	("help", "produce help message") //
	("gzins", boost::program_options::value<std::string>(), "list of all gzipped control beds") //
	("version", boost::program_options::value<std::string>(), "version") //
	("defaultreadlength", boost::program_options::value<int>(), "Default read length for the technology used in sequencing 101,151 etc.") //
	("outputrootname", boost::program_options::value<std::string>(), "outputrootname");
	boost::program_options::variables_map inputVariables { };
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), inputVariables);
	boost::program_options::notify(inputVariables);
	if (inputVariables.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}
	std::string gzInFilesList;
	if (inputVariables.count("gzins")) {
		gzInFilesList = inputVariables["gzins"].as<std::string>();
	} else {
		std::cerr << "No gzipped control bed list file given, exiting" << std::endl;
		return 1;
	}
	std::ifstream gzInFilesHandle { gzInFilesList };
	std::vector<std::string> gzListIn;
	for (std::string line; sophia::error_terminating_getline(gzInFilesHandle, line);) {
		gzListIn.push_back(line);
	}
	std::string version { };
	if (inputVariables.count("version")) {
		version = inputVariables["version"].as<std::string>();
	} else {
		std::cerr << "No input version given, exiting" << std::endl;
		return 1;
	}
	int defaultReadLength { 0 };
	if (inputVariables.count("defaultreadlength")) {
		defaultReadLength = inputVariables["defaultreadlength"].as<int>();
	} else {
		std::cerr << "Default read Length not given, exiting" << std::endl;
		return 1;
	}
	std::string outputRoot { };
	if (inputVariables.count("outputrootname")) {
		outputRoot = inputVariables["outputrootname"].as<std::string>();
	} else {
		std::cerr << "No output file root name given, exiting" << std::endl;
		return 1;
	}
	sophia::SuppAlignment::DEFAULTREADLENGTH = defaultReadLength;
	sophia::SuppAlignmentAnno::DEFAULTREADLENGTH = defaultReadLength;
	sophia::MrefEntry::NUMPIDS = gzListIn.size();
	sophia::MasterRefProcessor mRefProcessor { gzListIn, outputRoot, version, defaultReadLength };
}
