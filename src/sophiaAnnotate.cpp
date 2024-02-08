/*
 * sophiaAnnotate.cpp
 *
 *  Created on: 28 Apr 2016
 *      Author: umuttoprak
 */
#include <iostream>
#include <cmath>
#include <string>
#include "cxxopts.hpp"
#include "BreakpointReduced.h"
#include "AnnotationProcessor.h"
#include "SuppAlignment.h"
#include "SuppAlignmentAnno.h"
#include "SvEvent.h"
#include "strtk.hpp"
#include <boost/filesystem.hpp>
#include <vector>
#include "MrefEntryAnno.h"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "ChrConverter.h"
#include "HelperFunctions.h"

int main(int argc, char** argv) {
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(nullptr);
	cxxopts::Options options("SophiaAnnotate", "Annotates SOPHIA output");
	options.add_options() //
	("help", "produce help message") //
	("mref", "mref file", cxxopts::value<std::string>()) //
	("tumorresults", "_bps.bed.gz file from sophia for the tumor, or control for a no-tumor analysis", cxxopts::value<std::string>()) //
	("controlresults", "_bps.bed.gz file from sophia for the control", cxxopts::value<std::string>()) //
	("defaultreadlengthtumor", "Default read length for the technology used in sequencing 101,151 etc., tumor", cxxopts::value<int>()) //
	("defaultreadlengthcontrol", "Default read length for the technology used in sequencing 101,151 etc., tumor", cxxopts::value<int>()) //
	("pidsinmref", "Number of PIDs in the MREF", cxxopts::value<int>()) //
	("artifactlofreq", "PERCENTAGE frequency of artifact supports for bps to be considered as artifact_like (33)", cxxopts::value<int>()) //
	("artifacthifreq", "PERCENTAGE frequency of artifact supports for bps to be considered as artifacts (50)", cxxopts::value<int>()) //
	("clonalitylofreq", "PERCENTAGE clonality for bps to be considered as extreme_subclonal (10)", cxxopts::value<int>()) //
	("clonalitystrictlofreq", "PERCENTAGE clonality for bps to be considered as extreme_subclonal (20)", cxxopts::value<int>()) //
	("clonalityhifreq", "PERCENTAGE clonality for bps to be considered as homozygous (85)", cxxopts::value<int>()) //
	("bpfreq", "PERCENTAGE frequency of a BP for consideration as rare. (3)", cxxopts::value<int>()) //
	("germlineoffset", "Minimum offset a germline bp and a control bp. (5)", cxxopts::value<int>()) //
	("germlinedblimit", "Maximum occurrence of germline variants in the db. (5)", cxxopts::value<int>()) //
	("debugmode", "debugmode");
	options.parse(argc, argv);
	std::vector<std::vector<sophia::MrefEntryAnno>> mref { 85, std::vector<sophia::MrefEntryAnno> { } };
	if (!options.count("mref")) {
		std::cerr << "No mref file given, exiting" << std::endl;
		return 1;
	}

	std::string tumorResults;
	if (options.count("tumorresults")) {
		tumorResults = options["tumorresults"].as<std::string>();
	} else {
		std::cerr << "No input file given, exiting" << std::endl;
		return 1;
	}
	int pidsInMref { 0 };
	if (options.count("pidsinmref")) {
		pidsInMref = options["pidsinmref"].as<int>();
	} else {
		std::cerr << "number of PIDS in the MREF not given, exiting" << std::endl;
		return 1;
	}
	int defaultReadLengthTumor { 0 };
	if (options.count("defaultreadlengthtumor")) {
		defaultReadLengthTumor = options["defaultreadlengthtumor"].as<int>();
	} else {
		std::cerr << "Default read Length not given, exiting" << std::endl;
		return 1;
	}
	int artifactlofreq { 33 };
	if (options.count("artifactlofreq")) {
		artifactlofreq = options["artifactlofreq"].as<int>();
	}
	int artifacthifreq { 50 };
	if (options.count("artifacthifreq")) {
		artifacthifreq = options["artifacthifreq"].as<int>();
	}
	int clonalitylofreq { 5 };
	if (options.count("clonalitylofreq")) {
		clonalitylofreq = options["clonalitylofreq"].as<int>();
	}
	int clonalitystrictlofreq { 20 };
	if (options.count("clonalitystrictlofreq")) {
		clonalitystrictlofreq = options["clonalitystrictlofreq"].as<int>();
	}
	int clonalityhifreq { 85 };
	if (options.count("clonalityhifreq")) {
		clonalityhifreq = options["clonalityhifreq"].as<int>();
	}
	int bpFreq { 3 };
	if (options.count("bpfreq")) {
		bpFreq = options["bpfreq"].as<int>();
	}
	int germlineOffset { 5 };
	if (options.count("germlineoffset")) {
		germlineOffset = options["germlineoffset"].as<int>();
	}
	int germlineDbLimit { 5 };
	if (options.count("germlinedblimit")) {
		germlineDbLimit = options["germlinedblimit"].as<int>();
	}
	sophia::MrefEntryAnno::PIDSINMREF = pidsInMref;
	std::unique_ptr<std::ifstream> mrefInputHandle { std::make_unique<std::ifstream>(options["mref"].as<std::string>(), std::ios_base::in | std::ios_base::binary) };
	std::unique_ptr<boost::iostreams::filtering_istream> mrefGzHandle { std::make_unique<boost::iostreams::filtering_istream>() };
	mrefGzHandle->push(boost::iostreams::gzip_decompressor());
	mrefGzHandle->push(*mrefInputHandle);
	std::cerr << "m\n";
	std::string line { };
	while (sophia::error_terminating_getline(*mrefGzHandle, line)) {
		if (line.front() == '#') {
			continue;
		};
		auto chrIndex = sophia::ChrConverter::indexConverter[sophia::ChrConverter::readChromosomeIndex(line.cbegin(), '\t')];
		if (chrIndex < 0) {
			continue;
		}
		mref[chrIndex].emplace_back(line);
	}
	sophia::SvEvent::ARTIFACTFREQLOWTHRESHOLD = (artifactlofreq + 0.0) / 100;
	sophia::SvEvent::ARTIFACTFREQHIGHTHRESHOLD = (artifacthifreq + 0.0) / 100;
	sophia::BreakpointReduced::ARTIFACTFREQHIGHTHRESHOLD = sophia::SvEvent::ARTIFACTFREQHIGHTHRESHOLD;
	sophia::SvEvent::CLONALITYLOWTHRESHOLD = (clonalitylofreq + 0.0) / 100;
	sophia::SvEvent::CLONALITYSTRICTLOWTHRESHOLD = (clonalitystrictlofreq + 0.0) / 100;
	sophia::BreakpointReduced::CLONALITYSTRICTLOWTHRESHOLD = sophia::SvEvent::CLONALITYSTRICTLOWTHRESHOLD;
	sophia::SvEvent::CLONALITYHIGHTHRESHOLD = (clonalityhifreq + 0.0) / 100;
	sophia::SvEvent::BPFREQTHRESHOLD = pidsInMref * (bpFreq + 0.0) / 100;
	sophia::SvEvent::RELAXEDBPFREQTHRESHOLD = 3 * sophia::SvEvent::BPFREQTHRESHOLD;
	sophia::SvEvent::PIDSINMREFSTR = strtk::type_to_string<int>(pidsInMref);
	sophia::BreakpointReduced::PIDSINMREFSTR = sophia::SvEvent::PIDSINMREFSTR;
	sophia::BreakpointReduced::DEFAULTREADLENGTH = defaultReadLengthTumor;
	sophia::Breakpoint::DEFAULTREADLENGTH = defaultReadLengthTumor;
	sophia::SuppAlignment::DEFAULTREADLENGTH = defaultReadLengthTumor;
	sophia::SuppAlignmentAnno::DEFAULTREADLENGTH = defaultReadLengthTumor;
	sophia::SvEvent::HALFDEFAULTREADLENGTH = std::round(defaultReadLengthTumor / 2.0);
	sophia::SvEvent::GERMLINEOFFSETTHRESHOLD = germlineOffset;
	sophia::SvEvent::GERMLINEDBLIMIT = germlineDbLimit;
	sophia::SvEvent::ABRIDGEDOUTPUT = true;
	if (options.count("debugmode")) {
		sophia::SvEvent::DEBUGMODE = true;
	} else {
		sophia::SvEvent::DEBUGMODE = false;
	}
	sophia::AnnotationProcessor::ABRIDGEDOUTPUT = true;
	sophia::Breakpoint::BPSUPPORTTHRESHOLD = 3;
	if (options.count("controlresults")) {
		std::string controlResults { options["controlresults"].as<std::string>() };
		int defaultReadLengthControl { 0 };
		if (options.count("defaultreadlengthcontrol")) {
			defaultReadLengthControl = options["defaultreadlengthtumor"].as<int>();
		} else {
			std::cerr << "Default read Length not given, exiting" << std::endl;
			return 1;
		}
		auto lowQualControl = 0;
		auto pathogenInControl = false;
		{
			sophia::SvEvent::NOCONTROLMODE = true;
			sophia::AnnotationProcessor annotationProcessorControlCheck { controlResults, mref, defaultReadLengthControl, true, germlineDbLimit };
			lowQualControl = annotationProcessorControlCheck.getMassiveInvFilteringLevel();
			pathogenInControl = annotationProcessorControlCheck.isContaminationObserved();
			sophia::SvEvent::NOCONTROLMODE = false;
		}
		sophia::AnnotationProcessor annotationProcessor { tumorResults, mref, controlResults, defaultReadLengthTumor, defaultReadLengthControl, germlineDbLimit, lowQualControl, pathogenInControl };
		annotationProcessor.printFilteredResults(pathogenInControl, lowQualControl);
	} else {
		sophia::SvEvent::NOCONTROLMODE = true;
		sophia::AnnotationProcessor annotationProcessor { tumorResults, mref, defaultReadLengthTumor, false, germlineDbLimit };
		annotationProcessor.printFilteredResults(false, 0);
	}
	return 0;
}
