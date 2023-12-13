/*
 * sophiaAnnotate.cpp
 *
 *  Created on: 28 Apr 2016
 *      Author: umuttoprak
 */
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
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
#include "HelperFunctions.h"
#include "ChrConverter.h"
#include "Hg37ChrConverter.h"
#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"


int main(int argc, char** argv) {
    using namespace std;
    using namespace sophia;

	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);
	cxxopts::Options options("SophiaAnnotate", "Annotates SOPHIA output");
	options.add_options() //
	("help", "produce help message") //
	("mref", "mref file", cxxopts::value<string>()) //
	("assemblyname", "assembly name", cxxopts::value<string>()) //
	("tumorresults", "_bps.bed.gz file from sophia for the tumor, or control for a no-tumor analysis", cxxopts::value<string>()) //
	("controlresults", "_bps.bed.gz file from sophia for the control", cxxopts::value<string>()) //
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

	vector<vector<MrefEntryAnno>> mref { 85, vector<MrefEntryAnno> { } };
	if (!options.count("mref")) {
		cerr << "No mref file given, exiting" << endl;
		return 1;
	}

	string tumorResults;
	if (options.count("tumorresults")) {
		tumorResults = options["tumorresults"].as<string>();
	} else {
		cerr << "No input file given, exiting" << endl;
		return 1;
	}

	int pidsInMref { 0 };

	if (options.count("pidsinmref")) {
		pidsInMref = options["pidsinmref"].as<int>();
	} else {
		cerr << "number of PIDS in the MREF not given, exiting" << endl;
		return 1;
	}

	int defaultReadLengthTumor { 0 };
	if (options.count("defaultreadlengthtumor")) {
		defaultReadLengthTumor = options["defaultreadlengthtumor"].as<int>();
	} else {
		cerr << "Default read Length not given, exiting" << endl;
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

	unique_ptr<ChrConverter> chrConverter;
	if (!options.count("assemblyname") ||
	      options["assemblyname"].as<string>() == Hg37ChrConverter::assembly_name) {
	    chrConverter = unique_ptr<ChrConverter>(new Hg37ChrConverter());
    } else if (options["assemblyname"].as<string>() == Hg38ChrConverter::assembly_name) {
        chrConverter = unique_ptr<ChrConverter>(new Hg38ChrConverter());
    } else {
        cerr << "Unknown assembly name " << options["assemblyname"].as<string>() << ". I know "
             << Hg37ChrConverter::assembly_name << " and "
             << Hg38ChrConverter::assembly_name << endl;
        return 1;
    }

    // Initialize global application config.
    const GlobalAppConfig &config = GlobalAppConfig::init(move(chrConverter));

	MrefEntryAnno::PIDSINMREF = pidsInMref;
	unique_ptr<ifstream> mrefInputHandle {
	    make_unique<ifstream>(options["mref"].as<string>(), ios_base::in | ios_base::binary)
	};
	unique_ptr<boost::iostreams::filtering_istream> mrefGzHandle {
	    make_unique<boost::iostreams::filtering_istream>()
	};
	mrefGzHandle->push(boost::iostreams::gzip_decompressor());
	mrefGzHandle->push(*mrefInputHandle);
	cerr << "m\n";
	string line { };

    const ChrConverter &chrConverterTmp = config.getChrConverter();
	while (error_terminating_getline(*mrefGzHandle, line)) {
		if (line.front() == '#') {
			continue;
		};
		auto chrIndex =
		    chrConverterTmp.indexConverter[chrConverterTmp.readChromosomeIndex(line.cbegin(), '\t')];
		if (chrIndex < 0) {
			continue;
		}
		mref[chrIndex].emplace_back(line);
	}
	SvEvent::ARTIFACTFREQLOWTHRESHOLD = (artifactlofreq + 0.0) / 100;
	SvEvent::ARTIFACTFREQHIGHTHRESHOLD = (artifacthifreq + 0.0) / 100;
	BreakpointReduced::ARTIFACTFREQHIGHTHRESHOLD = SvEvent::ARTIFACTFREQHIGHTHRESHOLD;
	SvEvent::CLONALITYLOWTHRESHOLD = (clonalitylofreq + 0.0) / 100;
	SvEvent::CLONALITYSTRICTLOWTHRESHOLD = (clonalitystrictlofreq + 0.0) / 100;
	BreakpointReduced::CLONALITYSTRICTLOWTHRESHOLD = SvEvent::CLONALITYSTRICTLOWTHRESHOLD;
	SvEvent::CLONALITYHIGHTHRESHOLD = (clonalityhifreq + 0.0) / 100;
	SvEvent::BPFREQTHRESHOLD = pidsInMref * (bpFreq + 0.0) / 100;
	SvEvent::RELAXEDBPFREQTHRESHOLD = 3 * SvEvent::BPFREQTHRESHOLD;
	SvEvent::PIDSINMREFSTR = strtk::type_to_string<int>(pidsInMref);
	BreakpointReduced::PIDSINMREFSTR = SvEvent::PIDSINMREFSTR;
	BreakpointReduced::DEFAULTREADLENGTH = defaultReadLengthTumor;
	Breakpoint::DEFAULTREADLENGTH = defaultReadLengthTumor;
	SuppAlignment::DEFAULTREADLENGTH = defaultReadLengthTumor;
	SuppAlignmentAnno::DEFAULTREADLENGTH = defaultReadLengthTumor;
	SvEvent::HALFDEFAULTREADLENGTH = round(defaultReadLengthTumor / 2.0);
	SvEvent::GERMLINEOFFSETTHRESHOLD = germlineOffset;
	SvEvent::GERMLINEDBLIMIT = germlineDbLimit;
	SvEvent::ABRIDGEDOUTPUT = true;
	if (options.count("debugmode")) {
		SvEvent::DEBUGMODE = true;
	} else {
		SvEvent::DEBUGMODE = false;
	}
	AnnotationProcessor::ABRIDGEDOUTPUT = true;
	Breakpoint::BPSUPPORTTHRESHOLD = 3;
	if (options.count("controlresults")) {
		string controlResults { options["controlresults"].as<string>() };
		int defaultReadLengthControl { 0 };
		if (options.count("defaultreadlengthcontrol")) {
			defaultReadLengthControl = options["defaultreadlengthtumor"].as<int>();
		} else {
			cerr << "Default read Length not given, exiting" << endl;
			return 1;
		}
		auto lowQualControl = 0;
		auto pathogenInControl = false;
		{
			SvEvent::NOCONTROLMODE = true;
			AnnotationProcessor annotationProcessorControlCheck {
			    controlResults, mref, defaultReadLengthControl, true, germlineDbLimit
			};
			lowQualControl = annotationProcessorControlCheck.getMassiveInvFilteringLevel();
			pathogenInControl = annotationProcessorControlCheck.isContaminationObserved();
			SvEvent::NOCONTROLMODE = false;
		}
		AnnotationProcessor annotationProcessor {
		    tumorResults,
		    mref,
		    controlResults,
		    defaultReadLengthTumor,
		    defaultReadLengthControl,
		    germlineDbLimit,
		    lowQualControl,
		    pathogenInControl
		};
		annotationProcessor.printFilteredResults(pathogenInControl, lowQualControl);
	} else {
		SvEvent::NOCONTROLMODE = true;
		AnnotationProcessor annotationProcessor { tumorResults, mref, defaultReadLengthTumor, false, germlineDbLimit };
		annotationProcessor.printFilteredResults(false, 0);
	}

	return 0;
}
