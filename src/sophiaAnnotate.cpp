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
#include "BreakpointReduced.h"
#include "AnnotationProcessor.h"
#include "SuppAlignment.h"
#include "SuppAlignmentAnno.h"
#include "SvEvent.h"
#include "strtk-wrap.h"
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <vector>
#include "MrefEntryAnno.h"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "HelperFunctions.h"
#include "Hg37ChrConverter.h"
#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"


int main(int argc, char** argv) {
    using namespace std;
    using namespace sophia;

    string assemblyName = "hg37";
    int artifactlofreq { 33 };
    int artifacthifreq { 50 };
    int clonalitylofreq { 5 };
    int clonalitystrictlofreq { 20 };
    int clonalityhifreq { 85 };
    int bpFreq { 3 };
    int germlineOffset { 5 };
    int germlineDbLimit { 5 };
    int pidsInMref { 0 };
    int defaultReadLengthTumor { 0 };
    int defaultReadLengthControl { 0 };

    try {
        ios_base::sync_with_stdio(false);
        namespace po = boost::program_options;
        cin.tie(nullptr);
        po::options_description options("Allowed options for sophiaAnnotate");
        options.add_options()
            ("help",
                "produce help message")
            ("mref",
                po::value<string>(),
                "path to mref file")
            ("tumorresults",
                po::value<string>(),
                "path to _bps.bed.gz file from `sophia` for the tumor, or control for a no-tumor analysis")
            ("controlresults",
                po::value<string>(),
                "path to _bps.bed.gz file from `sophia` for the control")
            ("assemblyname",
                po::value<string>(&assemblyName)->default_value("hg37"),
                "assembly name")
            ("defaultreadlengthtumor",
                po::value<int>(&defaultReadLengthTumor)->default_value(defaultReadLengthTumor),
                "Default read length for the technology used in sequencing 101,151 etc., tumor")
            ("defaultreadlengthcontrol",
                po::value<int>(&defaultReadLengthControl)->default_value(defaultReadLengthControl),
                "Default read length for the technology used in sequencing 101,151 etc., tumor")
            ("pidsinmref",
                po::value<int>(&pidsInMref)->default_value(pidsInMref),
                "Number of PIDs in the MREF")
            ("artifactlofreq",
                po::value<int>(&artifactlofreq)->default_value(artifactlofreq),
                "PERCENTAGE frequency of artifact supports for bps to be considered as artifact_like")
            ("artifacthifreq",
                po::value<int>(&artifacthifreq)->default_value(artifacthifreq),
                "PERCENTAGE frequency of artifact supports for bps to be considered as artifacts")
            ("clonalitylofreq",
                po::value<int>(&clonalitylofreq)->default_value(clonalitylofreq),
                "PERCENTAGE clonality for bps to be considered as extreme_subclonal")
            ("clonalitystrictlofreq",
                po::value<int>(&clonalitystrictlofreq)->default_value(clonalitystrictlofreq),
                "PERCENTAGE clonality for bps to be considered as extreme_subclonal")
            ("clonalityhifreq",
                po::value<int>(&clonalityhifreq)->default_value(clonalityhifreq),
                "PERCENTAGE clonality for bps to be considered as homozygous")
            ("bpfreq",
                po::value<int>(&bpFreq)->default_value(bpFreq),
                "PERCENTAGE frequency of a BP for consideration as rare")
            ("germlineoffset",
                po::value<int>(&germlineOffset)->default_value(germlineOffset),
                "Minimum offset a germline bp and a control bp")
            ("germlinedblimit",
                po::value<int>(&germlineDbLimit)->default_value(germlineDbLimit),
                "Maximum occurrence of germline variants in the db")
            ("debugmode",
                "debugmode")
        ;

        po::variables_map inputVariables { };
        po::store(po::parse_command_line(argc, argv, options), inputVariables);
        po::notify(inputVariables);

        if (inputVariables.count("help")) {
            cout << options << endl;
            return 0;
        }

        unique_ptr<ChrConverter> chrConverter;
        if (!inputVariables.count("assemblyname") ||
              inputVariables["assemblyname"].as<string>() == Hg37ChrConverter::assemblyName) {
            chrConverter = unique_ptr<ChrConverter>(new Hg37ChrConverter());
        } else if (inputVariables["assemblyname"].as<string>() == Hg38ChrConverter::assemblyName) {
            chrConverter = unique_ptr<ChrConverter>(new Hg38ChrConverter());
        } else {
            cerr << "Unknown assembly name " << inputVariables["assemblyname"].as<string>()
                 << ". I know " << Hg37ChrConverter::assemblyName
                 << " and "     << Hg38ChrConverter::assemblyName << endl;
            return 1;
        }


        // Initialize global application config.
        // This will raise a warning with -Wextra or -Wdangling-reference, but that's probably a
        // false positive. See https://stackoverflow.com/a/77103062/8784544.
        GlobalAppConfig::init(move(chrConverter));

        vector<vector<MrefEntryAnno>> mref
            { GlobalAppConfig::getInstance().getChrConverter().
                nChromosomesCompressedMref(), vector<MrefEntryAnno> { } };
        if (!inputVariables.count("mref")) {
            cerr << "No mref file given, exiting" << endl;
            return 1;
        }

        string tumorResults;
        if (inputVariables.count("tumorresults")) {
            tumorResults = inputVariables["tumorresults"].as<string>();
        } else {
            cerr << "No input file given, exiting" << endl;
            return 1;
        }

        if (inputVariables.count("pidsinmref")) {
            pidsInMref = inputVariables["pidsinmref"].as<int>();
        } else {
            cerr << "number of PIDS in the MREF not given, exiting" << endl;
            return 1;
        }

        if (inputVariables.count("defaultreadlengthtumor")) {
            defaultReadLengthTumor = inputVariables["defaultreadlengthtumor"].as<int>();
        } else {
            cerr << "Default read Length not given, exiting" << endl;
            return 1;
        }

        if (inputVariables.count("artifactlofreq")) {
            artifactlofreq = inputVariables["artifactlofreq"].as<int>();
        }

        if (inputVariables.count("artifacthifreq")) {
            artifacthifreq = inputVariables["artifacthifreq"].as<int>();
        }

        if (inputVariables.count("clonalitylofreq")) {
            clonalitylofreq = inputVariables["clonalitylofreq"].as<int>();
        }

        if (inputVariables.count("clonalitystrictlofreq")) {
            clonalitystrictlofreq = inputVariables["clonalitystrictlofreq"].as<int>();
        }

        if (inputVariables.count("clonalityhifreq")) {
            clonalityhifreq = inputVariables["clonalityhifreq"].as<int>();
        }

        if (inputVariables.count("bpfreq")) {
            bpFreq = inputVariables["bpfreq"].as<int>();
        }

        if (inputVariables.count("germlineoffset")) {
            germlineOffset = inputVariables["germlineoffset"].as<int>();
        }

        if (inputVariables.count("germlinedblimit")) {
            germlineDbLimit = inputVariables["germlinedblimit"].as<int>();
        }

        MrefEntryAnno::PIDSINMREF = pidsInMref;
        unique_ptr<ifstream> mrefInputHandle
            { make_unique<ifstream>(inputVariables["mref"].as<string>(), ios_base::in | ios_base::binary) };
        unique_ptr<boost::iostreams::filtering_istream> mrefGzHandle
            { make_unique<boost::iostreams::filtering_istream>() };
        mrefGzHandle->push(boost::iostreams::gzip_decompressor());
        mrefGzHandle->push(*mrefInputHandle);
        cerr << "m\n";
        string line { };

        while (error_terminating_getline(*mrefGzHandle, line)) {
            if (line.front() == '#') {
                continue;
            };
            std::optional<ChrIndex> chrIndexO = GlobalAppConfig::getInstance().getChrConverter().
                compressedMrefIndexToIndex(
                    GlobalAppConfig::getInstance().getChrConverter().parseChrAndReturnIndex(
                    line.cbegin(), line.cend(), '\t'));
            ChrIndex chrIndex;
            if (!chrIndexO.has_value()) {
                continue;
            } else {
                chrIndex = chrIndexO.value();
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
        if (inputVariables.count("debugmode")) {
            SvEvent::DEBUGMODE = true;
        } else {
            SvEvent::DEBUGMODE = false;
        }
        AnnotationProcessor::ABRIDGEDOUTPUT = true;
        Breakpoint::BPSUPPORTTHRESHOLD = 3;
        if (inputVariables.count("controlresults")) {
            string controlResults { inputVariables["controlresults"].as<string>() };
            if (inputVariables.count("defaultreadlengthcontrol")) {
                defaultReadLengthControl = inputVariables["defaultreadlengthtumor"].as<int>();
            } else {
                cerr << "Default read Length not given, exiting" << endl;
                return 1;
            }
            auto lowQualControl = 0;
            auto pathogenInControl = false;
            {
                SvEvent::NOCONTROLMODE = true;
                AnnotationProcessor annotationProcessorControlCheck { controlResults, mref, defaultReadLengthControl, true, germlineDbLimit };
                lowQualControl = annotationProcessorControlCheck.getMassiveInvFilteringLevel();
                pathogenInControl = annotationProcessorControlCheck.isContaminationObserved();
                SvEvent::NOCONTROLMODE = false;
            }
            AnnotationProcessor annotationProcessor { tumorResults, mref, controlResults, defaultReadLengthTumor, defaultReadLengthControl, germlineDbLimit, lowQualControl, pathogenInControl };
            annotationProcessor.printFilteredResults(pathogenInControl, lowQualControl);
        } else {
            SvEvent::NOCONTROLMODE = true;
            AnnotationProcessor annotationProcessor { tumorResults, mref, defaultReadLengthTumor, false, germlineDbLimit };
            annotationProcessor.printFilteredResults(false, 0);
        }

        return 0;
    } catch (exception& e) {
        cerr << "error: " << e.what() << endl;
        return 1;
    }
}
