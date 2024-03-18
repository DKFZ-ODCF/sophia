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
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/exception/all.hpp>
#include "BreakpointReduced.h"
#include "AnnotationProcessor.h"
#include "SuppAlignment.h"
#include "SuppAlignmentAnno.h"
#include "SvEvent.h"
#include "strtk-wrap.h"
#include "MrefEntryAnno.h"
#include "HelperFunctions.h"


int main(int argc, char** argv) {
    using namespace sophia;

    std::string assemblyName = "classic_hg37";
    int artifactlofreq { 33 };
    int artifacthifreq { 50 };
    int clonalitylofreq { 5 };
    int clonalitystrictlofreq { 20 };
    int clonalityhifreq { 85 };
    int bpFreq { 3 };
    int germlineOffset { 5 };
    int germlineDbLimit { 5 };
    int pidsInMref { 0 };
    ChrSize defaultReadLengthTumor { 0 };
    ChrSize defaultReadLengthControl { 0 };

    try {
        std::ios_base::sync_with_stdio(false);
        namespace po = boost::program_options;
        std::cin.tie(nullptr);
        po::options_description options("Allowed options for sophiaAnnotate");
        options.add_options()
            ("help",
                "produce help message")
            ("mref",
                po::value<std::string>(),
                "path to mref file")
            ("tumorresults",
                po::value<std::string>(),
                "path to _bps.bed.gz file from `sophia` for the tumor, or control for a no-tumor analysis")
            ("controlresults",
                po::value<std::string>(),
                "path to _bps.bed.gz file from `sophia` for the control")
            ("assemblyname",
                po::value<std::string>(&assemblyName)->default_value(assemblyName),
                "assembly name (classic_hg37, hg38, ...)")
            ("defaultreadlengthtumor",
                po::value<ChrSize>(&defaultReadLengthTumor),
                "Default read length for the technology used in sequencing 101,151 etc., tumor")
            ("defaultreadlengthcontrol",
                po::value<ChrSize>(&defaultReadLengthControl),
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
            ("DEBUG_MODE",
                "DEBUG_MODE")
        ;

        po::variables_map inputVariables { };
        po::store(po::parse_command_line(argc, argv, options), inputVariables);
        po::notify(inputVariables);

        if (inputVariables.count("help")) {
            std::cout << options << std::endl;
            return 0;
        }

        std::optional<std::string> assemblyNameOpt { };
        if (inputVariables.count("assemblyname")) {
            assemblyNameOpt = inputVariables["assemblyname"].as<std::string>();
        }
        setApplicationConfig(assemblyNameOpt);

        CompressedMrefIndex vectorSize =
            GlobalAppConfig::getInstance().getChrConverter().nChromosomesCompressedMref();

        std::vector<std::vector<MrefEntryAnno>> mref { static_cast<unsigned int>(vectorSize), std::vector<MrefEntryAnno> { } };
        if (!inputVariables.count("mref")) {
            std::cerr << "No mref file given, exiting" << std::endl;
            return 1;
        }

        std::string tumorResults;
        if (inputVariables.count("tumorresults")) {
            tumorResults = inputVariables["tumorresults"].as<std::string>();
        } else {
            std::cerr << "No input file given, exiting" << std::endl;
            return 1;
        }

        if (inputVariables.count("pidsinmref")) {
            pidsInMref = inputVariables["pidsinmref"].as<int>();
        } else {
            std::cerr << "number of PIDS in the MREF not given, exiting" << std::endl;
            return 1;
        }

        if (inputVariables.count("defaultreadlengthtumor")) {
            defaultReadLengthTumor = inputVariables["defaultreadlengthtumor"].as<ChrSize>();
        } else {
            std::cerr << "Default read length not given, exiting" << std::endl;
            return 1;
        }
        if (defaultReadLengthTumor < 1) {
            std::cerr << "Default read length tumor " << std::to_string(defaultReadLengthTumor)
                 << " is invalid." << std::endl;
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

        MrefEntryAnno::PIDS_IN_MREF = pidsInMref;
        std::unique_ptr<std::ifstream> mrefInputHandle
            { std::make_unique<std::ifstream>(inputVariables["mref"].as<std::string>(),
                                              std::ios_base::in | std::ios_base::binary) };
        std::unique_ptr<boost::iostreams::filtering_istream> mrefGzHandle
            { std::make_unique<boost::iostreams::filtering_istream>() };

        mrefGzHandle->push(boost::iostreams::gzip_decompressor());
        mrefGzHandle->push(*mrefInputHandle);
        std::string line { };

        const ChrConverter &chrConverter = GlobalAppConfig::getInstance().getChrConverter();

        while (error_terminating_getline(*mrefGzHandle, line)) {
            if (line.front() == '#') {
                continue;
            };
            ChrIndex globalIndex;
            try {
                globalIndex = chrConverter.parseChrAndReturnIndex(line.cbegin(), line.cend(), '\t');
            } catch (const DomainError &e) {
                e <<
                    error_info_string("line = " + line);
                throw e;
            }
            CompressedMrefIndex chrIndex;
            if (!chrConverter.isCompressedMref(globalIndex)) {
                continue;
            } else {
                chrIndex = chrConverter.indexToCompressedMrefIndex(globalIndex);
                mref[static_cast<unsigned int>(chrIndex)].emplace_back(line);
            }
        }
        SvEvent::ARTIFACT_FREQ_LOW_THRESHOLD = (artifactlofreq + 0.0) / 100;
        SvEvent::ARTIFACT_FREQ_HIGH_THRESHOLD = (artifacthifreq + 0.0) / 100;
        BreakpointReduced::ARTIFACT_FREQ_HIGH_THRESHOLD = SvEvent::ARTIFACT_FREQ_HIGH_THRESHOLD;
        SvEvent::CLONALITY_LOW_THRESHOLD = (clonalitylofreq + 0.0) / 100;
        SvEvent::CLONALITY_STRICT_LOW_THRESHOLD = (clonalitystrictlofreq + 0.0) / 100;
        BreakpointReduced::CLONALITY_STRICT_LOW_THRESHOLD = SvEvent::CLONALITY_STRICT_LOW_THRESHOLD;
        SvEvent::CLONALITY_HIGH_THRESHOLD = (clonalityhifreq + 0.0) / 100;
        SvEvent::BP_FREQ_THRESHOLD = pidsInMref * (bpFreq + 0.0) / 100;
        SvEvent::RELAXED_BP_FREQ_THRESHOLD = 3 * SvEvent::BP_FREQ_THRESHOLD;
        SvEvent::PIDS_IN_MREF_STR = strtk::type_to_string<int>(pidsInMref);
        BreakpointReduced::PIDS_IN_MREF_STR = SvEvent::PIDS_IN_MREF_STR;
        BreakpointReduced::DEFAULT_READ_LENGTH = defaultReadLengthTumor;
        Breakpoint::DEFAULT_READ_LENGTH = defaultReadLengthTumor;
        SuppAlignment::DEFAULT_READ_LENGTH = defaultReadLengthTumor;
        SuppAlignmentAnno::DEFAULT_READ_LENGTH = defaultReadLengthTumor;
        SvEvent::HALF_DEFAULT_READ_LENGTH = round(defaultReadLengthTumor / 2.0);
        SvEvent::GERMLINE_OFFSET_THRESHOLD = germlineOffset;
        SvEvent::GERMLINE_DB_LIMIT = germlineDbLimit;
        SvEvent::ABRIDGED_OUTPUT = true;
        if (inputVariables.count("DEBUG_MODE")) {
            SvEvent::DEBUG_MODE = true;
        } else {
            SvEvent::DEBUG_MODE = false;
        }
        AnnotationProcessor::ABRIDGED_OUTPUT = true;
        Breakpoint::BP_SUPPORT_THRESHOLD = 3;
        if (inputVariables.count("controlresults")) {
            std::string controlResults { inputVariables["controlresults"].as<std::string>() };
            if (inputVariables.count("defaultreadlengthcontrol")) {
                defaultReadLengthControl = inputVariables["defaultreadlengthcontrol"].as<ChrSize>();
            } else {
                std::cerr << "Default read length control not given, exiting" << std::endl;
                return 1;
            }
            if (defaultReadLengthControl < 1) {
                std::cerr << "Default read length control " << std::to_string(defaultReadLengthControl)
                     << " is invalid." << std::endl;
                return 1;
            }

            auto lowQualControl = 0;
            auto pathogenInControl = false;
            {
                SvEvent::NO_CONTROL_MODE = true;
                AnnotationProcessor annotationProcessorControlCheck { controlResults, mref, defaultReadLengthControl, true, germlineDbLimit };
                lowQualControl = annotationProcessorControlCheck.getMassiveInvFilteringLevel();
                pathogenInControl = annotationProcessorControlCheck.isContaminationObserved();
                SvEvent::NO_CONTROL_MODE = false;
            }
            AnnotationProcessor annotationProcessor { tumorResults, mref, controlResults, defaultReadLengthTumor, defaultReadLengthControl, germlineDbLimit, lowQualControl, pathogenInControl };
            annotationProcessor.printFilteredResults(pathogenInControl, lowQualControl);
        } else {
            SvEvent::NO_CONTROL_MODE = true;
            AnnotationProcessor annotationProcessor { tumorResults, mref, defaultReadLengthTumor, false, germlineDbLimit };
            annotationProcessor.printFilteredResults(false, 0);
        }

        return 0;
    } catch (boost::exception &e) {
        std::cerr << "Error: " << boost::diagnostic_information(e) << std::endl;
        return 1;
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
