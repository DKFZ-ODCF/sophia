#include <string>
#include <fstream>
#include <iostream>
#include <utility>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include "Alignment.h"
#include "SuppAlignment.h"
#include "Breakpoint.h"
#include "SamSegmentMapper.h"
#include "HelperFunctions.h"
#include "Hg37ChrConverter.h"
#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"


std::pair<double, double> getIsizeParameters(const std::string &ISIZEFILE);
int main(int argc, char** argv) {
    using namespace sophia;

    int defaultReadLength = 0,
        clipSize = 10,
        baseQuality = 23,
        baseQualityLow = 12,
        lowQualClipSize = 5,
        isizeSigmaLevel = 5,
        bpSupport = 5;
    double properPairPercentage = 100.0;
    std::string assemblyName = "hg37";

    try {
        std::ios_base::sync_with_stdio(false);
        cin.tie(nullptr);
        namespace po = boost::program_options;
        po::options_description desc("Allowed options for sophia");
        desc.add_options()
            ("help", "printe help message")
            ("assemblyname",
                po::value<std::string>(&assemblyName)->default_value(assemblyName),
                ("assembly name (hg37, hg38)"))
            ("mergedisizes",
                po::value<std::string>(),
                "insertsize distribution file for the merged bam.")
            ("medianisize",
                po::value<double>(),
                "median insert size for the merged bam")
            ("stdisizepercentage",
                po::value<double>(),
               "percentage standard deviation of the insert size for the merged bam")
            ("defaultreadlength",
                po::value<int>(&defaultReadLength)->default_value(defaultReadLength),
                "Default read length for the technology used in sequencing 101,151 etc.")
            ("clipsize",
                po::value<int>(&clipSize)->default_value(clipSize),
                "Minimum length of soft/hard clips in the alignment.")
            ("basequality",
                po::value<int>(&baseQuality)->default_value(baseQuality),
                "Minimum median quality of split read overhangs")
            ("basequalitylow",
                po::value<int>(&baseQualityLow)->default_value(baseQualityLow),
               "If 5 consecutive bases in a split read overhang have lower quality than this strict threshold, it will be low-quality clipped")
            ("lowqualclipsize",
                po::value<int>(&lowQualClipSize)->default_value(lowQualClipSize),
                "Maximum length of a low quality split read overhang for discarding")
            ("isizesigma",
                po::value<int>(&isizeSigmaLevel)->default_value(isizeSigmaLevel),
                "The number of sds a s's mate has to be away to be called as discordant")
            ("bpsupport",
                po::value<int>(&bpSupport)->default_value(bpSupport),
                "Minimum number of reads supporting a discordant contig")
            ("properpairpercentage",
                po::value<double>(&properPairPercentage)->default_value(properPairPercentage),
                "Proper pair ratio as a percentage")
        ;
        double properPairRatio = properPairPercentage / 100.0;

        po::variables_map inputVariables { };
        po::store(po::parse_command_line(argc, argv, desc), inputVariables);
        po::notify(inputVariables);

        if (inputVariables.count("help")) {
            cout << desc << endl;
            return 0;
        }

        std::optional<std::string> assemblyNameOpt { };
        if (inputVariables.count("assemblyname")) {
            assemblyNameOpt = inputVariables["assemblyname"].as<string>();
        }
        setApplicationConfig(assemblyNameOpt);

        if (inputVariables.count("defaultreadlength")) {
            defaultReadLength = inputVariables["defaultreadlength"].as<int>();
        } else {
            cerr << "Default read Length not given, exiting" << endl;
            return 1;
        }

        if (inputVariables.count("clipsize")) {
            clipSize = inputVariables["clipsize"].as<int>();
        }

        if (inputVariables.count("basequality")) {
            baseQuality = inputVariables["basequality"].as<int>();
        }

        if (inputVariables.count("basequalitylow")) {
            baseQualityLow = inputVariables["basequalitylow"].as<int>();
        }

        if (inputVariables.count("lowqualclipsize")) {
            lowQualClipSize = inputVariables["lowqualclipsize"].as<int>();
        }

        if (inputVariables.count("isizesigma")) {
            isizeSigmaLevel = inputVariables["isizesigma"].as<int>();
        }

        if (inputVariables.count("bpsupport")) {
            bpSupport = inputVariables["bpsupport"].as<int>();
        }

        if (inputVariables.count("properpairpercentage")) {
            properPairRatio = inputVariables["properpairpercentage"].as<double>();
            properPairRatio /= 100;
            if (properPairRatio < 0.9) {
                Breakpoint::PROPERPAIRCOMPENSATIONMODE = true;
                Breakpoint::IMPROPERPAIRRATIO = 0.9 - properPairRatio;
            }
        }

        std::string mergedIsizeFile;
        if (inputVariables.count("mergedisizes")) {
            mergedIsizeFile = inputVariables["mergedisizes"].as<std::string>();
            auto isizeparams = getIsizeParameters(mergedIsizeFile);
            Alignment::ISIZEMAX = min(4000.0, isizeparams.first + isizeSigmaLevel * isizeparams.second);
            SuppAlignment::ISIZEMAX = Alignment::ISIZEMAX;
        } else {
            if (inputVariables.count("medianisize") && inputVariables.count("stdisizepercentage")) {
                auto medianIsize = inputVariables["medianisize"].as<double>();
                auto isizeStdPercentage = inputVariables["stdisizepercentage"].as<double>();
                Alignment::ISIZEMAX =
                    min(4000.0,
                        medianIsize + isizeSigmaLevel * medianIsize * isizeStdPercentage * 0.01);
                SuppAlignment::ISIZEMAX = Alignment::ISIZEMAX;
            } else {
                Alignment::ISIZEMAX = 2000.0;
                SuppAlignment::ISIZEMAX = 2000.0;
                cerr << "No insert size distribution file given, using a dummy default value of 2000 "
                     << "as the min insert size of a distant event"
                     << endl;
            }
        }

        Alignment::CLIPPEDNUCLEOTIDECOUNTTHRESHOLD = clipSize;
        Alignment::BASEQUALITYTHRESHOLD = baseQuality + 33;
        Alignment::BASEQUALITYTHRESHOLDLOW = baseQualityLow + 33;
        Alignment::LOWQUALCLIPTHRESHOLD = lowQualClipSize;
        Breakpoint::BPSUPPORTTHRESHOLD = bpSupport;
        Breakpoint::DEFAULTREADLENGTH = defaultReadLength;
        Breakpoint::DISCORDANTLOWQUALLEFTRANGE = static_cast<int>(round(defaultReadLength * 1.11));
        Breakpoint::DISCORDANTLOWQUALRIGHTRANGE = static_cast<int>(round(defaultReadLength * 0.51));

        SuppAlignment::DEFAULTREADLENGTH = defaultReadLength;
        ChosenBp::BPSUPPORTTHRESHOLD = bpSupport;
        cout << Breakpoint::COLUMNSSTR;
        SamSegmentMapper segmentRefMaster { defaultReadLength };
        segmentRefMaster.parseSamStream();

        return 0;
    } catch (boost::exception &e) {
        cerr << "Error: " << boost::diagnostic_information(e) << endl;
        return 1;
    } catch (std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}

std::pair<double, double> getIsizeParameters(const std::string &ISIZEFILE) {
	std::pair<double, double> isizeMedianStd { };
	std::ifstream infile { ISIZEFILE };
	std::string line;
	auto i = 0;
	while (sophia::error_terminating_getline(infile, line)) {
		boost::algorithm::trim_right(line);
		switch (i) {
		case 0:
			isizeMedianStd.first = boost::lexical_cast<double>(line);
			break;
		case 2:
			isizeMedianStd.second = boost::lexical_cast<double>(line);
			break;
		default:
			break;
		}
		++i;
	}
	return isizeMedianStd;
}
