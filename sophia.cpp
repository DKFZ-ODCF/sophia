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
#include "Alignment.h"
#include "SuppAlignment.h"
#include "Breakpoint.h"
#include "SamSegmentMapper.h"
#include "ChrConverter.h"
#include "HelperFunctions.h"

std::pair<double, double> getIsizeParameters(const std::string &ISIZEFILE);
int main(int argc, char** argv) {
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(nullptr);
	boost::program_options::options_description desc("Allowed options");
	desc.add_options() //
	("help", "produce help message") //
	("mergedisizes", boost::program_options::value<std::string>(), "insertsize distribution file for the merged bam") //
	("medianisize", boost::program_options::value<double>(), "median insert size for the merged bam") //
	("stdisizepercentage", boost::program_options::value<double>(), "percentage standard deviation of the insert size for the merged bam") //
	("defaultreadlength", boost::program_options::value<int>(), "Default read length for the technology used in sequencing 101,151 etc.") //
	("clipsize", boost::program_options::value<int>(), "Minimum length of soft/hard clips in the alignment. (10)") //
	("basequality", boost::program_options::value<int>(), "Minimum median quality of split read overhangs. (23)") //
	("basequalitylow", boost::program_options::value<int>(), "If 5 consecutive bases in a split read overhang have lower quality than this strict threshold, it will be low-qual. (12)") //
	("lowqualclipsize", boost::program_options::value<int>(), "Maximum length of a low qality split read overhang for discarding. (5)") //
	("isizesigma", boost::program_options::value<int>(), "The number of sds a s's mate has to be away to be called as discordant. (5)") //
	("bpsupport", boost::program_options::value<int>(), "Minimum number of reads supporting a discordant contig. (5)") //
	("properpairratio", boost::program_options::value<double>(), "Proper pair ratio as a percentage (100.0)");
	boost::program_options::variables_map inputVariables { };
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), inputVariables);
	boost::program_options::notify(inputVariables);
	if (inputVariables.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}
	int defaultReadLength { 0 };
	if (inputVariables.count("defaultreadlength")) {
		defaultReadLength = inputVariables["defaultreadlength"].as<int>();
	} else {
		std::cerr << "Default read Length not given, exiting" << std::endl;
		return 1;
	}
	auto clipSize = 10, baseQuality = 23, baseQualityLow = 12, lowQualClipSize = 5, isizeSigmaLevel = 5;
	auto bpSupport = 5;
	auto properPairRatio = 1.0;
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
	if (inputVariables.count("properpairratio")) {
		properPairRatio = inputVariables["properpairratio"].as<double>();
		properPairRatio /= 100;
		if (properPairRatio < 0.9) {
			sophia::Breakpoint::PROPERPAIRCOMPENSATIONMODE = true;
			sophia::Breakpoint::IMPROPERPAIRRATIO = 0.9 - properPairRatio;
		}
	}
	std::string mergedIsizeFile;
	if (inputVariables.count("mergedisizes")) {
		mergedIsizeFile = inputVariables["mergedisizes"].as<std::string>();
		auto isizeparams = getIsizeParameters(mergedIsizeFile);
		sophia::Alignment::ISIZEMAX = std::min(4000.0, isizeparams.first + isizeSigmaLevel * isizeparams.second);
		sophia::SuppAlignment::ISIZEMAX = sophia::Alignment::ISIZEMAX;
	} else {
		if (inputVariables.count("medianisize") && inputVariables.count("stdisizepercentage")) {
			auto medianIsize = inputVariables["medianisize"].as<double>();
			auto isizeStdPercentage = inputVariables["stdisizepercentage"].as<double>();
			sophia::Alignment::ISIZEMAX = std::min(4000.0, medianIsize + isizeSigmaLevel * medianIsize * isizeStdPercentage * 0.01);
			sophia::SuppAlignment::ISIZEMAX = sophia::Alignment::ISIZEMAX;
		} else {
			sophia::Alignment::ISIZEMAX = 2000.0;
			sophia::SuppAlignment::ISIZEMAX = 2000.0;
			std::cerr << "No insert size distribution file given, using a dummy default value of 2000 as the min insert size of a distant event" << std::endl;
		}
	}
	sophia::Alignment::CLIPPEDNUCLEOTIDECOUNTTHRESHOLD = clipSize;
	sophia::Alignment::BASEQUALITYTHRESHOLD = baseQuality + 33;
	sophia::Alignment::BASEQUALITYTHRESHOLDLOW = baseQualityLow + 33;
	sophia::Alignment::LOWQUALCLIPTHRESHOLD = lowQualClipSize;
	sophia::Breakpoint::BPSUPPORTTHRESHOLD = bpSupport;
	sophia::Breakpoint::DEFAULTREADLENGTH = defaultReadLength;
	sophia::Breakpoint::DISCORDANTLOWQUALLEFTRANGE = static_cast<int>(std::round(defaultReadLength * 1.11));
	sophia::Breakpoint::DISCORDANTLOWQUALRIGHTRANGE = static_cast<int>(std::round(defaultReadLength * 0.51));

	sophia::SuppAlignment::DEFAULTREADLENGTH = defaultReadLength;
	sophia::ChosenBp::BPSUPPORTTHRESHOLD = bpSupport;
	std::cout << sophia::Breakpoint::COLUMNSSTR;
	sophia::SamSegmentMapper segmentRefMaster { defaultReadLength };
	segmentRefMaster.parseSamStream();
	return 0;
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
