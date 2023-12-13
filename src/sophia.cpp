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
#include "ChrConverter.h"
#include "Hg37ChrConverter.h"
#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"


using namespace std;
using namespace sophia;

pair<double, double> getIsizeParameters(const string &ISIZEFILE);
int main(int argc, char** argv) {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);
	boost::program_options::options_description desc("Allowed options");
	desc.add_options() //
	("help", "produce help message") //
    ("assemblyname", boost::program_options::value<string>(), "assembly name (hg37, hg38)") //
	("mergedisizes", boost::program_options::value<string>(), "insertsize distribution file for the merged bam") //
	("medianisize", boost::program_options::value<double>(), "median insert size for the merged bam") //
	("stdisizepercentage", boost::program_options::value<double>(), "percentage standard deviation of the insert size for the merged bam") //
	("defaultreadlength", boost::program_options::value<int>(), "Default read length for the technology used in sequencing 101,151 etc.") //
	("clipsize", boost::program_options::value<int>(), "Minimum length of soft/hard clips in the alignment. (10)") //
	("basequality", boost::program_options::value<int>(), "Minimum median quality of split read overhangs. (23)") //
	("basequalitylow", boost::program_options::value<int>(), "If 5 consecutive bases in a split read overhang have lower quality than this strict threshold, it will be low-qual. (12)") //
	("lowqualclipsize", boost::program_options::value<int>(), "Maximum length of a low qality split read overhang for discarding. (5)") //
	("isizesigma", boost::program_options::value<int>(), "The number of sds a s's mate has to be away to be called as discordant. (5)") //
	("bpsupport", boost::program_options::value<int>(), "Minimum number of reads supporting a discordant contig. (5)") //
	("properpairpercentage", boost::program_options::value<double>(), "Proper pair ratio as a percentage (100.0)");

	boost::program_options::variables_map inputVariables { };
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), inputVariables);
	boost::program_options::notify(inputVariables);

	if (inputVariables.count("help")) {
		cout << desc << endl;
		return 0;
	}

	int defaultReadLength { 0 };
	if (inputVariables.count("defaultreadlength")) {
		defaultReadLength = inputVariables["defaultreadlength"].as<int>();
	} else {
		cerr << "Default read Length not given, exiting" << endl;
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

	if (inputVariables.count("properpairpercentage")) {
		properPairRatio = inputVariables["properpairpercentage"].as<double>();
		properPairRatio /= 100;
		if (properPairRatio < 0.9) {
			Breakpoint::PROPERPAIRCOMPENSATIONMODE = true;
			Breakpoint::IMPROPERPAIRRATIO = 0.9 - properPairRatio;
		}
	}

	string mergedIsizeFile;
	if (inputVariables.count("mergedisizes")) {
		mergedIsizeFile = inputVariables["mergedisizes"].as<string>();
		auto isizeparams = getIsizeParameters(mergedIsizeFile);
		Alignment::ISIZEMAX = min(4000.0, isizeparams.first + isizeSigmaLevel * isizeparams.second);
		SuppAlignment::ISIZEMAX = Alignment::ISIZEMAX;
	} else {
		if (inputVariables.count("medianisize") && inputVariables.count("stdisizepercentage")) {
			auto medianIsize = inputVariables["medianisize"].as<double>();
			auto isizeStdPercentage = inputVariables["stdisizepercentage"].as<double>();
			Alignment::ISIZEMAX = min(4000.0, medianIsize + isizeSigmaLevel * medianIsize * isizeStdPercentage * 0.01);
			SuppAlignment::ISIZEMAX = Alignment::ISIZEMAX;
		} else {
			Alignment::ISIZEMAX = 2000.0;
			SuppAlignment::ISIZEMAX = 2000.0;
			cerr << "No insert size distribution file given, using a dummy default value of 2000 as the min insert size of a distant event" << endl;
		}
	}

    unique_ptr<ChrConverter> chrConverter;
	if (!inputVariables.count("assemblyname") ||
	      inputVariables["assemblyname"].as<string>() == Hg37ChrConverter::assembly_name) {
	    chrConverter = unique_ptr<ChrConverter>(new Hg37ChrConverter());
    } else if (inputVariables["assemblyname"].as<string>() == Hg38ChrConverter::assembly_name) {
        chrConverter = unique_ptr<ChrConverter>(new Hg38ChrConverter());
    } else {
        cerr << "Unknown assembly name " << inputVariables["assemblyname"].as<string>() << ". I know "
             << Hg37ChrConverter::assembly_name << " and "
             << Hg38ChrConverter::assembly_name << endl;
        return 1;
    }

    // Initialize the global application configuration.
    GlobalAppConfig::init(move(chrConverter));

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
}

pair<double, double> getIsizeParameters(const string &ISIZEFILE) {
	pair<double, double> isizeMedianStd { };
	ifstream infile { ISIZEFILE };
	string line;
	auto i = 0;
	while (error_terminating_getline(infile, line)) {
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
