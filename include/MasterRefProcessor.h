/*
 * MasterRefProcessor.h
 *
 *  Created on: 27 Apr 2016
*      Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical Bioinformatics, Bioinformatics and Omics Data Analytics and currently Neuroblastoma Genomics)
 *      Copyright (C) 2018 Umut H. Toprak, Matthias Schlesner, Roland Eils and DKFZ Heidelberg
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *      LICENSE: GPL
 */

#ifndef MASTERREFPROCESSOR_H_
#define MASTERREFPROCESSOR_H_

#include <vector>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <memory>
#include <fstream>
#include <array>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <MrefEntry.h>
#include "SuppAlignment.h"
#include <BreakpointReduced.h>

namespace sophia {

    using namespace std;
    
class MasterRefProcessor {
public:
	MasterRefProcessor(const vector<string> &filesIn, const string &outputRootName, const string &version, const int defaultReadLengthIn);
	~MasterRefProcessor() = default;
private:

	unsigned long long processFile(const string &gzPath, short fileIndex);
	bool processBp(BreakpointReduced &bp, int chrIndex,  short fileIndex);
	const int NUMPIDS;
	const int DEFAULTREADLENGTH;
	unique_ptr<ofstream> mergedBpsOutput;
	vector<vector<MrefEntry>> mrefDb;
};

}
/* namespace sophiaMref */

#endif /* MASTERREFPROCESSOR_H_ */

