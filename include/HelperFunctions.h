/*
 * HelperFunctions.h
 *
 *  Created on: 23 May 2019
 *      Author: Philip R. Kensche, DKFZ Heidelberg (Omics IT and Data Management
 * Core Facility)
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

#ifndef HELPERFUNCIONS_H_
#define HELPERFUNCIONS_H_

#include <iostream>
#include <string>

namespace sophia {

using namespace std;

const int EXITCODE_IOERROR = 1;

istream &error_terminating_getline(istream &is, string &str);

} /* namespace sophia */

#endif /* HELPERFUNCIONS_H_ */
