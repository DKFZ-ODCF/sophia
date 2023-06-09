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

#include "HelperFunctions.h"

namespace sophia {

using namespace std;

istream &
error_terminating_getline(istream &is, string &str) {
    getline(is, str);

    if (is.bad()) {
        perror("Error reading line from file");
        exit(EXITCODE_IOERROR);
    }

    return is;
}

} /* namespace sophia */
