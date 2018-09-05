/*
 * Overhang.h
 *
 *  Created on: 18 Apr 2016
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

#ifndef OVERHANGRANGE_H_
#define OVERHANGRANGE_H_

namespace sophia {
class OverhangRange {
	friend class Alignment;
public:
	OverhangRange(bool encounteredMIn, int bpPosIn, int startPosOnReadIn, int lengthIn) :
					encounteredM { encounteredMIn },
					bpPos { bpPosIn },
					startPosOnRead { startPosOnReadIn },
					length { lengthIn } {
	}
	~OverhangRange() = default;
private:
	bool encounteredM;
	int bpPos;
	int startPosOnRead;
	int length;
};
}

#endif /* OVERHANGRANGE_H_ */