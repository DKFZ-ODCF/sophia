/*
 * Sdust.cpp
 *
 *  Created on: 16 Apr 2016
 *      Author: Umut H. Toprak, DKFZ Heidelberg (Divisions of Theoretical
 * Bioinformatics, Bioinformatics and Omics Data Analytics and currently
 * Neuroblastoma Genomics) Copyright (C) 2018 Umut H. Toprak, Matthias
 * Schlesner, Roland Eils and DKFZ Heidelberg
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

#include "Sdust.h"
#include <iterator>

namespace sophia {

using namespace std;

Sdust::Sdust(const vector<int> &overhangIn)
    : res{}, P{}, w{}, L{0}, rW{0}, rV{0}, cW{vector<int>(WINDOWSIZE, 0)},
      cV{vector<int>(WINDOWSIZE, 0)} {
    auto wStart = 0;
    for (auto wFinish = 2; wFinish < static_cast<int>(overhangIn.size());
         ++wFinish) {
        wStart = max(wFinish - WINDOWSIZE + 1, 0);
        saveMaskedRegions(wStart);
        auto t = triplet(overhangIn, wFinish - 2);
        shiftWindow(t);
        if ((rW * 10) > (L * SCORETHRESHOLD)) {
            findPerfectRegions(wStart, rV, cV);
        }
    }
    wStart = max(0, static_cast<int>(overhangIn.size()) - WINDOWSIZE + 1);
    while (!P.empty()) {
        saveMaskedRegions(wStart);
        ++wStart;
    }
}

void
Sdust::saveMaskedRegions(int wStart) {
    if (!P.empty() && P.rbegin()->startIndex < wStart) {
        if (!res.empty()) {
            auto interval = res.back();
            if (P.rbegin()->startIndex <= (interval.endIndex + 1)) {
                res[res.size() - 1].endIndex =
                    max(P.rbegin()->endIndex, interval.endIndex);
            } else {
                res.push_back(PerfectInterval{P.rbegin()->startIndex,
                                              P.rbegin()->endIndex, 0.0});
            }
        } else {
            res.push_back(PerfectInterval{P.rbegin()->startIndex,
                                          P.rbegin()->endIndex, 0.0});
        }
        for (;;) {
            if (!P.empty() && P.rbegin()->startIndex < wStart) {
                P.erase(prev(P.end()));
            } else {
                break;
            }
        }
    }
}

void
Sdust::findPerfectRegions(int wStart, int r, vector<int> c) {
    auto maxScore = 0.0;
    for (auto i = static_cast<int>(w.size()) - L - 1; i >= 0; --i) {
        auto t = w[i];
        addTripletInfo(r, c, t);
        auto newScore = r / (static_cast<int>(w.size()) - i - 1.0);
        if ((newScore * 10) > SCORETHRESHOLD) {
            auto cit = P.cbegin();
            while (cit != P.cend()) {
                if (cit->startIndex < i + wStart) {
                    break;
                }
                maxScore = max(maxScore, cit->score);
                ++cit;
            }
            if (newScore >= maxScore) {
                P.emplace_hint(cit, PerfectInterval{i + wStart,
                                                    static_cast<int>(w.size()) +
                                                        1 + wStart,
                                                    newScore});
            }
        }
    }
}

void
Sdust::shiftWindow(int t) {
    if (w.size() >= WINDOWSIZE - 2) {
        auto s = w.front();
        w.pop_front();
        removeTripletInfo(rW, cW, s);
        if (L > static_cast<int>(w.size())) {
            --L;
            removeTripletInfo(rV, cV, s);
        }
    }
    w.push_back(t);
    ++L;
    addTripletInfo(rW, cW, t);
    addTripletInfo(rV, cV, t);
    if ((cV[t] * 10) > (SCORETHRESHOLD * 2)) {
        int s{0};
        do {
            s = w[w.size() - L];
            removeTripletInfo(rV, cV, s);
            --L;
        } while (s != t);
    }
}

void
Sdust::addTripletInfo(int &r, vector<int> &c, int t) {
    r += c[t];
    ++c[t];
}

void
Sdust::removeTripletInfo(int &r, vector<int> &c, int t) {
    --c[t];
    r -= c[t];
}

int
Sdust::triplet(const vector<int> &overhangIn, int indexPos) {
    return 16 * overhangIn[indexPos] + 4 * overhangIn[indexPos + 1] +
           overhangIn[indexPos + 2];
}

} /* namespace sophia */
