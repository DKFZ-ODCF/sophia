/*
 * GlobalAppConfig.cpp
 *
 * Author: Philip R. Kensche Copyright (C) 2023 DKFZ Heidelberg
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
 *     LICENSE: GPL
 */

#include "GlobalAppConfig.h"
#include <mutex>
#include <stdexcept>
#include <memory>


namespace sophia {

    using namespace std;

    GlobalAppConfig* GlobalAppConfig::instance_ = nullptr;
    std::mutex GlobalAppConfig::mutex_ = std::mutex();

    const ChrConverter &GlobalAppConfig::getChrConverter() const {
        return *chrConverter;
    }

    GlobalAppConfig::GlobalAppConfig(unique_ptr<ChrConverter const> chrConverter):
        chrConverter(move(chrConverter)) {}

    GlobalAppConfig::~GlobalAppConfig() {}

    GlobalAppConfig &GlobalAppConfig::init(unique_ptr<ChrConverter const> chrConverter)
    {
        lock_guard<mutex> lock(mutex_);
        if (GlobalAppConfig::instance_ == nullptr) {
            GlobalAppConfig::instance_ = new GlobalAppConfig(move(chrConverter));
        } else {
            throw_with_trace(logic_error("GlobalAppConfig already initialized"));
        }
        return *GlobalAppConfig::instance_;
    }

    const GlobalAppConfig &GlobalAppConfig::getInstance() {
        if (GlobalAppConfig::instance_ == nullptr) {
            throw_with_trace(logic_error("GlobalAppConfig not initialized"));
        }
        return *GlobalAppConfig::instance_;
    }

}