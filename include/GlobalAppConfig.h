/*
 * GlobalAppConfig.h
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

#ifndef GLOBALAPPCONFIG_H
#define GLOBALAPPCONFIG_H

#include <mutex>
#include <memory>
#include "ChrConverter.h"


namespace sophia {

    /** Keep global application config in this singleton. This is mostly to avoid having to hand
        around configurations. */
    class GlobalAppConfig {

      private:
        static GlobalAppConfig *instance_;
        static std::mutex mutex_;

      protected:

        GlobalAppConfig(std::unique_ptr<ChrConverter const> chrConverter);

        ~GlobalAppConfig();

        /** The chromosome converter. */
        const std::unique_ptr<ChrConverter const> chrConverter;

      public:

        const ChrConverter &getChrConverter() const;

        /** Prevent copying. */
        GlobalAppConfig(GlobalAppConfig &other) = delete;

        /** Prevent assignment. */
        void operator=(const GlobalAppConfig &) = delete;

        /** Factory method. */
        static GlobalAppConfig &init(std::unique_ptr<ChrConverter const> chrConverter);

        /** Getter. */
        static const GlobalAppConfig &getInstance();

    };

} /* namespace sophia */

#endif /* GLOBALAPPCONFIG_H */