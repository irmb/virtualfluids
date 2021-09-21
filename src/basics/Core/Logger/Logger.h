//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Logger.h
//! \ingroup Logger
//! \author Stephan Lenz
//=======================================================================================
#ifndef Logger_H
#define Logger_H

#include "basics_export.h"

#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace logging
{
class BASICS_EXPORT Logger
{
protected:
    Logger(std::ostream *stream);

public:
    virtual ~Logger();

    enum Level { INFO_LOW = 3, INFO_INTERMEDIATE = 2, INFO_HIGH = 1, WARNING = 0, LOGGER_ERROR = -1 };

    enum TimeStamp { ENABLE, DISABLE };

    static void setStream(std::ostream *stream);
    static void addStream(std::ostream *stream);
    static void resetStreams();

    static void timeStamp(TimeStamp timeStamp);

    static void setDebugLevel(const Level &level = Level::LOGGER_ERROR);
    static void enablePrintedRankNumbers(bool printRankNumbers);

    virtual Logger &operator<<(const Level &level)       = 0;
    virtual Logger &operator<<(const std::string &log)   = 0;
    virtual Logger &operator<<(const int &log)           = 0;
    virtual Logger &operator<<(const unsigned int &log)  = 0;
    virtual Logger &operator<<(const unsigned long &log) = 0;
    virtual Logger &operator<<(const float &log)         = 0;
    virtual Logger &operator<<(const double &log)        = 0;

protected:
    void addStreamToList(std::ostream *stream);
    void resetStreamList();

    std::vector<std::ostream *> streams;

    static Level globalLogLevel;
    static Level localLogLevel;
    static bool printRankNumber;
    static bool timeStampEnabled;
};
extern BASICS_EXPORT std::shared_ptr<Logger> out;
} // namespace logging

#endif
