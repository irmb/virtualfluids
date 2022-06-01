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
//! \file LoggerImp.cpp
//! \ingroup Logger
//! \author Stephan Lenz
//=======================================================================================
#include "LoggerImp.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <sstream>

logging::LoggerImp::LoggerImp(std::ostream *stream) : logging::Logger(stream)
{
    levelString[Level::WARNING]           = "[WARNING]          ";
    levelString[Level::LOGGER_ERROR]      = "[ERROR]            ";
    levelString[Level::INFO_LOW]          = "[INFO_LOW]         ";
    levelString[Level::INFO_INTERMEDIATE] = "[INFO_INTERMEDIATE]";
    levelString[Level::INFO_HIGH]         = "[INFO_HIGH]        ";
}

logging::LoggerImp::~LoggerImp() = default;

logging::Logger &logging::LoggerImp::operator<<(const Level &level)
{
    localLogLevel = level;
    return *this;
}

logging::Logger &logging::LoggerImp::operator<<(const std::string &message) { return this->log(message); }

logging::Logger &logging::LoggerImp::operator<<(const int &message) { return this->log(std::to_string(message)); }

logging::Logger &logging::LoggerImp::operator<<(const unsigned int &message)
{
    return this->log(std::to_string(message));
}

logging::Logger &logging::LoggerImp::operator<<(const unsigned long &message)
{
    return this->log(std::to_string(message));
}

logging::Logger &logging::LoggerImp::operator<<(const float &message) { return this->log(std::to_string(message)); }

logging::Logger &logging::LoggerImp::operator<<(const double &message) { return this->log(std::to_string(message)); }

logging::Logger &logging::LoggerImp::log(const std::string &message)
{
    if (shouldBeLogged()) {
        std::string modifiedMessage = message;
        addDebugInformation(modifiedMessage);
        for (auto stream : streams)
            *stream << modifiedMessage << std::flush;
    }
    std::size_t found = message.find(std::string("\n"));
    if (found != std::string::npos)
        newLoggingLine = true;
    else
        newLoggingLine = false;

    return *this;
}

bool logging::LoggerImp::shouldBeLogged() { return localLogLevel <= globalLogLevel; }

void logging::LoggerImp::addDebugInformation(std::string &message)
{
    if (newLoggingLine) {
        std::stringstream os;
        os << levelString[localLogLevel] << getTimeStamp() << " " << message;
        message = os.str();
    }
}

std::string logging::LoggerImp::getTimeStamp()
{
    if (!timeStampEnabled)
        return "";

    const auto now = std::chrono::system_clock::now();
    time_t tt      = std::chrono::system_clock::to_time_t(now);
    // const tm utc_tm = *gmtime(&tt);
    const tm local_tm = *localtime(&tt);

    std::stringstream os;
    os << " [" << std::setw(2) << std::setfill('0') << local_tm.tm_hour << ":";
    os << std::setw(2) << std::setfill('0') << local_tm.tm_min << ":";
    os << std::setw(2) << std::setfill('0') << local_tm.tm_sec << "]";
    return os.str();
}

std::string logging::LoggerImp::getRankString()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return printRankNumber ? "[" + std::to_string(rank) + "] " : "";
}
