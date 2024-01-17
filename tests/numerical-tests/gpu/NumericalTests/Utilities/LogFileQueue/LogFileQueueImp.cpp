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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup NumericalTests
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#include "LogFileQueueImp.h"

#include <helper_functions.h>

#include <ctime>
#include <iomanip>

#include <basics/DataTypes.h>

#include "Utilities/LogFileWriter/LogFileWriter.h"

std::shared_ptr<LogFileQueueImp> LogFileQueueImp::getNewInstance(std::string basicLogFilePath)
{
    return std::shared_ptr<LogFileQueueImp>(new LogFileQueueImp(basicLogFilePath));
}

void LogFileQueueImp::writeLogFiles()
{
    for (uint i = 0; i < logFileWriter.size(); i++){
        logFileWriter.at(i)->writeLogFile(basicLogFilePath);
    }
}

void LogFileQueueImp::addLogFileWriter(std::shared_ptr<LogFileWriter> aLogFileWriter)
{
    logFileWriter.push_back(aLogFileWriter);
}

LogFileQueueImp::LogFileQueueImp(std::string basicLogFilePath)
{
    logFileWriter.resize(0);

    std::ostringstream oss;
    oss << basicLogFilePath << "/NumericalTestLogFiles/";
    this->basicLogFilePath = oss.str();
}

std::string LogFileQueueImp::calcDateAndTime()
{
    std::ostringstream oss;
    now = time(NULL);
    nowLocal = *localtime(&now);
    oss << std::setfill('0') << nowLocal.tm_year + 1900 << std::setw(2) << nowLocal.tm_mon + 1 << std::setw(2) << nowLocal.tm_mday << "_" << std::setw(2) << nowLocal.tm_hour << std::setw(2) << nowLocal.tm_min << std::setw(2) << nowLocal.tm_sec;
    return oss.str();
}

//! \}
