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
#ifndef LOG_FILE_WRITER_IMP_H
#define LOG_FILE_WRITER_IMP_H

#include "LogFileWriter.h"

#include <fstream>
#include <vector>
#include <memory>

class BasicSimulationInfo;
class BasicTestLogFileInformation;
class SimulationLogFileInformation;
class LogFileHead;
class LogFileInformation;
class LogFileTimeInformation;
class TestLogFileInformation;

class LogFileWriterImp : public LogFileWriter
{
public:
    static std::shared_ptr<LogFileWriterImp> getNewInstance(std::shared_ptr<LogFileHead> logFileHead, std::shared_ptr<BasicSimulationInfo> basicSimInfo, std::shared_ptr<BasicTestLogFileInformation> basicTestInfo, std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles,
                                                        std::shared_ptr<LogFileTimeInformation> logFileTimeInfo,
                                                        std::shared_ptr<SimulationLogFileInformation> simLogInfo, 
                                                        std::string kernel, double viscosity);
    void writeLogFile(std::string basicFilePath);
    

private:
    LogFileWriterImp(std::shared_ptr<LogFileHead> logFileHead, std::shared_ptr<BasicSimulationInfo> basicSimInfo, std::shared_ptr<BasicTestLogFileInformation> basicTestInfo, std::vector<std::shared_ptr<TestLogFileInformation> > testLogFiles, std::shared_ptr<LogFileTimeInformation> logFileTimeInfo, std::shared_ptr<SimulationLogFileInformation> simLogInfo, std::string kernel, double viscosity);
    std::string calcDateAndTime();
    std::string buildFilePath(std::string basicFilePath);

    std::fstream logFile;
    std::string logFilePath;
    time_t now;
    struct tm nowLocal;
    std::string kernelName;
    double viscosity;
    std::vector<std::shared_ptr<LogFileInformation> > logFileInfo;
    std::shared_ptr<SimulationLogFileInformation> simLogInfo;
};
#endif 

//! \}
