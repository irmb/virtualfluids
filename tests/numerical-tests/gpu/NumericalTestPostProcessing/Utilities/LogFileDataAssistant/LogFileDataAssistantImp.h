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
//! \addtogroup NumericalTestPostProcessing
//! \ingroup numerical_tests
//! \{
//=======================================================================================
#ifndef LOGFILE_DATA_ASSISTANT_IMP_H
#define LOGFILE_DATA_ASSISTANT_IMP_H

#include "LogFileDataAssistant.h"

class LogFileDataGroupImp;
class LogFileDataAssistantStrategy;

class LogFileDataAssistantImp : public LogFileDataAssistant
{
public:
    static std::shared_ptr<LogFileDataAssistant> getNewInstance();


    std::vector<std::shared_ptr<LogFileDataGroup> > findDataCombination(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy, DataCombination combination);
    

protected:
    LogFileDataAssistantImp();
    
    std::vector<std::shared_ptr<LogFileDataGroup> > findEqualSimulationsForDifferentKernels(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy);
    std::vector<std::shared_ptr<LogFileDataGroup> > findEqualKernelSimulationsForDifferentViscosities(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy);

    std::vector<std::shared_ptr<LogFileData> > getSimulationGroupLogFileData(std::string simName, std::vector<std::shared_ptr<LogFileData> > allLogFileData);
    std::vector<std::vector<std::shared_ptr<LogFileData> > > sortLogFileDataAfterKernels(std::vector<std::shared_ptr<LogFileData> > logFileData);

    bool checkEqualSimulationData(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);
    bool checkEqualViscosity(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);
    bool checkEqualKernel(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2);

    bool checkBasicSimulationIsInLogFiles(std::vector<std::shared_ptr<LogFileData> > allLogFileData, std::string simName);

    std::vector<std::shared_ptr<LogFileDataGroup>> castToLogFileDataGroup(std::vector<std::shared_ptr<LogFileDataGroupImp>> data);

    bool equalDouble(double num1, double num2);

};
#endif
//! \}
