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
#include "LogFileDataAssistantImp.h"

#include "Utilities/AlmostEquals.h"
#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroupImp.h"
#include "Utilities/LogFileDataAssistant/LogFileDataAssistantStrategy/LogFileDataAssistantStrategy.h"


#include <iostream>

std::vector<std::vector<std::shared_ptr<LogFileData>>> LogFileDataAssistantImp::sortLogFileDataAfterKernels(std::vector<std::shared_ptr<LogFileData>> logFileData)
{
    std::vector<std::string> kernelNames;
    for (int i = 0; i < logFileData.size(); i++) {
        if (i == 0)
            kernelNames.push_back(logFileData.at(i)->getKernel());
        else
        {
            bool newKernel = true;
            for (int j = 0; j < kernelNames.size(); j++) {
                if (kernelNames.at(i) == logFileData.at(i)->getKernel())
                    newKernel = false;
            }
            if (newKernel)
                kernelNames.push_back(logFileData.at(i)->getKernel());
        }
    }

    std::vector<std::vector<std::shared_ptr<LogFileData> > > logFileDataAfterKernels;
    logFileDataAfterKernels.resize(kernelNames.size());
    for (int i = 0; i < kernelNames.size(); i++) {
        for (int j = 0; j < logFileData.size(); j++) {
            if (kernelNames.at(i) == logFileData.at(j)->getKernel())
                logFileDataAfterKernels.at(i).push_back(logFileData.at(j));
        }
    }
    return logFileDataAfterKernels;
}

bool LogFileDataAssistantImp::checkEqualSimulationData(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
    if (logFileData1->getNumberOfTimeSteps() != logFileData2->getNumberOfTimeSteps())
        return false;
    if (logFileData1->getBasisTimeStepLength() != logFileData2->getBasisTimeStepLength())
        return false;

    return true;
}

bool LogFileDataAssistantImp::checkEqualViscosity(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
    if (!equalDouble(logFileData1->getViscosity(), logFileData2->getViscosity()))
        return false;

    return true;
}

bool LogFileDataAssistantImp::checkEqualKernel(std::shared_ptr<LogFileData> logFileData1, std::shared_ptr<LogFileData> logFileData2)
{
    if (logFileData1->getKernel() != logFileData2->getKernel())
        return false;

    return true;
}

bool LogFileDataAssistantImp::checkBasicSimulationIsInLogFiles(std::vector<std::shared_ptr<LogFileData>> allLogFileData, std::string simName)
{
    bool isInLogFileData = false;
    for (int i = 0; i < allLogFileData.size(); i++) {
        if (allLogFileData.at(i)->getSimName() == simName)
            return true;
    }
    return false;
}

std::vector<std::shared_ptr<LogFileDataGroup>> LogFileDataAssistantImp::castToLogFileDataGroup(std::vector<std::shared_ptr<LogFileDataGroupImp>> data)
{
    std::vector<std::shared_ptr<LogFileDataGroup>> casted;

    for (int i = 0; i < data.size(); i++)
        casted.push_back(data.at(i));
    return casted;
}

bool LogFileDataAssistantImp::equalDouble(double num1, double num2)
{
    const FloatingPoint<double> lhs(num1), rhs(num2);

    if (lhs.AlmostEquals(rhs))
        return true;
    return false;
}

std::vector<std::shared_ptr<LogFileData>> LogFileDataAssistantImp::getSimulationGroupLogFileData(std::string simName, std::vector<std::shared_ptr<LogFileData>> allLogFileData)
{
    std::vector<std::shared_ptr<LogFileData>> simGroupLogFileData;

    for (int i = 0; i < allLogFileData.size(); i++) {
        if (allLogFileData.at(i)->getSimName() == simName)
            simGroupLogFileData.push_back(allLogFileData.at(i));
    }

    return simGroupLogFileData;
}

std::shared_ptr<LogFileDataAssistant> LogFileDataAssistantImp::getNewInstance()
{
    return std::shared_ptr<LogFileDataAssistant>(new LogFileDataAssistantImp());
}

std::vector<std::shared_ptr<LogFileDataGroup>> LogFileDataAssistantImp::findDataCombination(std::vector<std::shared_ptr<LogFileData>> allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy, DataCombination combination)
{
    std::vector<std::shared_ptr<LogFileDataGroup>> myLogFileDataGroup;
    if (checkBasicSimulationIsInLogFiles(allLogFileData, strategy->getSimulationName())) {
        if (allLogFileData.size() > 1) {
            switch (combination)
            {
            case EqualSimulationsForDifferentKernels:
                myLogFileDataGroup = findEqualSimulationsForDifferentKernels(allLogFileData, strategy);
                break;
            case EqualKernelSimulationsForDifferentViscosities:
                myLogFileDataGroup = findEqualKernelSimulationsForDifferentViscosities(allLogFileData, strategy);
                break;
            default:
                break;
            }
        }
        else {
            std::shared_ptr<LogFileDataGroupImp> newGroup = LogFileDataGroupImp::getNewInstance();
            newGroup->addLogFileData(allLogFileData.at(0));
            myLogFileDataGroup.push_back(newGroup);
        }
    }
    return myLogFileDataGroup;
}

std::vector<std::shared_ptr<LogFileDataGroup>> LogFileDataAssistantImp::findEqualSimulationsForDifferentKernels(std::vector<std::shared_ptr<LogFileData>> allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy)
{
    std::vector<std::shared_ptr<LogFileData>> myLogFileData = getSimulationGroupLogFileData(strategy->getSimulationName(), allLogFileData);

    std::vector<std::shared_ptr<LogFileDataGroupImp>  > kernelGroups;
    kernelGroups.push_back(LogFileDataGroupImp::getNewInstance());
    kernelGroups.at(0)->addLogFileData(myLogFileData.at(0));

    for (int i = 1; i < myLogFileData.size(); i++) {
        bool added = false;
        for (int j = 0; j < kernelGroups.size(); j++) {
            if (checkEqualSimulationData(myLogFileData.at(i), kernelGroups.at(j)->getLogFileData(0))) {
                if (checkEqualViscosity(myLogFileData.at(i), kernelGroups.at(j)->getLogFileData(0))) {
                    if (strategy->checkSimulationParameter(myLogFileData.at(i), kernelGroups.at(j)->getLogFileData(0))) {
                        kernelGroups.at(j)->addLogFileData(myLogFileData.at(i));
                        added = true;
                    }
                }
            }
        }
        if (!added) {
            std::shared_ptr<LogFileDataGroupImp> newGroup = LogFileDataGroupImp::getNewInstance();
            newGroup->addLogFileData(myLogFileData.at(i));
            kernelGroups.push_back(newGroup);
        }
    }

    return castToLogFileDataGroup(kernelGroups);
}

std::vector<std::shared_ptr<LogFileDataGroup> > LogFileDataAssistantImp::findEqualKernelSimulationsForDifferentViscosities(std::vector<std::shared_ptr<LogFileData>> allLogFileData, std::shared_ptr<LogFileDataAssistantStrategy> strategy)
{
    std::vector<std::shared_ptr<LogFileData>> myLogFileData = getSimulationGroupLogFileData(strategy->getSimulationName(), allLogFileData);

    std::vector<std::shared_ptr<LogFileDataGroupImp>  > viscosityGroups;
    viscosityGroups.push_back(LogFileDataGroupImp::getNewInstance());
    viscosityGroups.at(0)->addLogFileData(myLogFileData.at(0));

    for (int i = 1; i < myLogFileData.size(); i++) {
        bool added = false;
        for (int j = 0; j < viscosityGroups.size(); j++) {
            if (checkEqualSimulationData(myLogFileData.at(i), viscosityGroups.at(j)->getLogFileData(0))) {
                if (checkEqualKernel(myLogFileData.at(i), viscosityGroups.at(j)->getLogFileData(0))) {
                    if (strategy->checkSimulationParameter(myLogFileData.at(i), viscosityGroups.at(j)->getLogFileData(0))) {
                        viscosityGroups.at(j)->addLogFileData(myLogFileData.at(i));
                        added = true;
                    }
                }
            }
        }
        if (!added) {
            std::shared_ptr<LogFileDataGroupImp> newGroup = LogFileDataGroupImp::getNewInstance();
            newGroup->addLogFileData(myLogFileData.at(i));
            viscosityGroups.push_back(newGroup);
        }
    }

    return castToLogFileDataGroup(viscosityGroups);
}

LogFileDataAssistantImp::LogFileDataAssistantImp()
{
    
}
//! \}
