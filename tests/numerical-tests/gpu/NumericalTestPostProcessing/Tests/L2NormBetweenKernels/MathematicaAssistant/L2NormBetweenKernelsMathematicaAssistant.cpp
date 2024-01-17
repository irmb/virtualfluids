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
#include "L2NormBetweenKernelsMathematicaAssistant.h"

#include "Tests/L2NormBetweenKernels/LogFileData/L2NormBetweenKernelsLogFileData.h"
#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"

#include <sstream> 


std::shared_ptr<L2NormBetweenKernelsMathematicaAssistant> L2NormBetweenKernelsMathematicaAssistant::getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
    return std::shared_ptr<L2NormBetweenKernelsMathematicaAssistant>(new L2NormBetweenKernelsMathematicaAssistant(functionFactory));
}

void L2NormBetweenKernelsMathematicaAssistant::makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
    std::shared_ptr<SortedDataL2NormBetweenKernels> mySortedData = sortLogFileData(logFileData);

    makeL2NormMathematicaOutput(aMathmaticaFile, mySortedData);
    makeL2NormBetweenKernelsMathematicaOutput(aMathmaticaFile, mySortedData);
}

L2NormBetweenKernelsMathematicaAssistant::L2NormBetweenKernelsMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory) : MathematicaAssistantImp(functionFactory)
{
}

bool L2NormBetweenKernelsMathematicaAssistant::checkTestParameter(std::shared_ptr<L2NormBetweenKernelsLogFileData> logFileData1, std::shared_ptr<L2NormBetweenKernelsLogFileData> logFileData2)
{
    if (logFileData1->getBasicKernel() != logFileData2->getBasicKernel())
        return false;
    if (logFileData1->getNormalizeData() != logFileData2->getNormalizeData())
        return false;
    if (logFileData1->getTimeStep() != logFileData2->getTimeStep())
        return false;
    if (logFileData1->getDataToCalculate() != logFileData2->getDataToCalculate())
        return false;


    return true;
}

std::shared_ptr<SortedDataL2NormBetweenKernels> L2NormBetweenKernelsMathematicaAssistant::sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData)
{
    std::vector<std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > > testLogFileData;
    std::vector<std::vector<std::string> > basicListNames;
    for (int i = 0; i < logFileData->getLogFileData(0)->getL2NormBetweenKernelsLogFileData().size(); i++) {
        std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > aTestLogFileDataGroup;
        aTestLogFileDataGroup.push_back(logFileData->getLogFileData(0)->getL2NormBetweenKernelsLogFileData().at(i));
        std::vector<std::string> aListNameGroup;
        aListNameGroup.push_back(logFileData->getLogFileData(0)->getSimulationSigniture());
        basicListNames.push_back(aListNameGroup);
    }
    for (int i = 0; i < logFileData->getGroupSize(); i++) {
        for (int j = 0; j < logFileData->getLogFileData(i)->getL2NormBetweenKernelsLogFileData().size(); j++) {
            std::string dataToCalc = logFileData->getLogFileData(i)->getL2NormBetweenKernelsLogFileData().at(j)->getDataToCalculate();
            bool added = false;
            for (int k = 0; k < testLogFileData.size(); k++) {
                if (checkTestParameter(logFileData->getLogFileData(i)->getL2NormBetweenKernelsLogFileData().at(j), testLogFileData.at(k).at(0))) {
                    testLogFileData.at(k).push_back(logFileData->getLogFileData(i)->getL2NormBetweenKernelsLogFileData().at(j));
                    basicListNames.at(k).push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
                    added = true;
                }
            }
            if (!added) {
                std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > aTestLogFileDataGroup;
                aTestLogFileDataGroup.push_back(logFileData->getLogFileData(i)->getL2NormBetweenKernelsLogFileData().at(j));
                testLogFileData.push_back(aTestLogFileDataGroup);
                std::vector<std::string> aListNameGroup;
                aListNameGroup.push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
                basicListNames.push_back(aListNameGroup);
            }
        }
    }
    std::shared_ptr<SortedDataL2NormBetweenKernels> mySortedData = std::shared_ptr<SortedDataL2NormBetweenKernels>(new SortedDataL2NormBetweenKernels);
    mySortedData->basicListNames = basicListNames;
    mySortedData->testLogFileData = testLogFileData;

    return mySortedData;
}

void L2NormBetweenKernelsMathematicaAssistant::makeL2NormMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2NormBetweenKernels> sortedData)
{
    for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
        std::vector<std::vector<double> > gridLengths;
        std::vector<std::vector<double> > l2Norm;
        std::vector<std::string> aListNamesList;
        std::vector<std::string> timeSteps;
        for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            l2Norm.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormForBasicKernel());
            l2Norm.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormForDivergentKernel());
            aListNamesList.push_back(sortedData->basicListNames.at(i).at(j));

            std::string divergentKernel = sortedData->basicListNames.at(i).at(j);
            int sizeOfString = sortedData->testLogFileData.at(i).at(j)->getDivergentKernel().size();
            divergentKernel.erase(divergentKernel.begin(), divergentKernel.begin() + sizeOfString);
            divergentKernel = sortedData->testLogFileData.at(i).at(j)->getBasicKernel() + divergentKernel;
            aListNamesList.push_back(divergentKernel);

            std::ostringstream timeStep;
            timeStep << "TimeStep" << sortedData->testLogFileData.at(i).at(j)->getTimeStep();
            timeSteps.push_back(timeStep.str());
            timeSteps.push_back(timeStep.str());
        }
        
        std::vector<std::string> finalListNames = finalizeListNames(aListNamesList, "L2Norm", sortedData->testLogFileData.at(i).at(0)->getDataToCalculate(), sortedData->testLogFileData.at(i).at(0)->getNormalizeData(), timeSteps);

        addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNames, gridLengths, l2Norm, "L[dx]", "L2Norm[-]");
    }
}

void L2NormBetweenKernelsMathematicaAssistant::makeL2NormBetweenKernelsMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2NormBetweenKernels> sortedData)
{
    for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
        std::vector<std::vector<double> > gridLengths;
        std::vector<std::vector<double> > l2NormBK;
        std::vector<std::string> aListNamesList;
        std::vector<std::string> timeSteps;
        std::vector<std::string> basicKernel;
        for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            l2NormBK.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormBetweenKernels());
            aListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
            basicKernel.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicKernel());
            std::ostringstream timeStep;
            timeStep << "TimeStep" << sortedData->testLogFileData.at(i).at(j)->getTimeStep();
            timeSteps.push_back(timeStep.str());
        }

        std::vector<std::string> finalListNames = finalizeListNames(aListNamesList, "L2NormBetweenKernels", sortedData->testLogFileData.at(i).at(0)->getDataToCalculate(), sortedData->testLogFileData.at(i).at(0)->getNormalizeData(), timeSteps, basicKernel);

        addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNames, gridLengths, l2NormBK, "L[dx]", "L2Norm[-]");
    }
}

L2NormBetweenKernelsMathematicaAssistant::L2NormBetweenKernelsMathematicaAssistant()
{
}

//! \}
