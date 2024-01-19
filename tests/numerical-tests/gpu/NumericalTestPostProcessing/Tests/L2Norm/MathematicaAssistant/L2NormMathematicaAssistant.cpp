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
#include "L2NormMathematicaAssistant.h"

#include "Tests/L2Norm/LogFileData/L2NormLogFileData.h"
#include "Utilities/LogFileData/LogFileData.h"
#include "Utilities/LogFileData/LogFileDataGroup/LogFileDataGroup.h"

#include <sstream> 


std::shared_ptr<L2NormMathematicaAssistant> L2NormMathematicaAssistant::getNewInstance(std::shared_ptr<MathematicaFunctionFactory> functionFactory)
{
    return std::shared_ptr<L2NormMathematicaAssistant>(new L2NormMathematicaAssistant(functionFactory));
}

void L2NormMathematicaAssistant::makeMathematicaOutput(std::shared_ptr<LogFileDataGroup> logFileData, std::shared_ptr<MathematicaFile> aMathmaticaFile)
{
    std::shared_ptr<SortedDataL2Norm> mySortedData = sortLogFileData(logFileData);

    makeL2NormDiffMathematicaOutput(aMathmaticaFile, mySortedData);
    makeL2NormBasicTimeStepMathematicaOutput(aMathmaticaFile, mySortedData);
    makeL2NormDivergentTimeStepMathematicaOutput(aMathmaticaFile, mySortedData);
}

L2NormMathematicaAssistant::L2NormMathematicaAssistant(std::shared_ptr<MathematicaFunctionFactory> functionFactory) : MathematicaAssistantImp(functionFactory)
{
}

bool L2NormMathematicaAssistant::checkTestParameter(std::shared_ptr<L2NormLogFileData> logFileData1, std::shared_ptr<L2NormLogFileData> logFileData2)
{
    if (logFileData1->getDataToCalc() != logFileData2->getDataToCalc())
        return false;
    if (logFileData1->getNormalizeData() != logFileData2->getNormalizeData())
        return false;
    if (logFileData1->getBasicTimeStep() != logFileData2->getBasicTimeStep())
        return false;
    if (logFileData1->getDivergentTimeStep() != logFileData2->getDivergentTimeStep())
        return false;


    return true;
}

std::shared_ptr<SortedDataL2Norm> L2NormMathematicaAssistant::sortLogFileData(std::shared_ptr<LogFileDataGroup> logFileData)
{
    std::vector<std::vector<std::shared_ptr<L2NormLogFileData> > > testLogFileData;
    std::vector<std::vector<std::string> > basicListNames;
    for (int i = 0; i < logFileData->getLogFileData(0)->getL2NormLogFileData().size(); i++) {
        std::vector<std::shared_ptr<L2NormLogFileData> > aTestLogFileDataGroup;
        aTestLogFileDataGroup.push_back(logFileData->getLogFileData(0)->getL2NormLogFileData().at(i));
        std::vector<std::string> aListNameGroup;
        aListNameGroup.push_back(logFileData->getLogFileData(0)->getSimulationSigniture());
        basicListNames.push_back(aListNameGroup);
    }
    for (int i = 0; i < logFileData->getGroupSize(); i++) {
        for (int j = 0; j < logFileData->getLogFileData(i)->getL2NormLogFileData().size(); j++) {
            std::string dataToCalc = logFileData->getLogFileData(i)->getL2NormLogFileData().at(j)->getDataToCalc();
            bool added = false;
            for (int k = 0; k < testLogFileData.size(); k++) {
                if (checkTestParameter(logFileData->getLogFileData(i)->getL2NormLogFileData().at(j), testLogFileData.at(k).at(0))) {
                    testLogFileData.at(k).push_back(logFileData->getLogFileData(i)->getL2NormLogFileData().at(j));
                    basicListNames.at(k).push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
                    added = true;
                }
            }
            if (!added) {
                std::vector<std::shared_ptr<L2NormLogFileData> > aTestLogFileDataGroup;
                aTestLogFileDataGroup.push_back(logFileData->getLogFileData(i)->getL2NormLogFileData().at(j));
                testLogFileData.push_back(aTestLogFileDataGroup);
                std::vector<std::string> aListNameGroup;
                aListNameGroup.push_back(logFileData->getLogFileData(i)->getSimulationSigniture());
                basicListNames.push_back(aListNameGroup);
            }
        }
    }
    std::shared_ptr<SortedDataL2Norm> mySortedData = std::shared_ptr<SortedDataL2Norm>(new SortedDataL2Norm);
    mySortedData->basicListNames = basicListNames;
    mySortedData->testLogFileData = testLogFileData;

    return mySortedData;
}

void L2NormMathematicaAssistant::makeL2NormDiffMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData)
{
    for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
        std::vector<std::vector<double> > gridLengths;
        std::vector<std::vector<double> > l2NormDiff;
        std::vector<std::string> aBasicListNamesList;
        for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            l2NormDiff.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormDiff());
            aBasicListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
        }
        std::vector<std::string> finalListNames = finalizeListNames(aBasicListNamesList, "L2NormDiff", sortedData->testLogFileData.at(i).at(0)->getDataToCalc(), sortedData->testLogFileData.at(i).at(0)->getNormalizeData());
        
        addSecondOrderOfAccuracyRef(gridLengths, l2NormDiff, finalListNames);
        addFourthOrderOfAccuracyRef(gridLengths, l2NormDiff, finalListNames);
        addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNames, gridLengths, l2NormDiff, "L[dx]", "L2NormDiff[-]");
    }
}

void L2NormMathematicaAssistant::makeL2NormAllTimeStepsMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData)
{
    for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
        std::vector<std::vector<double> > gridLengths;
        std::vector<std::vector<double> > l2Norm;
        std::vector<std::string> aBasicListNamesList;
        std::vector<std::string> basicTimeSteps;
        for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            l2Norm.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormForBasicTimeStep());
            l2Norm.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormForDivergentTimeStep());
            aBasicListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
            aBasicListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
            std::ostringstream aBasicTimeStep;
            aBasicTimeStep << "TimeStep" << sortedData->testLogFileData.at(i).at(j)->getBasicTimeStep();
            basicTimeSteps.push_back(aBasicTimeStep.str());
            std::ostringstream aDivTimeStep;
            aDivTimeStep << "TimeStep" << sortedData->testLogFileData.at(i).at(j)->getDivergentTimeStep();
            basicTimeSteps.push_back(aDivTimeStep.str());
        }

        std::vector<std::string> finalListNamesBasic = finalizeListNames(aBasicListNamesList, "L2NormAllTimeSteps", sortedData->testLogFileData.at(i).at(0)->getDataToCalc(), sortedData->testLogFileData.at(i).at(0)->getNormalizeData(), basicTimeSteps);

        addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNamesBasic, gridLengths, l2Norm, "L[dx]", "L2Norm[-]");
    }
}

void L2NormMathematicaAssistant::makeL2NormBasicTimeStepMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData)
{
    for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
        std::vector<std::vector<double> > gridLengths;
        std::vector<std::vector<double> > l2NormBasic;
        std::vector<std::string> aBasicListNamesList;
        for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            l2NormBasic.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormForBasicTimeStep());
            aBasicListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
        }
        std::ostringstream basicTimeStep;
        basicTimeStep << "TimeStep" << sortedData->testLogFileData.at(i).at(0)->getBasicTimeStep();
        std::vector<std::string> finalListNamesBasic = finalizeListNames(aBasicListNamesList, "L2Norm", sortedData->testLogFileData.at(i).at(0)->getDataToCalc(), sortedData->testLogFileData.at(i).at(0)->getNormalizeData(), basicTimeStep.str());
        
        addSecondOrderOfAccuracyRef(gridLengths, l2NormBasic, finalListNamesBasic);
        addFourthOrderOfAccuracyRef(gridLengths, l2NormBasic, finalListNamesBasic);
        addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNamesBasic, gridLengths, l2NormBasic, "L[dx]", "L2Norm[-]");
    }
}

void L2NormMathematicaAssistant::makeL2NormDivergentTimeStepMathematicaOutput(std::shared_ptr<MathematicaFile> aMathmaticaFile, std::shared_ptr<SortedDataL2Norm> sortedData)
{
    for (int i = 0; i < sortedData->testLogFileData.size(); i++) {
        std::vector<std::vector<double> > gridLengths;
        std::vector<std::vector<double> > l2NormDivergent;
        std::vector<std::string> aBasicListNamesList;
        for (int j = 0; j < sortedData->testLogFileData.at(i).size(); j++) {
            gridLengths.push_back(sortedData->testLogFileData.at(i).at(j)->getBasicGridLengths());
            l2NormDivergent.push_back(sortedData->testLogFileData.at(i).at(j)->getL2NormForDivergentTimeStep());
            aBasicListNamesList.push_back(sortedData->basicListNames.at(i).at(j));
        }
        std::ostringstream divTimeStep;
        divTimeStep << "TimeStep" << sortedData->testLogFileData.at(i).at(0)->getDivergentTimeStep();
        std::vector<std::string> finalListNamesDiv = finalizeListNames(aBasicListNamesList, "L2Norm", sortedData->testLogFileData.at(i).at(0)->getDataToCalc(), sortedData->testLogFileData.at(i).at(0)->getNormalizeData(), divTimeStep.str());
        
        addSecondOrderOfAccuracyRef(gridLengths, l2NormDivergent, finalListNamesDiv);
        addFourthOrderOfAccuracyRef(gridLengths, l2NormDivergent, finalListNamesDiv);
        addListLogLogPlotToMathematicaFile(aMathmaticaFile, finalListNamesDiv, gridLengths, l2NormDivergent, "L[dx]", "L2Norm[-]");
    }
}

L2NormMathematicaAssistant::L2NormMathematicaAssistant()
{
}

//! \}
