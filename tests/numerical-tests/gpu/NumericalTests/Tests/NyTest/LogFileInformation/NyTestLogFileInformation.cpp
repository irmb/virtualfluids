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
#include "NyTestLogFileInformation.h"

#include "Tests/NyTest/NyTest.h"
#include "Tests/NyTest/NyTestParameterStruct.h"

#include <iomanip>
#include <sstream>


std::shared_ptr<NyTestLogFileInformation> NyTestLogFileInformation::getNewInstance(std::shared_ptr<NyTestParameterStruct> testPara)
{
    return std::shared_ptr<NyTestLogFileInformation>(new NyTestLogFileInformation(testPara));
}

std::string NyTestLogFileInformation::getOutput()
{
    std::ostringstream headName;
    headName <<"Ny Test";
    makeCenterHead(headName.str());

    oss << "StartTimeStepCalculation_NyTest=" << startTimeStepCalculation << std::endl;
    oss << "EndTimeStepCalculation_NyTest=" << endTimeStepCalculation << std::endl;
    oss << "DataToCalc_NyTest=\"";
    for (int i = 0; i < testGroups.size(); i++) {
        oss << testGroups.at(i).at(0)->getDataToCalculate();
        if (i < testGroups.size() - 1)
            oss << " ";
        else
            oss << "\"" << std::endl;
    }
    oss << std::endl;

    std::ostringstream failMessageNy;
    failMessageNy << "FailTests_Ny_NyTest=\"";
    std::ostringstream failMessageOOA;
    failMessageOOA << "FailTests_OOA_NyTest=\"";
    for (int i = 0; i < testGroups.size(); i++) {
        fillMyData(testGroups.at(i));
        for (int j = 0; j < lxForErase.size(); j++) {
            if (status.at(j) == passed || status.at(j) == failed) {
                oss << "Ny_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << ny.at(j) << std::endl;
                oss << "NyDiff_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << nyDiff.at(j) << std::endl;
            }
            else
                failMessageNy << lxForErase.at(j) << "_" << dataToCalc.at(j) << " ";
        }
        oss << std::endl;
        for (int j = 0; j < orderOfAccuracyNyDiff.size(); j++) {
            if (status.at(j) == passed || status.at(j) == failed) {
                oss << "OrderOfAccuracy_NyDiff_" << lx.at(2 * j) << "_" << lx.at(2 * j + 1) << "_" << dataToCalc.at(j) << "=" << orderOfAccuracyNyDiff.at(j) << std::endl;
            }
            else
                failMessageOOA << lx.at(2 * j) << "_" << lx.at(2 * j + 1) << "_" << dataToCalc.at(j) << " ";
        }
        oss << std::endl;
    }
    std::string failNy = failMessageNy.str();
    if (failNy.back() == ' ')
        failNy = failNy.substr(0, failNy.size() - 1);
    failMessageNy.str(std::string());
    failMessageNy << failNy << "\"";
    oss << failMessageNy.str() << std::endl << std::endl;

    std::string failOOA = failMessageOOA.str();
    if (failOOA.back() == ' ')
        failOOA = failOOA.substr(0, failOOA.size() - 1);
    failMessageOOA.str(std::string());
    failMessageOOA << failOOA << "\"";
    oss << failMessageOOA.str() << std::endl << std::endl;

    return oss.str();
}

void NyTestLogFileInformation::addTestGroup(std::vector<std::shared_ptr<NyTest> > tests)
{
    testGroups.push_back(tests);
}

void NyTestLogFileInformation::fillMyData(std::vector<std::shared_ptr<NyTest> > testGroup)
{
    lxForErase.resize(0);
    lx.resize(0);
    ny.resize(0);
    nyDiff.resize(0);
    orderOfAccuracyNyDiff.resize(0);
    dataToCalc.resize(0);
    status.resize(0);
    for (int i = 0; i < testGroup.size(); i++) {
        status.push_back(testGroup.at(i)->getTestStatus());
        status.push_back(testGroup.at(i)->getTestStatus());

        std::vector<int> myLx = testGroup.at(i)->getLx();
        std::vector<double> myNy;
        std::vector<double> myNyDiff;

        if (testGroup.at(i)->getTestStatus() == simulationCrashed || testGroup.at(i)->getTestStatus() == test_error) {
            for (int i = 0; i < myLx.size(); i++) {
                myNy.push_back((double)0.0);
                myNyDiff.push_back((double)0.0);
            }
                
        }
        else {
            myNy = testGroup.at(i)->getNy();
            myNyDiff = testGroup.at(i)->getNyDiff();
        }
        lx.insert(lx.end(), myLx.begin(), myLx.end());
        lxForErase.insert(lxForErase.end(), myLx.begin(), myLx.end());
        ny.insert(ny.end(), myNy.begin(), myNy.end());
        nyDiff.insert(nyDiff.end(), myNyDiff.begin(), myNyDiff.end());
        orderOfAccuracyNyDiff.push_back(testGroup.at(i)->getOrderOfAccuracyNyDiff());
        dataToCalc.push_back(testGroup.at(i)->getDataToCalculate());
        dataToCalc.push_back(testGroup.at(i)->getDataToCalculate());
    }

    for (int i = 0; i < lxForErase.size(); i++) 
        for (int j = i + 1; j < lxForErase.size(); j++)
            if (lxForErase.at(i) == lxForErase.at(j))
                lxForErase.at(j) = -1;
    
    for (int i = lxForErase.size() - 1; i >= 0; i--) {
        if (lxForErase.at(i) == -1) {
            ny.erase(ny.begin() + i);
            nyDiff.erase(nyDiff.begin() + i);
            lxForErase.erase(lxForErase.begin() + i);
        }
    }

    
}

NyTestLogFileInformation::NyTestLogFileInformation(std::shared_ptr<NyTestParameterStruct> testPara)
{
    startTimeStepCalculation = testPara->startTimeStepCalculation;
    endTimeStepCalculation = testPara->endTimeStepCalculation;
}
//! \}
