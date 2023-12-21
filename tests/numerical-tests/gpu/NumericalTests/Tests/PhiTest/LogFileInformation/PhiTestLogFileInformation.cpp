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
#include "PhiTestLogFileInformation.h"

#include "Tests/PhiTest/PhiTest.h"
#include "Tests/PhiTest/PhiTestParameterStruct.h"

#include <iomanip>
#include <sstream>


std::shared_ptr<PhiTestLogFileInformation> PhiTestLogFileInformation::getNewInstance(std::shared_ptr<PhiTestParameterStruct> testPara)
{
    return std::shared_ptr<PhiTestLogFileInformation>(new PhiTestLogFileInformation(testPara));
}

std::string PhiTestLogFileInformation::getOutput()
{
    std::ostringstream headName;
    headName <<" Phi Test";
    makeCenterHead(headName.str());

    oss << "StartTimeStepCalculation_PhiTest=" << startTimeStepCalculation << std::endl;
    oss << "EndTimeStepCalculation_PhiTest=" << endTimeStepCalculation << std::endl;
    oss << "DataToCalc_PhiTest=\"";
    for (int i = 0; i < testGroups.size(); i++) {
        oss << testGroups.at(i).at(0)->getDataToCalculate();
        if (i < testGroups.size() - 1)
            oss << " ";
        else
            oss << "\"" << std::endl;
    }
    oss << std::endl;

    std::ostringstream failMessagePhi;
    failMessagePhi << "FailTests_Phi_PhiTest=\"";
    std::ostringstream failMessageOOA;
    failMessageOOA << "FailTests_OOA_PhiTest=\"";
    for (int i = 0; i < testGroups.size(); i++) {
        fillMyData(testGroups.at(i));
        for (int j = 0; j < lxForErase.size(); j++) {
            if (status.at(j) == passed || status.at(j) == failed) {
                oss << "PhiDiff_" << lxForErase.at(j) << "_" << dataToCalc.at(j) << "=" << phiDiff.at(j) << std::endl;
            }
            else 
                failMessagePhi << lxForErase.at(j) << "_" << dataToCalc.at(j) << " ";
        }
        oss << std::endl;
        for (int j = 0; j < orderOfAccuracy.size(); j++) {
            if (status.at(j) == passed || status.at(j) == failed) {
                oss << "OrderOfAccuracy_PhiDiff_" << lx.at(2 * j) << "_" << lx.at(2 * j + 1) << "_" << dataToCalc.at(j) << "=" << orderOfAccuracy.at(j) << std::endl;
            }
            else
                failMessageOOA << lx.at(2 * j) << "_" << lx.at(2 * j + 1) << "_" << dataToCalc.at(j) << " ";
        }
        oss << std::endl;
            
    }
    std::string failPhi = failMessagePhi.str();
    if (failPhi.back() == ' ')
        failPhi = failPhi.substr(0, failPhi.size() - 1);
    failMessagePhi.str(std::string());
    failMessagePhi << failPhi << "\"";
    oss << failMessagePhi.str() << std::endl << std::endl;

    std::string failOOA = failMessageOOA.str();
    if (failOOA.back() == ' ')
        failOOA = failOOA.substr(0, failOOA.size() - 1);
    failMessageOOA.str(std::string());
    failMessageOOA << failOOA << "\"";
    oss << failMessageOOA.str() << std::endl << std::endl;

    return oss.str();
}

void PhiTestLogFileInformation::addTestGroup(std::vector<std::shared_ptr<PhiTest> > tests)
{
    testGroups.push_back(tests);
}

void PhiTestLogFileInformation::fillMyData(std::vector<std::shared_ptr<PhiTest> > testGroup)
{
    lxForErase.resize(0);
    lx.resize(0);
    phiDiff.resize(0);
    orderOfAccuracy.resize(0);
    dataToCalc.resize(0);
    status.resize(0);
    for (int i = 0; i < testGroup.size(); i++) {
        status.push_back(testGroup.at(i)->getTestStatus());
        status.push_back(testGroup.at(i)->getTestStatus());

        std::vector<int> myLx = testGroup.at(i)->getLx();
        std::vector<double> myPhiDiff;
        if (testGroup.at(i)->getTestStatus() == simulationCrashed || testGroup.at(i)->getTestStatus() == test_error) {
            for (int i = 0; i < myLx.size(); i++)
                myPhiDiff.push_back((double)0.0);
        }
        else
            myPhiDiff = testGroup.at(i)->getPhiDiff();

        lx.insert(lx.end(), myLx.begin(), myLx.end());
        lxForErase.insert(lxForErase.end(), myLx.begin(), myLx.end());
        phiDiff.insert(phiDiff.end(), myPhiDiff.begin(), myPhiDiff.end());
        orderOfAccuracy.push_back(testGroup.at(i)->getOrderOfAccuracy());
        dataToCalc.push_back(testGroup.at(i)->getDataToCalculate());
        dataToCalc.push_back(testGroup.at(i)->getDataToCalculate());
        
    }

    for (int i = 0; i < lxForErase.size(); i++) 
        for (int j = i + 1; j < lxForErase.size(); j++)
            if (lxForErase.at(i) == lxForErase.at(j))
                lxForErase.at(j) = -1;
    
    for (int i = lxForErase.size() - 1; i >= 0; i--) {
        if (lxForErase.at(i) == -1) {
            phiDiff.erase(phiDiff.begin() + i);
            lxForErase.erase(lxForErase.begin() + i);
        }
    }

    
}

PhiTestLogFileInformation::PhiTestLogFileInformation(std::shared_ptr<PhiTestParameterStruct> testPara)
{
    startTimeStepCalculation = testPara->startTimeStepCalculation;
    endTimeStepCalculation = testPara->endTimeStepCalculation;
}
//! \}
