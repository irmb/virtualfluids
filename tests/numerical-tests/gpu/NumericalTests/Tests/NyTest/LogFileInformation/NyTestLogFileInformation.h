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
#ifndef NY_TEST_LOGFILE_INFORMATION_H
#define NY_TEST_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"
#include "Utilities/Test/TestStatus.h"

#include <memory>
#include <vector>

class NyTest;
struct NyTestParameterStruct;

class NyTestLogFileInformation : public TestLogFileInformation
{
public:
    static std::shared_ptr<NyTestLogFileInformation> getNewInstance(std::shared_ptr<NyTestParameterStruct> testPara);

    std::string getOutput();
    void addTestGroup(std::vector<std::shared_ptr<NyTest> > tests);

private:
    NyTestLogFileInformation() {};
    NyTestLogFileInformation(std::shared_ptr<NyTestParameterStruct> testPara);

    void fillMyData(std::vector<std::shared_ptr<NyTest> > testGroup);

    std::vector<std::vector<std::shared_ptr<NyTest> > > testGroups;
    unsigned int startTimeStepCalculation, endTimeStepCalculation;
    std::vector<int> lx;
    std::vector<int> lxForErase;
    std::vector<double> ny, nyDiff;
    std::vector<double> orderOfAccuracyNyDiff;
    std::vector<std::string> dataToCalc;
    std::vector<TestStatus> status;
};
#endif
//! \}
