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
#include "BasicTestLogFileInformation.h"

std::shared_ptr<BasicTestLogFileInformation> BasicTestLogFileInformation::getNewInstance()
{
    return std::shared_ptr<BasicTestLogFileInformation>(new BasicTestLogFileInformation());;
}

std::string BasicTestLogFileInformation::getOutput()
{
    if (!outputBuild) {
        buildOutput();
        outputBuild = true;
    }
    return oss.str();
}

void BasicTestLogFileInformation::addTest(std::string testName, bool testRun)
{
    this->testName.push_back(testName);
    this->testRun.push_back(testRun);
}

BasicTestLogFileInformation::BasicTestLogFileInformation()
{
    testName.resize(0);
    testRun.resize(0);
    outputBuild = false;
}

void BasicTestLogFileInformation::buildOutput()
{
    makeCenterHead("Basic Test Information");

    for (int i = 0; i < testName.size(); i++)
        oss << testName.at(i) << "=" << std::boolalpha << testRun.at(i) << std::endl;
    oss << std::endl;

    outputBuild = true;
}

//! \}
