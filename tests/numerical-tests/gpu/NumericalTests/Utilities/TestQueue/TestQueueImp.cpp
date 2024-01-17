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
#include "TestQueueImp.h"
#include <algorithm>

#include "Utilities/ColorConsoleOutput/ColorConsoleOutput.h"
#include "Utilities/Test/Test.h"

#include <basics/DataTypes.h>

TestSuiteResult TestQueueImp::run()
{
    for (const auto& test : tests)
        test->run();

    makeFinalOutput();

    return TestSuiteResult(std::clamp(numberOfFailedTest, 0, 1));
}

void TestQueueImp::makeFinalOutput()
{
    calcTestNumbers();
    colorOutput->makeFinalTestOutputHead(numberOfTests, numberOfExecutedTest, numberOfPassedTest, numberOfFailedTest,
                                         numberOfErrorTest, numberOfNotExecutedTest);
    for (uint i = 0; i < tests.size(); i++)
        tests.at(i)->makeConsoleOutput();
    colorOutput->makeFinalTestOutputFoot(numberOfTests, numberOfExecutedTest, numberOfPassedTest, numberOfFailedTest,
                                         numberOfErrorTest, numberOfNotExecutedTest);
}

int TestQueueImp::getNumberOfFailedTests() const noexcept
{
    return numberOfFailedTest;
}

std::shared_ptr<TestQueueImp> TestQueueImp::getNewInstance(std::shared_ptr<ColorConsoleOutput> colorOutput)
{
    return std::shared_ptr<TestQueueImp>(new TestQueueImp(colorOutput));
}

void TestQueueImp::addTest(std::shared_ptr<Test> test)
{
    tests.push_back(test);
}

TestQueueImp::TestQueueImp(std::shared_ptr<ColorConsoleOutput> colorOutput) : colorOutput(colorOutput)
{
    tests.resize(0);
}

void TestQueueImp::calcTestNumbers()
{
    numberOfTests = tests.size();
    numberOfExecutedTest = 0;
    numberOfPassedTest = 0;
    numberOfFailedTest = 0;
    numberOfErrorTest = 0;
    numberOfNotExecutedTest = 0;

    for (uint i = 0; i < tests.size(); i++) {
        switch (tests.at(i)->getTestStatus()) {
            case passed:
                numberOfPassedTest++;
                numberOfExecutedTest++;
                break;
            case failed:
                numberOfFailedTest++;
                numberOfExecutedTest++;
                break;
            case test_error:
                numberOfErrorTest++;
                break;
            case simulationCrashed:
                numberOfNotExecutedTest++;
                break;
            default:
                break;
        }
    }
}

//! \}
