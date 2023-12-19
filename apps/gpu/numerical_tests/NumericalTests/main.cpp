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
#include <gmock/gmock.h>
#include <mpi.h>

#include "Utilities/ConfigFileReaderNT/ConfigFileReaderNT.h"
#include "Utilities/LogFileQueue/LogFileQueue.h"
#include "Utilities/NumericalTestFactory/NumericalTestFactoryImp.h"
#include "Utilities/TestQueue/TestQueue.h"
#include "Utilities/VirtualFluidSimulationFactory/VirtualFluidSimulationFactory.h"

// validation
#include "Utilities/Calculator/FFTCalculator/FFTCalculator.h"
#include <fstream>
#include <iostream>

static TestSuiteResult startNumericalTests(const std::string &configFile, const std::string &pathNumericalTests)
{
    auto configData = vf::gpu::tests::readConfigFile(configFile, pathNumericalTests);

    std::shared_ptr<NumericalTestFactoryImp> numericalTestFactory = NumericalTestFactoryImp::getNewInstance(configData);

    std::shared_ptr<TestQueue> testQueue = numericalTestFactory->getTestQueue();
    std::shared_ptr<LogFileQueue> logFileQueue = numericalTestFactory->getLogFileQueue();

    auto result = testQueue->run();
    logFileQueue->writeLogFiles();

    return result;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    auto tests_passed = TestSuiteResult::FAILED;

    if (argc == 3) {
        tests_passed = startNumericalTests(argv[1], argv[2]);
    }
    else
        std::cout << "Configuration file must be set!: lbmgm <config file> <path to grid data>" << std::endl << std::flush;

    MPI_Finalize();

    return tests_passed;
}

//! \}
