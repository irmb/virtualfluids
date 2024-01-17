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
#ifndef L2NORM_LOGFILE_INFORMATION_H
#define L2NORM_LOGFILE_INFORMATION_H

#include "Utilities/LogFileInformation/TestLogFileInformation/TestLogFileInformation.h"

#include <memory>
#include <vector>

class L2NormTest;
struct L2NormTestParameterStruct;

class L2NormInformation : public TestLogFileInformation
{
public:
    static std::shared_ptr<L2NormInformation> getNewInstance(std::vector<std::shared_ptr<L2NormTest> > tests, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::string> dataToCalcTests);

    std::string getOutput();

private:
    L2NormInformation() {};
    L2NormInformation(std::vector<std::shared_ptr<L2NormTest> > tests, std::shared_ptr<L2NormTestParameterStruct> testParameter, std::vector<std::string> dataToCalcTests);

    std::vector<std::shared_ptr<L2NormTest> > tests;

    unsigned int basicTimeStep, divergentTimeStep;
    std::vector<std::string> dataToCalc;
    std::vector<std::string> normalizeData;
};
#endif
//! \}
