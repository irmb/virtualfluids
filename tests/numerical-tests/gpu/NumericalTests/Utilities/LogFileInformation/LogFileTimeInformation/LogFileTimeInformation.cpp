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
#include "LogFileTimeInformation.h"

#include "Utilities/SimulationInfo/SimulationInfo.h"

#include <iomanip>

std::shared_ptr<LogFileTimeInformation> LogFileTimeInformation::getNewInstance(std::vector<std::shared_ptr<SimulationInfo> > simInfo, bool fileWriting)
{
    return std::shared_ptr<LogFileTimeInformation>(new LogFileTimeInformation(simInfo, fileWriting));
}

std::string LogFileTimeInformation::getOutput()
{
    makeCenterHead("Simulation Time Information");
    oss << "VTKFileWriting=" << std::boolalpha << fileWriting <<std::endl << std::endl;
    for (int i = 0; i < simInfo.size(); i++) {
        oss << simInfo.at(i)->getRunTimeOutput() << std::endl;
    }
    return oss.str();
}

LogFileTimeInformation::LogFileTimeInformation(std::vector<std::shared_ptr<SimulationInfo> > simInfo, bool fileWriting) : simInfo(simInfo), fileWriting(fileWriting)
{

}
//! \}
