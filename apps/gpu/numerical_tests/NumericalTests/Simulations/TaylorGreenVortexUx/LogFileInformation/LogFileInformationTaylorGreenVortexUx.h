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
#ifndef LOGFILE_INFORMATION_TAYLOR_GREEN_UX_H
#define LOGFILE_INFORMATION_TAYLOR_GREEN_UX_H

#include "Utilities/LogFileInformation/LogFileInformationImp.h"
#include "Utilities/LogFileInformation/SimulationLogFileInformation/SimulationLogFileInformation.h"

#include <memory>
#include <vector>

struct TaylorGreenVortexUxParameterStruct;
struct GridInformationStruct;

class LogFileInformationTaylorGreenUx : public SimulationLogFileInformation
{
public:
    static std::shared_ptr<LogFileInformationTaylorGreenUx> getNewInstance(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);
    
    std::string getOutput();
    std::vector<std::string> getFilePathExtension();
    

private:
    LogFileInformationTaylorGreenUx() {};
    LogFileInformationTaylorGreenUx(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct);

    double ux;
    double amplitude;
    std::vector<double> lx;
    int l0;
};
#endif 
//! \}
