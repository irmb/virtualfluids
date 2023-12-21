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
#include "LogFileInformationTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

std::shared_ptr<LogFileInformationTaylorGreenUx> LogFileInformationTaylorGreenUx::getNewInstance(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
    return std::shared_ptr<LogFileInformationTaylorGreenUx>(new LogFileInformationTaylorGreenUx(simParaStruct, gridInfoStruct));
}

std::string LogFileInformationTaylorGreenUx::getOutput()
{
    makeCenterHead("TaylorGreenVortex Ux Information");
    oss << "SimulationName=TaylorGreenVortexUx" << std::endl;
    oss << "Lx=\"";
    for (int i = 0; i < lx.size(); i++) {
        oss << lx.at(i);
        if (i < lx.size() - 1)
            oss << " ";
        else
            oss << "\"" << std::endl << std::endl;
    }

    for (int i = 0; i < lx.size(); i++) {
        oss << "ux_" << lx.at(i) << "=" << ux / (lx.at(i) / l0) << std::endl;
        oss << "Amplitude_"<< lx.at(i) << "=" << amplitude / (lx.at(i) / l0) << std::endl;
        oss << "l0_" << lx.at(i) << "=" << l0 << std::endl;
        oss << std::endl;
    }
    
    return oss.str();
}

std::vector<std::string> LogFileInformationTaylorGreenUx::getFilePathExtension()
{
    std::vector<std::string> myFilePath;
    myFilePath.push_back("TaylorGreenVortexUx");
    std::ostringstream oss;
    oss << "ux_" << ux << "_Amplitude_" << amplitude;
    myFilePath.push_back(oss.str());
    return myFilePath;
}

LogFileInformationTaylorGreenUx::LogFileInformationTaylorGreenUx(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::vector<std::shared_ptr<GridInformationStruct> > gridInfoStruct)
{
    this->ux = simParaStruct->ux;
    this->amplitude = simParaStruct->amplitude;
    this->l0 = simParaStruct->l0;

    for(int i = 0; i < gridInfoStruct.size(); i++)
        lx.push_back(gridInfoStruct.at(i)->lx);
}

//! \}
