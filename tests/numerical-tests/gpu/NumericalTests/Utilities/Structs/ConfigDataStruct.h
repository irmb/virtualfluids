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
#ifndef CONFIG_DATA_STRUCT_H
#define CONFIG_DATA_STRUCT_H

#include <memory>
#include <string>
#include <vector>

#include "Simulations/ShearWave/ShearWaveParameterStruct.h"
#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"
#include "Tests/L2NormTest/L2NormTestParameterStruct.h"
#include "Tests/L2NormTestBetweenKernels/L2NormTestBetweenKernelsParameterStruct.h"
#include "Tests/NyTest/NyTestParameterStruct.h"
#include "Tests/PhiTest/PhiTestParameterStruct.h"
#include "Utilities/Structs/BasicSimulationParameterStruct.h"
#include "Utilities/Structs/VectorWriterInformationStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"
#include "Utilities/Structs/LogFileParameterStruct.h"

struct ConfigDataStruct
{
    std::vector<double> viscosity;
    std::vector<std::string> kernelsToTest;

    std::vector<std::shared_ptr<TaylorGreenVortexUxParameterStruct> > taylorGreenVortexUxParameter;
    std::vector<std::shared_ptr<GridInformationStruct> > taylorGreenVortexUxGridInformation;

    std::vector<std::shared_ptr<TaylorGreenVortexUzParameterStruct> > taylorGreenVortexUzParameter;
    std::vector<std::shared_ptr<GridInformationStruct> > taylorGreenVortexUzGridInformation;

    std::vector<std::shared_ptr<ShearWaveParameterStruct> > shearWaveParameter;
    std::vector<std::shared_ptr<GridInformationStruct> > shearWaveGridInformation;

    


    bool writeAnalyticalToVTK;
    unsigned int ySliceForCalculation;
    
    std::string logFilePath;

    int numberOfSimulations;

    std::shared_ptr<PhiTestParameterStruct> phiTestParameter;
    std::shared_ptr<NyTestParameterStruct> nyTestParameter;
    std::shared_ptr<L2NormTestParameterStruct> l2NormTestParameter;
    std::shared_ptr<L2NormTestBetweenKernelsParameterStruct> l2NormTestBetweenKernelsParameter;

    std::shared_ptr<VectorWriterInformationStruct> vectorWriterInfo;

    std::shared_ptr<LogFileParameterStruct> logFilePara;
};

#endif 
//! \}
