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
#ifndef SIMULATION_PARAMETER_IMP_H
#define SIMULATION_PARAMETER_IMP_H

#include "SimulationParameter.h"

#include "Calculation/Calculation.h"
struct GridInformationStruct;
struct BasicSimulationParameterStruct;

class SimulationParameterImp : public SimulationParameter
{
public:
    double getViscosity();
    std::string getGridPath();
    std::string getFilePath();
    unsigned int getNumberOfGridLevels();
    unsigned int getEndTime();
    unsigned int getTimeStepLength();
    unsigned int getLx();
    unsigned int getLz();
    unsigned int getL0();
    std::vector<unsigned int> getDevices();
    double getMaxVelocity();
    
    std::shared_ptr<KernelConfiguration> getKernelConfiguration();

protected:
    SimulationParameterImp() {};
    SimulationParameterImp(std::string kernelName, double viscosity, std::shared_ptr<BasicSimulationParameterStruct> basicSimPara, std::shared_ptr<GridInformationStruct> gridInfo);

    void generateFileDirectionInMyStystem(std::string filePath);

    unsigned int timeStepLength;
    std::string filePath;
    double maxVelocity;
    real lx, l0, lz;

private:
    real viscosity;
    unsigned int numberOfTimeSteps, basisTimeStepLength;
    std::string gridPath;
    std::vector<unsigned int> devices;
    unsigned int maxLevel, numberOfGridLevels;
    std::shared_ptr<KernelConfiguration> kernelConfig;
};

#endif

//! \}
