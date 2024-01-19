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
#include "AnalyticalResultsTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<AnalyticalResults> AnalyticalResultsTaylorGreenUz::getNewInstance(double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct)
{
    return std::shared_ptr<AnalyticalResults>(new AnalyticalResultsTaylorGreenUz(viscosity, simParaStruct));
}

void AnalyticalResultsTaylorGreenUz::calc(std::shared_ptr<SimulationResults> simResults)
{
    AnalyticalResultsImp::init(simResults);

    for (int i = 0; i < numberOfTimeSteps; i++) {
        for (int j = 0; j < numberOfNodes; j++) {
            vx.at(i).at(j) = (amplitude*exp( time.at(i)*viscosity*((-(real)4.0 * pow(M_PI, (real)2.0)) / pow(xNodes, (real)2.0) - ((real)4.0 * pow(M_PI, (real)2.0)) / pow(zNodes, (real)2.0)))*l0*cos(((real)2.0 * M_PI*((l0*time.at(i)*uz) / xNodes + z.at(i).at(j))) / zNodes)*sin(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)) / xNodes;
            vy.at(i).at(j) = (real)0.0;
            vz.at(i).at(j) = (l0*uz) / zNodes - (amplitude*exp(time.at(i)*viscosity*((-(real)4.0 * pow(M_PI, (real)2.0)) / pow(xNodes, (real)2.0) - ((real)4.0 * pow(M_PI, (real)2.0)) / pow(zNodes, (real)2.0)))*l0*zNodes*cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((real)2.0 * M_PI*(z.at(i).at(j) + (l0*time.at(i)*uz) / zNodes)) / zNodes)) / pow(xNodes, (real)2.0);
            press.at(i).at(j) = (amplitude*pow(l0, (real)2.0)*rho0*(amplitude*pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)*pow(cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes), (real)2.0) - (real)4.0 * exp(((real)4.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*uz*pow(xNodes, (real)2.0)*(pow(xNodes, (real)2.0) - pow(zNodes, (real)2.0))*cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)) - amplitude*pow(zNodes, (real)4.0)*pow(sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)), (real)2.0))) / ((real)2.0*exp(((real)8.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*pow(xNodes, (real)4.0)*pow(zNodes, (real)2.0));
            rho.at(i).at(j) = (amplitude*pow(l0, (real)2.0)*rho0*(amplitude*pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)*pow(cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes), (real)2.0) - (real)4.0 * exp(((real)4.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*uz*pow(xNodes, (real)2.0)*(pow(xNodes, (real)2.0) - pow(zNodes, (real)2.0))*cos(((real)2.0 * M_PI*x.at(i).at(j)) / xNodes)*sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)) - amplitude*pow(zNodes, (real)4.0)*pow(sin(((real)2.0 * M_PI*(l0*time.at(i)*uz + z.at(i).at(j)*zNodes)) / pow(zNodes, (real)2.0)), (real)2.0))) / ((real)2.0*exp(((real)8.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*pow(xNodes, (real)4.0)*pow(zNodes, (real)2.0));
        }
    }
    calculated = true;
}

AnalyticalResultsTaylorGreenUz::AnalyticalResultsTaylorGreenUz(double viscosity, std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct) : AnalyticalResultsImp(simParaStruct->l0)
{
    this->viscosity = viscosity;
    this->uz = simParaStruct->uz;
    this->amplitude = simParaStruct->amplitude;
    this->l0 = simParaStruct->l0;
    this->rho0 = simParaStruct->rho0;
}
//! \}
