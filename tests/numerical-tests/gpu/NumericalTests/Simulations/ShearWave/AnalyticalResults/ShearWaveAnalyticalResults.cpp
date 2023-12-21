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
#include "ShearWaveAnalyticalResults.h"

#include "Simulations/ShearWave/ShearWaveParameterStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

std::shared_ptr<AnalyticalResults> ShearWaveAnalyticalResults::getNewInstance(double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct)
{
    return std::shared_ptr<AnalyticalResults>(new ShearWaveAnalyticalResults(viscosity, simParaStruct));
}

void ShearWaveAnalyticalResults::calc(std::shared_ptr<SimulationResults> simResults)
{
    AnalyticalResultsImp::init(simResults);

    for (int i = 0; i < numberOfTimeSteps; i++) {
        for (int j = 0; j < numberOfNodes; j++) {
            vx.at(i).at(j) = (l0*u0) / xNodes;
            vy.at(i).at(j) = (l0*v0*cos(((real)2.0 * M_PI*z.at(i).at(j)) / zNodes)*sin(((real)2.0 * M_PI*(x.at(i).at(j) + (l0*time.at(i)*u0) / xNodes)) / xNodes)) / (exp(time.at(i)*viscosity*(((real)4.0 * pow(M_PI, (real)2.0)) / pow(xNodes, (real)2.0) + ((real)4.0 * pow(M_PI, (real)2.0)) / pow(zNodes, (real)2.0)))*xNodes);
            vz.at(i).at(j) = (real)0.0;
            press.at(i).at(j) = (pow(l0, (real)2.0)*rho0*v0*sin(((real)2.0 * M_PI*z.at(i).at(j)) / zNodes)*(-(real)4.0 * exp(((real)4.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*u0*zNodes*cos(((real)2.0 * M_PI*(l0*time.at(i)*u0 + x.at(i).at(j)*xNodes)) / pow(xNodes, (real)2.0)) + v0*xNodes*pow(sin(((real)2.0 * M_PI*(l0*time.at(i)*u0 + x.at(i).at(j)*xNodes)) / pow(xNodes, (real)2.0)), (real)2.0)*sin(((real)2.0 * M_PI*z.at(i).at(j)) / zNodes))) / ((real)2.0*exp(((real)8.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*pow(xNodes, (real)3.0));
            rho.at(i).at(j) = (pow(l0, (real)2.0)*rho0*v0*sin(((real)2.0 * M_PI*z.at(i).at(j)) / zNodes)*(-(real)4.0 * exp(((real)4.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*u0*zNodes*cos(((real)2.0 * M_PI*(l0*time.at(i)*u0 + x.at(i).at(j)*xNodes)) / pow(xNodes, (real)2.0)) + v0*xNodes*pow(sin(((real)2.0 * M_PI*(l0*time.at(i)*u0 + x.at(i).at(j)*xNodes)) / pow(xNodes, (real)2.0)), (real)2.0)*sin(((real)2.0 * M_PI*z.at(i).at(j)) / zNodes))) / ((real)2.0*exp(((real)8.0 * pow(M_PI, (real)2.0)*time.at(i)*viscosity*(pow(xNodes, (real)2.0) + pow(zNodes, (real)2.0))) / (pow(xNodes, (real)2.0)*pow(zNodes, (real)2.0)))*pow(xNodes, (real)3.0));
        }
    }
    calculated = true;
}

ShearWaveAnalyticalResults::ShearWaveAnalyticalResults(double viscosity, std::shared_ptr<ShearWaveParameterStruct> simParaStruct)
    : AnalyticalResultsImp(simParaStruct->l0)
{
    this->viscosity = viscosity;
    this->u0 = simParaStruct->ux;
    this->v0 = simParaStruct->uz;
    this->l0 = simParaStruct->l0;
    this->rho0 = simParaStruct->rho0;
}
//! \}
