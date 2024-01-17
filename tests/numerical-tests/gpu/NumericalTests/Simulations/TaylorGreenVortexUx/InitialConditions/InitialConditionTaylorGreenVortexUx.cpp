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
#include "InitialConditionTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

InitialConditionTaylorGreenUx::InitialConditionTaylorGreenUx(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct>  gridInfoStruct)
{
    this->amp = simParaStruct->amplitude;
    this->l0 = simParaStruct->l0;
    this->lx = gridInfoStruct->lx;
    this->lz = gridInfoStruct->lz;
    this->rho = simParaStruct->rho0;
    this->ux = simParaStruct->ux;
}

std::shared_ptr<InitialConditionTaylorGreenUx> InitialConditionTaylorGreenUx::getNewInstance(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct>  gridInfoStruct)
{
    return std::shared_ptr<InitialConditionTaylorGreenUx>(new InitialConditionTaylorGreenUx(simParaStruct, gridInfoStruct));
}

real InitialConditionTaylorGreenUx::getInitVX(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vx = (ux* l0 / lx + (amp * l0 * cos((real)2.0 * M_PI * z / lz) * sin((real)2.0 * M_PI * x / lx) / lx));
        return vx;
    }
    else
        return (real)0.0;

}

real InitialConditionTaylorGreenUx::getInitVY(int i, int level)
{
    return (real) 0.0;
}

real InitialConditionTaylorGreenUx::getInitVZ(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vz = (-amp * l0 * lz * cos((real)2.0 * M_PI * x / lx) * sin((real)2.0 * M_PI * z / lz) / (lx*lx));
        return vz;
    }
    else
        return (real) 0.0;
}

real InitialConditionTaylorGreenUx::getInitROH(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real press = (amp*pow(l0, (real)2.0)*rho*(amp*pow(lz, (real)2.0)*pow(cos(((real)2.0 * M_PI*z) / lz), (real)2.0) - (real)4.0 * (pow(lx, (real)2.0) - pow(lz, (real)2.0))*ux*cos(((real)2.0 * M_PI*z) / lz)*sin(((real)2.0 * M_PI*x) / lx) - amp*pow(lx, (real)2.0)*pow(sin(((real)2.0 * M_PI*x) / lx), (real)2.0))) / ((real)2.0*pow(lx, (real)4.0));
        real rho = (real)3.0 * press;
        return rho;
    }
    else
        return (real) 0.0;
}

real InitialConditionTaylorGreenUx::getInitPRESS(int i, int level)
{
    //nicht ben�tigt, da Druck aus Dichte berechnet wird
    return (real) 0.0;
}
//! \}
