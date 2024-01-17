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
#include "InitialConditionShearWave.h"

#include "Simulations/ShearWave/ShearWaveParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>


InitialConditionShearWave::InitialConditionShearWave(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct)
{
    this->l0 = simParaStruct->l0;
    this->lx = gridInfoStruct->lx;
    this->lz = gridInfoStruct->lz;
    this->rho = simParaStruct->rho0;
    this->u0 = simParaStruct->ux;
    this->v0 = simParaStruct->uz;
}

std::shared_ptr<InitialConditionShearWave> InitialConditionShearWave::getNewInstance(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct)
{
    return std::shared_ptr<InitialConditionShearWave>(new InitialConditionShearWave(simParaStruct, gridInfoStruct));
}

real InitialConditionShearWave::getInitVX(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vx = l0 * u0 / lx;
        return vx;
    }
    else
        return (real)0.0;

}

real InitialConditionShearWave::getInitVY(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vy = v0 * l0 / lx * cos((real)2.0 * M_PI * z / lz) * sin((real)2.0 * M_PI * x / lx);
        return vy;
    }
    else
        return (real) 0.0;
}

real InitialConditionShearWave::getInitVZ(int i, int level)
{
    return (real) 0.0;
}

real InitialConditionShearWave::getInitROH(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real press = (l0*l0 * v0 * rho * sin(((real)2.0 * M_PI * z) / lz) * ((real)-4.0 * lz * u0 * cos(((real)2.0 * M_PI * x) / lx) + lx * v0 * sin(((real)2.0 * M_PI * x) / lx)*sin(((real)2.0 * M_PI * x) / lx) * sin(((real)2.0 * M_PI * z) / lz))) / ((real)2.0 * lx*lx*lx);
        real rho = (real)3.0 * press;
        return 0.0;
    }
    else
        return (real) 0.0;
}

real InitialConditionShearWave::getInitPRESS(int i, int level)
{
    //nicht ben�tigt, da Druck aus Dichte berechnet wird
    return (real) 0.0;
}
//! \}
