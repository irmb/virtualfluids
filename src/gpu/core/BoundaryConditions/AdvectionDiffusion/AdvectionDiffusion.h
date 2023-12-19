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
//! \addtogroup gpu_BoundaryConditions BoundaryConditions
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#ifndef AdvectionDiffusion_H
#define AdvectionDiffusion_H

#include "Calculation/Calculation.h"

#include <cuda.h>
#include <cuda_runtime.h>

struct LBMSimulationParameter;
class Parameter;

//////////////////////////////////////////////////////////////////////////
//! \brief defines the behavior of a slip-AD boundary condition
void AdvectionDiffusionSlipVelocityCompressible(
    uint numberOfThreads,
    real * normalX,
    real * normalY,
    real * normalZ,
    real * distributions,
    real * distributionsAD,
    int* QindexArray,
    real * Qarrays,
    uint numberOfBCnodes,
    real omegaDiffusivity,
    uint * neighborX,
    uint * neighborY,
    uint * neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep);
    
void AdvectionDiffusionDirichlet(
    unsigned int numberOfThreads,
    real* DD, 
    real* DD27,
    real* temp,
    real diffusivity,
    int* k_Q, 
    real* QQ,
    unsigned int numberOfBCnodes, 
    real om1, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep);

void AdvectionDiffusionBounceBack(
    unsigned int numberOfThreads,
    real* DD, 
    real* DD27,
    real* temp,
    real diffusivity,
    int* k_Q, 
    real* QQ,
    unsigned int numberOfBCnodes, 
    real om1, 
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes, 
    bool isEvenTimestep);

#endif

//! \}
