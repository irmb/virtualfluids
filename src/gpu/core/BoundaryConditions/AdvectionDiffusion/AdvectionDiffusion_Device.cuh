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
#ifndef AdvectionDiffusion_Device_H
#define AdvectionDiffusion_Device_H

#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Calculation/Calculation.h"

//////////////////////////////////////////////////////////////////////////

template <BoundaryConditionFactory::AdvectionDiffusionSlipVelocityBC bcType>
__global__ void
AdvectionDiffusionSlipVelocity_Device(real* distributions, AdvectionDiffusionSlipVelocityBoundaryConditions bcParameters,
                                      const real* density, const real* velocityX, const real* velocityY,
                                      const real* velocityZ, const real* turbulentDiffusivity, real diffusivity,
                                      real omegaDiffusivity, const uint* neighborX, const uint* neighborY,
                                      const uint* neighborZ, unsigned long long numberOfLBnodes, bool isEvenTimestep);

template <BoundaryConditionFactory::AdvectionDiffusionDirichletBC bcType>
__global__ void AdvectionDiffusionDirichlet_Device(real* distributionsConcentration,
                                                   AdvectionDiffusionDirichletBoundaryConditions bcParameters,
                                                   const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                                   const real* velocityX, const real* velocityY, const real* velocityZ,
                                                   unsigned long long numberOfLBnodes, real relaxationFrequency,
                                                   bool isEvenTimestep);

template <BoundaryConditionFactory::AdvectionDiffusionNeumannBC bcType>
__global__ void
AdvectionDiffusionNeumann_Device(real* distributionsConcentration, AdvectionDiffusionNeumannBoundaryConditions bcParameters,
                                 const uint* neighborX, const uint* neighborY, const uint* neighborZ, const real* velocityX,
                                 const real* velocityY, const real* velocityZ, unsigned long long numberOfLBnodes,
                                 real relaxationFrequency, bool isEvenTimestep);

__global__ void AdvectionDiffusionBounceBack_Device(real* distributions,
                                                    AdvectionDiffusionNoSlipBoundaryConditions bcParameters,
                                                    const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                                    unsigned long long numberOfLBnodes, bool isEvenTimestep);

#endif

//! \}
