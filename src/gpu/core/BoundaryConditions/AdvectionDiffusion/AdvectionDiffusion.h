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
//! \author Martin Schoenherr, Henry Korb
//=======================================================================================
#ifndef AdvectionDiffusion_H
#define AdvectionDiffusion_H

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

#include <cuda.h>
#include <cuda_runtime.h>

struct LBMSimulationParameter;
class Parameter;

void AdvectionDiffusionBounceBack(LBMSimulationParameter* parameterDevice,
                                  AdvectionDiffusionNoSlipBoundaryConditions bcParameters);
void AdvectionDiffusionFluxTurbulentViscosityCompressible(LBMSimulationParameter* parameterDevice,
                                                          AdvectionDiffusionFluxBoundaryConditions bcParameters);
void AdvectionDiffusionFluxCompressible(LBMSimulationParameter* parameterDevice,
                                        AdvectionDiffusionFluxBoundaryConditions bcParameters);
void AdvectionDiffusionFluxBounceBack(LBMSimulationParameter* parameterDevice,
                                      AdvectionDiffusionFluxBoundaryConditions bcParameters);
void AdvectionDiffusionDirichletAntiBounceBackSlip(LBMSimulationParameter* parameterDevice,
                                                   AdvectionDiffusionDirichletBoundaryConditions bcParameters);
void AdvectionDiffusionDirichletInterpolatedSlip(LBMSimulationParameter* parameterDevice,
                                                 AdvectionDiffusionDirichletBoundaryConditions bcParameters);
void AdvectionDiffusionDirichletAntiBounceBackNoSlip(LBMSimulationParameter* parameterDevice,
                                                     AdvectionDiffusionDirichletBoundaryConditions bcParameters);
void AdvectionDiffusionDirichletInterpolatedNoSlip(LBMSimulationParameter* parameterDevice,
                                                   AdvectionDiffusionDirichletBoundaryConditions bcParameters);
void AdvectionDiffusionNeumannAntiBounceBackSlip(LBMSimulationParameter* parameterDevice,
                                                 AdvectionDiffusionNeumannBoundaryConditions bcParameters);
void AdvectionDiffusionNeumannInterpolatedSlip(LBMSimulationParameter* parameterDevice,
                                               AdvectionDiffusionNeumannBoundaryConditions bcParameters);
void AdvectionDiffusionNeumannAntiBounceBackNoSlip(LBMSimulationParameter* parameterDevice,
                                                   AdvectionDiffusionNeumannBoundaryConditions bcParameters);
void AdvectionDiffusionNeumannInterpolatedNoSlip(LBMSimulationParameter* parameterDevice,
                                                 AdvectionDiffusionNeumannBoundaryConditions bcParameters);

#endif

//! \}
