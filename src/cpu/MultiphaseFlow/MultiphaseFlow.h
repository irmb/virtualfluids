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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file MultiphaseFlow.h
//! \ingroup MultiphaseFlow
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef MultiphaseFlow_h
#define MultiphaseFlow_h

#include "MultiphaseFlow/BoundaryConditions/MultiphaseNoSlipBCStrategy.h"
#include "MultiphaseFlow/BoundaryConditions/MultiphaseNonReflectingOutflowBCStrategy.h"
#include "MultiphaseFlow/BoundaryConditions/MultiphaseVelocityBC.h"
#include "MultiphaseFlow/BoundaryConditions/MultiphaseVelocityBCStrategy.h"
#include "MultiphaseFlow/BoundaryConditions/MultiphaseSlipBCStrategy.h"
#include "MultiphaseFlow/BoundaryConditions/MultiphasePressureBCStrategy.h"
          
#include "MultiphaseFlow/LBM/MultiphaseCumulantLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphasePressureFilterCompressibleAirLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphasePressureFilterLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphaseScaleDistributionLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphaseScratchCumulantLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphaseSharpInterfaceLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphaseSimpleVelocityBaseExternalPressureLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphaseTwoPhaseFieldsCumulantLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphaseTwoPhaseFieldsPressureFilterLBMKernel.h"
#include "MultiphaseFlow/LBM/MultiphaseTwoPhaseFieldsVelocityCumulantLBMKernel.h"
          
#include "MultiphaseFlow/SimulationObservers/WriteMultiphaseQuantitiesSimulationObserver.h"
#include "MultiphaseFlow/SimulationObservers/WriteSharpInterfaceQuantitiesSimulationObserver.h"

#include "MultiphaseFlow/Visitors/MultiphaseSetKernelBlockVisitor.h"
#include "MultiphaseFlow/Visitors/MultiphaseBoundaryConditionsBlockVisitor.h"
#include "MultiphaseFlow/Visitors/MultiphaseInitDistributionsBlockVisitor.h"
#include "MultiphaseFlow/Visitors/MultiphaseVelocityFormInitDistributionsBlockVisitor.h"

#endif