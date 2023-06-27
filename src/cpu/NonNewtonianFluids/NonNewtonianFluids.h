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
//! \file NonNewtonianFluids.h
//! \ingroup NonNewtonianFluids
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef NonNewtonianFluids_h
#define NonNewtonianFluids_h

#include "NonNewtonianFluids/BoundaryConditions/ThixotropyDensityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyVelocityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyNonReflectingOutflowBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/ThixotropyVelocityWithDensityBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyBinghamModelNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyHerschelBulkleyModelNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyPowellEyringModelNoSlipBCStrategy.h"
#include "NonNewtonianFluids/BoundaryConditions/RheologyBinghamModelVelocityBCStrategy.h"

#include "NonNewtonianFluids/SimulationObservers/CalculateTorqueSimulationObserver.h"
#include "NonNewtonianFluids/SimulationObservers/WriteThixotropyQuantitiesSimulationObserver.h"

#include "NonNewtonianFluids/LBM/ThixotropyLBMKernel.h"
#include "NonNewtonianFluids/LBM/ThixotropyExpLBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyBinghamModelLBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyHerschelBulkleyModelLBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyInterpolator.h"
#include "NonNewtonianFluids/LBM/Rheology.h"
#include "NonNewtonianFluids/LBM/RheologyK17LBMKernel.h"
#include "NonNewtonianFluids/LBM/RheologyPowellEyringModelLBMKernel.h"

#endif