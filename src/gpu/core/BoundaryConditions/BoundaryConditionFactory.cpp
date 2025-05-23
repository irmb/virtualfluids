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
#include "BoundaryConditionFactory.h"

#include <variant>

#include <GridGenerator/grid/BoundaryConditions/BoundaryCondition.h>

#include "BoundaryConditions/NoSlip/NoSlip.h"
#include "BoundaryConditions/Outflow/Outflow.h"
#include "BoundaryConditions/Precursor/Precursor.h"
#include "BoundaryConditions/Pressure/Pressure.h"
#include "BoundaryConditions/Slip/Slip.h"
#include "BoundaryConditions/Stress/Stress.h"
#include "BoundaryConditions/Velocity/Velocity.h"
#include "BoundaryConditions/AdvectionDiffusion/AdvectionDiffusion.h"
#include "Parameter/Parameter.h"

void BoundaryConditionFactory::setVelocityBoundaryCondition(VelocityBC boundaryConditionType)
{
    this->velocityBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setNoSlipBoundaryCondition(NoSlipBC boundaryConditionType)
{
    this->noSlipBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setSlipBoundaryCondition(SlipBC boundaryConditionType)
{
    this->slipBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setPressureBoundaryCondition(PressureBC boundaryConditionType)
{
    this->pressureBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setGeometryBoundaryCondition(std::variant<VelocityBC, NoSlipBC, SlipBC> boundaryConditionType)
{
    this->geometryBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setStressBoundaryCondition(StressBC boundaryConditionType)
{
    this->stressBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setPrecursorBoundaryCondition(PrecursorBC boundaryConditionType)
{
    this->precursorBoundaryCondition = boundaryConditionType;
}
void BoundaryConditionFactory::setAdvectionDiffusionNoSlipBoundaryCondition(AdvectionDiffusionNoSlipBC boundaryConditionType)
{
    this->advectionDiffusionNoSlipBoundaryCondition = boundaryConditionType;
}
void BoundaryConditionFactory::setAdvectionDiffusionSlipVelocityBoundaryCondition(AdvectionDiffusionSlipVelocityBC boundaryConditionType)
{
    this->advectionDiffusionSlipVelocityBoundaryCondition = boundaryConditionType;
}
void BoundaryConditionFactory::setAdvectionDiffusionDirichletBoundaryCondition(AdvectionDiffusionDirichletBC boundaryConditionType)
{
    this->advectionDiffusionDirichletBoundaryCondition = boundaryConditionType;
}
void BoundaryConditionFactory::setAdvectionDiffusionNeumannBoundaryCondition(AdvectionDiffusionNeumannBC boundaryConditionType)
{
    this->advectionDiffusionNeumannBoundaryCondition = boundaryConditionType;
}
BoundaryConditionKernel BoundaryConditionFactory::getVelocityBoundaryConditionPost(bool isGeometryBC) const
{
    const VelocityBC& boundaryCondition =
        isGeometryBC ? std::get<VelocityBC>(this->geometryBoundaryCondition) : this->velocityBoundaryCondition;

    // for descriptions of the boundary conditions refer to the header
    switch (boundaryCondition) {
        case VelocityBC::VelocityBounceBack:
            return VelocityBounceBack;
            break;
        case VelocityBC::VelocityInterpolatedIncompressible:
            return VelocityInterpolatedIncompressible;
            break;
        case VelocityBC::VelocityInterpolatedCompressible:
            return VelocityInterpolatedCompressible;
            break;
        case VelocityBC::VelocityWithPressureInterpolatedCompressible:
            return VelocityWithPressureInterpolatedCompressible;
            break;
        default:
            return nullptr;
    }
}

BoundaryConditionKernel BoundaryConditionFactory::getNoSlipBoundaryConditionPost(bool isGeometryBC) const
{
    const NoSlipBC& boundaryCondition =
        isGeometryBC ? std::get<NoSlipBC>(this->geometryBoundaryCondition) : this->noSlipBoundaryCondition;

    // for descriptions of the boundary conditions refer to the header
    switch (boundaryCondition) {
        case NoSlipBC::NoSlipDelayBounceBack:
            return [](LBMSimulationParameter*, QforBoundaryConditions*) {};
            break;
        case NoSlipBC::NoSlipBounceBack:
            return NoSlipBounceBack;
            break;
        case NoSlipBC::NoSlipInterpolatedIncompressible:
            return NoSlipInterpolatedIncompressible;
            break;
        case NoSlipBC::NoSlipInterpolatedCompressible:
            return NoSlipInterpolatedCompressible;
            break;
        default:
            return nullptr;
    }
}

BoundaryConditionKernel BoundaryConditionFactory::getSlipBoundaryConditionPost(bool isGeometryBC) const
{
    const SlipBC& boundaryCondition =
        isGeometryBC ? std::get<SlipBC>(this->geometryBoundaryCondition) : this->slipBoundaryCondition;

    // for descriptions of the boundary conditions refer to the header
    switch (boundaryCondition) {
        case SlipBC::SlipCompressible:
            return SlipCompressible;
            break;
        case SlipBC::SlipTurbulentViscosityCompressible:
            return SlipTurbulentViscosityCompressible;
            break;
        default:
            return nullptr;
    }
}

std::variant<BoundaryConditionKernel, DirectionalBoundaryConditionKernel>
BoundaryConditionFactory::getPressureBoundaryConditionPre() const
{
    // for descriptions of the boundary conditions refer to the header
    switch (this->pressureBoundaryCondition) {
        case PressureBC::PressureNonEquilibriumIncompressible:
            return (DirectionalBoundaryConditionKernel)PressureNonEquilibriumIncompressible;
            break;
        case PressureBC::PressureNonEquilibriumCompressible:
            return (DirectionalBoundaryConditionKernel)PressureNonEquilibriumCompressible;
            break;
        case PressureBC::OutflowNonReflective:
            return (DirectionalBoundaryConditionKernel)OutflowNonReflecting;
            break;
        case PressureBC::OutflowNonReflectivePressureCorrection:
            return (DirectionalBoundaryConditionKernel)OutflowNonReflectingPressureCorrection;
        default:
            return (BoundaryConditionKernel) nullptr;
    }
}

bool BoundaryConditionFactory::hasDirectionalPressureBoundaryCondition() const
{
    return std::holds_alternative<DirectionalBoundaryConditionKernel>(getPressureBoundaryConditionPre());
}

PrecursorBoundaryConditionKernel BoundaryConditionFactory::getPrecursorBoundaryConditionPost() const
{
    switch (this->precursorBoundaryCondition) {
        case PrecursorBC::PrecursorNonReflectiveCompressible:
            return PrecursorNonReflectiveCompressible;
            break;
        case PrecursorBC::PrecursorDistributions:
            return PrecursorDistributions;
            break;
        default:
            return nullptr;
    }
}

BoundaryConditionWithParameterKernel BoundaryConditionFactory::getStressBoundaryConditionPost() const
{
    switch (this->stressBoundaryCondition) {
        case StressBC::StressBounceBackCompressible:
            return StressBounceBackCompressible;
            break;
        case StressBC::StressBounceBackPressureCompressible:
            return StressBounceBackPressureCompressible;
            break;
        case StressBC::StressCompressible:
            return StressCompressible;
            break;
        default:
            return nullptr;
    }
}

BoundaryConditionKernel BoundaryConditionFactory::getGeometryBoundaryConditionPost() const
{
    if (std::holds_alternative<VelocityBC>(this->geometryBoundaryCondition))
        return this->getVelocityBoundaryConditionPost(true);
    if (std::holds_alternative<NoSlipBC>(this->geometryBoundaryCondition))
        return this->getNoSlipBoundaryConditionPost(true);
    if (std::holds_alternative<SlipBC>(this->geometryBoundaryCondition))
        return this->getSlipBoundaryConditionPost(true);
    return nullptr;
}

AdvectionDiffusionNoSlipBoundaryConditionKernel BoundaryConditionFactory::getAdvectionDiffusionNoSlipBoundaryConditionPost() const
{
    if(this->advectionDiffusionNoSlipBoundaryCondition == AdvectionDiffusionNoSlipBC::NoSlipDelayedBounceBack)
        return [](LBMSimulationParameter*, AdvectionDiffusionNoSlipBoundaryConditions) {};
    if(this->advectionDiffusionNoSlipBoundaryCondition == AdvectionDiffusionNoSlipBC::NoSlipBounceBack)
        return AdvectionDiffusionBounceBack;
    return nullptr;
}

AdvectionDiffusionSlipVelocityBoundaryConditionKernel BoundaryConditionFactory::getAdvectionDiffusionSlipVelocityBoundaryConditionPost() const
{
    if(this->advectionDiffusionSlipVelocityBoundaryCondition == AdvectionDiffusionSlipVelocityBC::SlipVelocityBounceBack)
        return AdvectionDiffusionSlipVelocityBounceBack;
    if(this->advectionDiffusionSlipVelocityBoundaryCondition == AdvectionDiffusionSlipVelocityBC::SlipVelocityCompressible)
        return AdvectionDiffusionSlipVelocityCompressible;
    if(this->advectionDiffusionSlipVelocityBoundaryCondition == AdvectionDiffusionSlipVelocityBC::SlipVelocityTurbulentViscosityCompressible)
        return AdvectionDiffusionSlipVelocityTurbulentViscosityCompressible;
    return nullptr;
}

AdvectionDiffusionDirichletBoundaryConditionKernel BoundaryConditionFactory::getAdvectionDiffusionDirichletBoundaryConditionPost() const
{
    switch(this->advectionDiffusionDirichletBoundaryCondition)
    {
        case AdvectionDiffusionDirichletBC::DirichletAntiBounceBackNoSlip:
            return AdvectionDiffusionDirichletAntiBounceBackNoSlip;
        case AdvectionDiffusionDirichletBC::DirichletInterpolatedNoSlip:
            return AdvectionDiffusionDirichletInterpolatedNoSlip;
        case AdvectionDiffusionDirichletBC::DirichletAntiBounceBackSlip:
            return AdvectionDiffusionDirichletAntiBounceBackSlip;
        case AdvectionDiffusionDirichletBC::DirichletInterpolatedSlip:
            return AdvectionDiffusionDirichletInterpolatedSlip;
        default:
            return nullptr;
    };
}

AdvectionDiffusionNeumannBoundaryConditionKernel BoundaryConditionFactory::getAdvectionDiffusionNeumannBoundaryConditionPost() const
{
    switch(this->advectionDiffusionNeumannBoundaryCondition)
    {
        case AdvectionDiffusionNeumannBC::NeumannAntiBounceBackNoSlip:
            return AdvectionDiffusionNeumannAntiBounceBackNoSlip;
        case AdvectionDiffusionNeumannBC::NeumannInterpolatedNoSlip:
            return AdvectionDiffusionNeumannInterpolatedNoSlip;
        case AdvectionDiffusionNeumannBC::NeumannAntiBounceBackSlip:
            return AdvectionDiffusionNeumannAntiBounceBackSlip;
        case AdvectionDiffusionNeumannBC::NeumannInterpolatedSlip:
            return AdvectionDiffusionNeumannInterpolatedSlip;
        default:
            return nullptr;
    };
}



//! \}
