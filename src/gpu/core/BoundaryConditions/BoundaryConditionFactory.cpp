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
//! \author Martin Schoenherr
//=======================================================================================
#include "BoundaryConditionFactory.h"

#include <variant>

#include <GridGenerator/grid/BoundaryConditions/BoundaryCondition.h>

#include "BoundaryConditions/Outflow/Outflow.h"
#include "BoundaryConditions/Pressure/Pressure.h"
#include "BoundaryConditions/NoSlip/NoSlip.h"
#include "BoundaryConditions/Velocity/Velocity.h"
#include "BoundaryConditions/Slip/Slip.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

void BoundaryConditionFactory::setVelocityBoundaryCondition(VelocityBC boundaryConditionType)
{
    this->velocityBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setNoSlipBoundaryCondition(const NoSlipBC boundaryConditionType)
{
    this->noSlipBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setSlipBoundaryCondition(const SlipBC boundaryConditionType)
{
    this->slipBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setPressureBoundaryCondition(const PressureBC boundaryConditionType)
{
    this->pressureBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setGeometryBoundaryCondition(
    const std::variant<VelocityBC, NoSlipBC, SlipBC> boundaryConditionType)
{
    this->geometryBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setStressBoundaryCondition(const StressBC boundaryConditionType)
{
    this->stressBoundaryCondition = boundaryConditionType;
}

void BoundaryConditionFactory::setPrecursorBoundaryCondition(const PrecursorBC boundaryConditionType)
{
    this->precursorBoundaryCondition = boundaryConditionType;
}

boundaryCondition BoundaryConditionFactory::getVelocityBoundaryConditionPost(bool isGeometryBC) const
{
    const VelocityBC &boundaryCondition =
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

boundaryCondition BoundaryConditionFactory::getNoSlipBoundaryConditionPost(bool isGeometryBC) const
{
    const NoSlipBC &boundaryCondition =
        isGeometryBC ? std::get<NoSlipBC>(this->geometryBoundaryCondition) : this->noSlipBoundaryCondition;

    // for descriptions of the boundary conditions refer to the header
    switch (boundaryCondition) {
        case NoSlipBC::NoSlipDelayBounceBack:
            return [](LBMSimulationParameter *, QforBoundaryConditions *) {};
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

boundaryCondition BoundaryConditionFactory::getSlipBoundaryConditionPost(bool isGeometryBC) const
{
    const SlipBC &boundaryCondition =
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

boundaryCondition BoundaryConditionFactory::getPressureBoundaryConditionPre() const
{
    // for descriptions of the boundary conditions refer to the header
    switch (this->pressureBoundaryCondition) {
        case PressureBC::PressureNonEquilibriumIncompressible:
            return PressureNonEquilibriumIncompressible;
            break;
        case PressureBC::PressureNonEquilibriumCompressible:
            return PressureNonEquilibriumCompressible;
            break;
        case PressureBC::OutflowNonReflective:
            return OutflowNonReflecting;
            break;
        case PressureBC::OutflowNonReflectivePressureCorrection:
            return OutflowNonReflectingPressureCorrection;
        default:
            return nullptr;
    }
}

precursorBoundaryConditionFunc BoundaryConditionFactory::getPrecursorBoundaryConditionPost() const
{
    switch (this->precursorBoundaryCondition) {
        case PrecursorBC::VelocityPrecursor:
            return QPrecursorDevCompZeroPress;
            break;
        case PrecursorBC::DistributionsPrecursor:
            return PrecursorDevDistributions;
            break;
        default:
            return nullptr;
    }
}

boundaryConditionWithParameter BoundaryConditionFactory::getStressBoundaryConditionPost() const
{
    switch (this->stressBoundaryCondition) {
        case StressBC::StressBounceBack:
            return BBStressDev27;
            break;
        case StressBC::StressPressureBounceBack:
            return BBStressPressureDev27;
            break;
        case StressBC::StressCompressible:
            return QStressDevComp27;
            break;
        default:
            return nullptr;
    }
}

boundaryCondition BoundaryConditionFactory::getGeometryBoundaryConditionPost() const
{
    if (std::holds_alternative<VelocityBC>(this->geometryBoundaryCondition))
        return this->getVelocityBoundaryConditionPost(true);
    else if (std::holds_alternative<NoSlipBC>(this->geometryBoundaryCondition))
        return this->getNoSlipBoundaryConditionPost(true);
    else if (std::holds_alternative<SlipBC>(this->geometryBoundaryCondition))
        return this->getSlipBoundaryConditionPost(true);
    return nullptr;
}
