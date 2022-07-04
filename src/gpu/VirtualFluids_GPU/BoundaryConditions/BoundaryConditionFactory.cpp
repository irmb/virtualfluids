#include "BoundaryConditionFactory.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "grid/BoundaryConditions/BoundaryCondition.h"

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

boundaryCondition BoundaryConditionFactory::getVelocityBoundaryConditionPost() const
{
    // for descriptions of the boundary conditions refer to the header
    switch (this->velocityBoundaryCondition) {
        case VelocityBC::VelocitySimpleBounceBackCompressible:
            return QVelDevicePlainBB27;
            break;
        case VelocityBC::VelocityIncompressible:
            return QVelDev27;
            break;
        case VelocityBC::VelocityCompressible:
            return QVelDevComp27;
            break;
        case VelocityBC::VelocityAndPressureCompressible:
            return QVelDevCompZeroPress27;
            break;
        default:
            return nullptr;
    }
}

boundaryCondition BoundaryConditionFactory::getNoSlipBoundaryConditionPost() const
{
    // for descriptions of the boundary conditions refer to the header
    switch (this->noSlipBoundaryCondition) {
        case NoSlipBC::NoSlipBounceBack:
            return BBDev27;
            break;
        case NoSlipBC::NoSlipIncompressible:
            return QDev27;
            break;
        case NoSlipBC::NoSlipCompressible:
            return QDevComp27;
            break;
        default:
            return nullptr;
    }
}

boundaryCondition BoundaryConditionFactory::getSlipBoundaryConditionPost() const
{
    // for descriptions of the boundary conditions refer to the header
    switch (this->slipBoundaryCondition) {
        case SlipBC::SlipIncompressible:
            return QSlipDev27;
            break;
        case SlipBC::SlipCompressible:
            return QSlipDevComp27;
            break;
        case SlipBC::SlipCompressibleTurbulentViscosity:
            return QSlipDevCompTurbulentViscosity27;
            break;
        default:
            return nullptr;
    }
}

boundaryCondition BoundaryConditionFactory::getPressureBoundaryConditionPre() const
{
    // for descriptions of the boundary conditions refer to the header
    switch (this->pressureBoundaryCondition) {
        case PressureBC::PressureEquilibrium:
            return QPressDev27;
            break;
        case PressureBC::PressureEquilibrium2:
            return QPressDevEQZ27;
            break;
        case PressureBC::PressureNonEquilibriumIncompressible:
            return QPressDevIncompNEQ27;
            break;
        case PressureBC::PressureNonEquilibriumCompressible:
            return QPressDevNEQ27;
            break;
        case PressureBC::OutflowNonReflective:
            return QPressNoRhoDev27;
            break;
        default:
            return nullptr;
    }
}

// boundaryCondition BoundaryConditionFactory::getGeometryBoundaryConditionPost() const{
//     this->getNoSlipBoundaryConditionPost();
// }