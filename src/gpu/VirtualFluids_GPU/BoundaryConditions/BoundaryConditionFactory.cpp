#include "BoundaryConditionFactory.h"
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

// boundaryCondition BoundaryConditionFactory::getGeometryBoundaryConditionPost() const{
//     this->getNoSlipBoundaryConditionPost();
// }