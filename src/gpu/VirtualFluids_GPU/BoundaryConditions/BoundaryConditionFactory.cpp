#include "BoundaryConditionFactory.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

void BoundaryConditionFactory::setVelocityBoundaryCondition(VelocityBC boundaryConditionType)
{
    this->velocityBoundaryCondition = boundaryConditionType;
}

boundaryCondition BoundaryConditionFactory::getVelocityBoundaryConditionPost() const
{
    // for descriptions of the boundary conditions refer to the header ( VelocityBC)
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