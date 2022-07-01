#include "BoundaryConditionFactory.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

void BoundaryConditionFactory::setVelocityBoundaryCondition(VelocityBC boundaryConditionType)
{
    this->velocityBoundaryCondition = boundaryConditionType;
}

boundaryCondition BoundaryConditionFactory::getVelocityBoundaryConditionPost() const
{
    switch (this->velocityBoundaryCondition) {
        case VelocityBC::VelocitySimpleBounceBackCompressible:
            // VelocitySimpleBounceBackCompressible = plain bounce back velocity boundary condition
            return QVelDevicePlainBB27;
            break;
        case VelocityBC::VelocityIncompressible:
            // VelocityIncompressible = interpolated velocity boundary condition, based on subgrid distances
            return QVelDev27;
            break;
        case VelocityBC::VelocityCompressible:
            // VelocityCompressible = interpolated velocity boundary condition, based on subgrid distances
            return QVelDevComp27;
            break;
        case VelocityBC::VelocityAndPressureCompressible:
            // VelocityAndPressureCompressible = interpolated velocity boundary condition, based on subgrid distances.
            // Also sets the pressure to the bulk pressure.
            return QVelDevCompZeroPress27;
            break;
        default:
            return nullptr;
    }
}