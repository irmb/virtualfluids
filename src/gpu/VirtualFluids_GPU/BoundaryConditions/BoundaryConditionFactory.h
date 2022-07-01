#ifndef BC_FACTORY
#define BC_FACTORY

#include <functional>
#include <map>
#include <string>

#include "LBM/LB.h"

class LBMSimulationParameter;

using boundaryCondition = std::function<void(LBMSimulationParameter *, QforBoundaryConditions *)>;

class BoundaryConditionFactory
{
public:
    enum class VelocityBC {
        VelocitySimpleBounceBackCompressible,
        VelocityIncompressible,
        VelocityCompressible,
        VelocityAndPressureCompressible
    };
    // enum class OutflowBoundaryConditon {};  // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    void setVelocityBoundaryCondition(const VelocityBC boundaryConditionType);
    // void setOutflowBoundaryCondition(...); // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    boundaryCondition getVelocityBoundaryConditionPost() const;

private:
    VelocityBC velocityBoundaryCondition;
    // OutflowBoundaryConditon outflowBC // TODO: https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16
};

#endif
