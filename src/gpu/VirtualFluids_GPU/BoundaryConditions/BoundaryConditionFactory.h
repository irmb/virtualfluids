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
    //! \brief An enumeration for selecting a velocity boundary condition
    enum class VelocityBC {
        //! - VelocitySimpleBounceBackCompressible = plain bounce back velocity boundary condition
        VelocitySimpleBounceBackCompressible,
        //! - VelocityIncompressible = interpolated velocity boundary condition, based on subgrid distances
        VelocityIncompressible,
        //! - VelocityCompressible = interpolated velocity boundary condition, based on subgrid distances
        VelocityCompressible,
        //! - VelocityAndPressureCompressible = interpolated velocity boundary condition, based on subgrid distances.
        //! Also sets the pressure to the bulk pressure.
        VelocityAndPressureCompressible
    };

    // enum class OutflowBoundaryCondition {};  // TODO:
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
