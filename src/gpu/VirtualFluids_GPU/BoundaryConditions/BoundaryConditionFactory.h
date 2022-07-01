#ifndef BC_FACTORY
#define BC_FACTORY

#include <functional>
#include <map>
#include <string>

#include "LBM/LB.h"

struct LBMSimulationParameter;

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

    //! \brief An enumeration for selecting a no-slip boundary condition
    enum class NoSlipBC {
        //! - NoSlipBounceBack = bounce back no-slip boundary condition
        NoSlipBounceBack,
        //! - NoSlipIncompressible = interpolated no-slip boundary condition, based on subgrid distances
        NoSlipIncompressible,
        //! - NoSlipCompressible = interpolated no-slip boundary condition, based on subgrid distances
        NoSlipCompressible
    };

    // enum class OutflowBoundaryCondition {};  // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    void setVelocityBoundaryCondition(const VelocityBC boundaryConditionType);
    void setNoSlipBoundaryCondition(const NoSlipBC boundaryConditionType);
    // void setGeometryBoundaryCondition(const std::variant<VelocityBC, NoSlipBC, SlipBC> boundaryConditionType);

    // void setOutflowBoundaryCondition(...); // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    boundaryCondition getVelocityBoundaryConditionPost() const;
    boundaryCondition getNoSlipBoundaryConditionPost() const;

private:
    VelocityBC velocityBoundaryCondition;
    NoSlipBC noSlipBoundaryCondition;

    // OutflowBoundaryConditon outflowBC // TODO: https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16
};

#endif
