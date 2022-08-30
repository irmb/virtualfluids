#ifndef BC_FACTORY
#define BC_FACTORY

#include <functional>
#include <map>
#include <string>
#include <variant>

#include "LBM/LB.h"
#include "Parameter/Parameter.h"
#include "grid/BoundaryConditions/Side.h"

struct LBMSimulationParameter;
class Parameter;

using boundaryCondition = std::function<void(LBMSimulationParameter *, QforBoundaryConditions *)>;
using boundaryConditionPara = std::function<void(Parameter *, QforBoundaryConditions *, const int level)>;

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
        //! Also sets the pressure to the bulk pressure. Can be combined with OutflowNonReflective
        VelocityAndPressureCompressible
    };

    //! \brief An enumeration for selecting a no-slip boundary condition
    enum class NoSlipBC {
        //! - NoSlipImplicitBounceBack = implicit bounce back by Esoteric Twist
        NoSlipImplicitBounceBack,
        //! - NoSlipBounceBack = bounce back no-slip boundary condition
        NoSlipBounceBack,
        //! - NoSlipIncompressible = interpolated no-slip boundary condition, based on subgrid distances
        NoSlipIncompressible,
        //! - NoSlipCompressible = interpolated no-slip boundary condition, based on subgrid distances
        NoSlipCompressible,
        //! - NoSlipCompressible = interpolated no-slip boundary condition, based on subgrid distances
        //! Also uses the third order moments.
        NoSlip3rdMomentsCompressible
    };

    //! \brief An enumeration for selecting a slip boundary condition
    enum class SlipBC {
        //! - SlipIncompressible = interpolated slip boundary condition, based on subgrid distances
        SlipIncompressible,
        //! - SlipCompressible = interpolated slip boundary condition, based on subgrid distances
        SlipCompressible,
        //! - SlipCompressible = interpolated slip boundary condition, based on subgrid distances.
        //! With turbulent viscosity -> para->setUseTurbulentViscosity(true) has to be set to true
        SlipCompressibleTurbulentViscosity
    };

    //! \brief An enumeration for selecting a pressure boundary condition
    enum class PressureBC {
        //! - PressureEquilibrium = pressure boundary condition based on equilibrium
        PressureEquilibrium, // incorrect pressure :(
        //! - PressureEquilibrium2 = pressure boundary condition based on equilibrium (potentially better?! than PressureEquilibrium)
        PressureEquilibrium2, // is broken --> nan :(
        //! - PressureNonEquilibriumIncompressible = pressure boundary condition based on non-equilibrium
        PressureNonEquilibriumIncompressible,
        //! - PressureNonEquilibriumCompressible = pressure boundary condition based on non-equilibrium
        PressureNonEquilibriumCompressible,
        //! - OutflowNonReflective = outflow boundary condition, should be combined with VelocityAndPressureCompressible
        OutflowNonReflective,
        //! - OutflowNonreflectivePressureCorrection = like OutflowNonReflective, but also reduces pressure overshoot
        OutflowNonReflectivePressureCorrection
    };

    //! \brief An enumeration for selecting a stress boundary condition
    enum class StressBC {
        //! - StressCompressible
        StressCompressible,
        //! - StressBounceBack
        StressBounceBack
    };

    // enum class OutflowBoundaryCondition {};  // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    void setVelocityBoundaryCondition(const BoundaryConditionFactory::VelocityBC boundaryConditionType);
    void setNoSlipBoundaryCondition(const BoundaryConditionFactory::NoSlipBC boundaryConditionType);
    void setSlipBoundaryCondition(const BoundaryConditionFactory::SlipBC boundaryConditionType);
    void setPressureBoundaryCondition(const BoundaryConditionFactory::PressureBC boundaryConditionType);
    void setStressBoundaryCondition(const BoundaryConditionFactory::StressBC boundaryConditionType);
    //!param boundaryConditionType: a velocity, no-slip or slip boundary condition
    //! \details suggestions for boundaryConditionType:
    //!
    //! - velocity: VelocityIncompressible, VelocityCompressible, VelocityAndPressureCompressible
    //!
    //! - no-slip:  NoSlipBounceBack, NoSlipIncompressible, NoSlipCompressible, NoSlip3rdMomentsCompressible
    //!
    //! - slip:     SlipIncompressible
    void setGeometryBoundaryCondition(const std::variant<VelocityBC, NoSlipBC, SlipBC> boundaryConditionType);

    // void setOutflowBoundaryCondition(...); // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    boundaryCondition getVelocityBoundaryConditionPost(bool isGeometryBC = false) const;
    boundaryCondition getNoSlipBoundaryConditionPost(bool isGeometryBC = false) const;
    boundaryCondition getSlipBoundaryConditionPost(bool isGeometryBC = false) const;
    boundaryCondition getPressureBoundaryConditionPre() const;
    boundaryCondition getGeometryBoundaryConditionPost() const;

    boundaryConditionPara getStressBoundaryConditionPost() const;

private:
    VelocityBC velocityBoundaryCondition;
    NoSlipBC noSlipBoundaryCondition = NoSlipBC::NoSlipImplicitBounceBack;
    SlipBC slipBoundaryCondition;
    PressureBC pressureBoundaryCondition;
    std::variant<VelocityBC, NoSlipBC, SlipBC> geometryBoundaryCondition  = NoSlipBC::NoSlipImplicitBounceBack;
    StressBC stressBoundaryCondition;


    // OutflowBoundaryConditon outflowBC // TODO: https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16
};

#endif
