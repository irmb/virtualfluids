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
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file BoundaryConditionFactory.h
//! \ingroup Factories
//! \author Anna Wellmann
//=======================================================================================
#ifndef BC_FACTORY
#define BC_FACTORY

#include <functional>
#include <map>
#include <string>
#include <variant>

#include "LBM/LB.h"
#include "Parameter/Parameter.h"
#include "gpu/GridGenerator/grid/BoundaryConditions/Side.h"


struct LBMSimulationParameter;
class Parameter;

using boundaryCondition = std::function<void(LBMSimulationParameter *, QforBoundaryConditions *)>;
using boundaryConditionWithParameter = std::function<void(Parameter *, QforBoundaryConditions *, const int level)>;
using precursorBoundaryConditionFunc = std::function<void(LBMSimulationParameter *, QforPrecursorBoundaryConditions *, real timeRatio, real velocityRatio)>;

class BoundaryConditionFactory
{
public:
    //! \brief An enumeration for selecting a velocity boundary condition
    enum class VelocityBC {
        //! - VelocitySimpleBounceBackCompressible = plain bounce back velocity boundary condition
        VelocityBounceBack,
        //! - VelocityIncompressible = interpolated velocity boundary condition, based on subgrid distances
        VelocityInterpolatedIncompressible,
        //! - VelocityCompressible = interpolated velocity boundary condition, based on subgrid distances
        VelocityInterpolatedCompressible,
        //! - VelocityAndPressureCompressible = interpolated velocity boundary condition, based on subgrid distances.
        //! Also sets the pressure to the bulk pressure. Can be combined with OutflowNonReflective
        VelocityWithPressureInterpolatedCompressible,
        //! - NotSpecified =  the user did not set a boundary condition
        NotSpecified
    };

    //! \brief An enumeration for selecting a no-slip boundary condition
    enum class NoSlipBC {
        //! - NoSlipDelayBounceBack = implicit bounce back by Esoteric Twist
        NoSlipDelayBounceBack,
        //! - NoSlipBounceBack = explicit bounce back
        NoSlipBounceBack,
        //! - NoSlipInterpolatedIncompressible = interpolated no-slip boundary condition, based on subgrid distances
        NoSlipInterpolatedIncompressible,
        //! - NoSlipInterpolatedCompressible = interpolated no-slip boundary condition, based on subgrid distances
        NoSlipInterpolatedCompressible,
    };

    //! \brief An enumeration for selecting a slip boundary condition
    enum class SlipBC {
        //! - SlipIncompressible = interpolated slip boundary condition, based on subgrid distances
        SlipIncompressible,
        //! - SlipCompressible = interpolated slip boundary condition, based on subgrid distances
        SlipCompressible,
        //! - SlipBounceBack = simple bounce-back slip boundary condition.
        SlipBounceBack,
        //! With turbulent viscosity -> para->setUseTurbulentViscosity(true) has to be set to true
        SlipCompressibleTurbulentViscosity,
        //! With turbulent viscosity -> para->setUseTurbulentViscosity(true) has to be set to true
        SlipPressureCompressibleTurbulentViscosity,
        //! - NotSpecified =  the user did not set a boundary condition
        NotSpecified
    };

    //! \brief An enumeration for selecting a pressure boundary condition
    enum class PressureBC {
        //! - PressureNonEquilibriumIncompressible = pressure boundary condition based on non-equilibrium
        PressureNonEquilibriumIncompressible,
        //! - PressureNonEquilibriumCompressible = pressure boundary condition based on non-equilibrium
        PressureNonEquilibriumCompressible,
        //! - OutflowNonReflective = outflow boundary condition, should be combined with VelocityAndPressureCompressible
        OutflowNonReflective,
        //! - OutflowNonreflectivePressureCorrection = like OutflowNonReflective, but also reduces pressure overshoot
        OutflowNonReflectivePressureCorrection,
        //! - NotSpecified =  the user did not set a boundary condition
        NotSpecified
    };

    //! \brief An enumeration for selecting a stress boundary condition
    enum class StressBC {
        //! - StressCompressible
        StressCompressible,
        //! - StressBounceBack
        StressBounceBack,
        //! - StressPressureBounceBack
        StressPressureBounceBack,
        //! - NotSpecified =  the user did not set a boundary condition
        NotSpecified
    };

    // enum class OutflowBoundaryCondition {};  // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    enum class PrecursorBC {
        //! - VelocityPrecursor
        VelocityPrecursor,
        //! - DisitributionsPrecursor
        DistributionsPrecursor,
        //! - NotSpecified =  the user did not set a boundary condition
        NotSpecified
    };

    void setVelocityBoundaryCondition(const BoundaryConditionFactory::VelocityBC boundaryConditionType);
    void setNoSlipBoundaryCondition(const BoundaryConditionFactory::NoSlipBC boundaryConditionType);
    void setSlipBoundaryCondition(const BoundaryConditionFactory::SlipBC boundaryConditionType);
    void setPressureBoundaryCondition(const BoundaryConditionFactory::PressureBC boundaryConditionType);
    void setStressBoundaryCondition(const BoundaryConditionFactory::StressBC boundaryConditionType);
    void setPrecursorBoundaryCondition(const BoundaryConditionFactory::PrecursorBC boundaryConditionType);
    //! \brief set a boundary condition for the geometry
    //! param boundaryConditionType: a velocity, no-slip or slip boundary condition
    //! \details suggestions for boundaryConditionType:
    //!
    //! - velocity: VelocityIncompressible, VelocityCompressible, VelocityAndPressureCompressible
    //!
    //! - no-slip: NoSlipBounceBack, NoSlipIncompressible, NoSlipCompressible, NoSlip3rdMomentsCompressible
    //!
    //! - slip: only use a slip boundary condition which sets the normals
    void setGeometryBoundaryCondition(const std::variant<VelocityBC, NoSlipBC, SlipBC> boundaryConditionType);

    // void setOutflowBoundaryCondition(...); // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    [[nodiscard]] virtual boundaryCondition getVelocityBoundaryConditionPost(bool isGeometryBC = false) const;
    [[nodiscard]] boundaryCondition getNoSlipBoundaryConditionPost(bool isGeometryBC = false) const;
    [[nodiscard]] boundaryCondition getSlipBoundaryConditionPost(bool isGeometryBC = false) const;
    [[nodiscard]] boundaryCondition getPressureBoundaryConditionPre() const;
    [[nodiscard]] boundaryCondition getGeometryBoundaryConditionPost() const;

    [[nodiscard]] boundaryConditionWithParameter getStressBoundaryConditionPost() const;
    [[nodiscard]] precursorBoundaryConditionFunc getPrecursorBoundaryConditionPost() const;

private:
    VelocityBC velocityBoundaryCondition = VelocityBC::NotSpecified;
    NoSlipBC noSlipBoundaryCondition = NoSlipBC::NoSlipDelayBounceBack;
    SlipBC slipBoundaryCondition = SlipBC::NotSpecified;
    PressureBC pressureBoundaryCondition = PressureBC::NotSpecified;
    std::variant<VelocityBC, NoSlipBC, SlipBC> geometryBoundaryCondition = NoSlipBC::NoSlipDelayBounceBack;
    StressBC stressBoundaryCondition = StressBC::NotSpecified;
    PrecursorBC precursorBoundaryCondition = PrecursorBC::NotSpecified;

    // OutflowBoundaryConditon outflowBC // TODO: https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16
};

#endif
