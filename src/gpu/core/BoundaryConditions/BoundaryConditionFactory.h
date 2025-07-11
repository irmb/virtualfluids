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
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_BoundaryConditions BoundaryConditions
//! \ingroup gpu_core core
//! \{
//! \author Anna Wellmann
//=======================================================================================
#ifndef BC_FACTORY
#define BC_FACTORY

#include <functional>
#include <map>
#include <string>
#include <variant>

#include <GridGenerator/grid/BoundaryConditions/Side.h>

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

struct LBMSimulationParameter;
class Parameter;

using BoundaryConditionKernel = std::function<void(LBMSimulationParameter*, QforBoundaryConditions*)>;
using DirectionalBoundaryConditionKernel = std::function<void(LBMSimulationParameter*, QforDirectionalBoundaryCondition*)>;
using BoundaryConditionWithParameterKernel = std::function<void(Parameter*, QforBoundaryConditions*, const int level)>;
using PrecursorBoundaryConditionKernel =
    std::function<void(LBMSimulationParameter*, QforPrecursorBoundaryConditions*, real timeRatio, real velocityRatio)>;

using AdvectionDiffusionNoFluxBoundaryConditionKernel = std::function<void(LBMSimulationParameter *, AdvectionDiffusionNoFluxBoundaryConditions)>;
using AdvectionDiffusionFluxBoundaryConditionKernel = std::function<void(LBMSimulationParameter *, AdvectionDiffusionFluxBoundaryConditions)>;
using AdvectionDiffusionDirichletBoundaryConditionKernel = std::function<void(LBMSimulationParameter *, AdvectionDiffusionDirichletBoundaryConditions)>;
using AdvectionDiffusionNeumannBoundaryConditionKernel = std::function<void(LBMSimulationParameter *, AdvectionDiffusionNeumannBoundaryConditions)>;

class BoundaryConditionFactory
{
public:
    virtual ~BoundaryConditionFactory() = default;
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
        //! SlipBounceBack = slip boundary condition based on bounce back
        SlipBounceBack,
        //! - SlipCompressible = interpolated slip boundary condition, based on subgrid distances
        SlipCompressible,
        //! With turbulent viscosity -> para->setUseTurbulentViscosity(true) has to be set to true
        SlipTurbulentViscosityCompressible,
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
        //! - StressBounceBackCompressible
        StressBounceBackCompressible,
        //! - StressBounceBackWithPressureCompressible
        StressBounceBackWithPressureCompressible,
        //! - StressInterpolatedCompressible
        StressInterpolatedCompressible,
        //! - NotSpecified =  the user did not set a boundary condition
        NotSpecified
    };

    // enum class OutflowBoundaryCondition {};  // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    enum class PrecursorBC {
        //! - VelocityPrecursor
        PrecursorNonReflectiveCompressible,
        //! - DistributionsPrecursor
        PrecursorDistributions,
        //! - NotSpecified =  the user did not set a boundary condition
        NotSpecified
    };

    //! \brief Equivalent to an adiabatic boundary condition, best used in combination with NoSlip
    enum class AdvectionDiffusionNoFluxBC { 
        //! NoFluxBounceBackDelayed = implicit bounce back
        NoFluxDelayedBounceBack, 
        //! NoFluxBounceBack = simple bounce back
        NoFluxBounceBack,
    };

    //! \brief Can set flux, best used in combination with Slip or velocity. 
    // Works well at high velocities.
    enum class AdvectionDiffusionFluxBC {
        FluxTurbulentViscosityCompressible,
        FluxCompressible,
        FluxBounceBack,
        NotSpecified
    };

    //! \brief Sets constant value at boundary via Anti bounce back rule
    enum class AdvectionDiffusionDirichletBC {
        //! - Interpolated Dirichlet boundary condition, uses subgrid distances, must be used with Slip
        DirichletInterpolatedSlip,
        //! - Bounce Back Dirichlet boundary condition, does not use subgrid distances, must be used with Slip
        DirichletAntiBounceBackSlip,
        //! - Interpolated Dirichlet boundary condition, uses subgrid distances, must be used with NoSlip
        DirichletInterpolatedNoSlip,
        //! - BounceBack Dirichlet boundary condition, does not use subgrid distances, must be used with NoSlip
        DirichletAntiBounceBackNoSlip,
        NotSpecified
    };

    //! \brief Sets gradient at boundary via anti bounce back rule. Only works well at low velocities.
    enum class AdvectionDiffusionNeumannBC {
        //! - Interpolated Neumann boundary condition, uses subgrid distances, must be used with Slip
        NeumannInterpolatedSlip,
        //! - Interpolated Neumann boundary condition, does not use subgrid distances, must be used with Slip
        NeumannAntiBounceBackSlip,
        //! - Interpolated Neumann boundary condition, uses subgrid distances, must be used with NoSlip
        NeumannInterpolatedNoSlip,
        //! - Interpolated Neumann boundary condition, does not use subgrid distances, must be used with NoSlip
        NeumannAntiBounceBackNoSlip,
        NotSpecified
    };

    //! \brief Enum to differentiate between setting heatlfux or surfaceTemperature in SurfaceLayer BC
    enum class SurfaceLayerBC {
        SurfaceHeatFlux,
        SurfaceTemperature,
        NotSpecified
    };

    void setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC boundaryConditionType);
    void setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC boundaryConditionType);
    void setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC boundaryConditionType);
    void setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC boundaryConditionType);
    void setStressBoundaryCondition(BoundaryConditionFactory::StressBC boundaryConditionType);
    void setPrecursorBoundaryCondition(BoundaryConditionFactory::PrecursorBC boundaryConditionType);
    //! \brief set a boundary condition for the geometry
    //! param boundaryConditionType: a velocity, no-slip or slip boundary condition
    //! \details suggestions for boundaryConditionType:
    //!
    //! - velocity: VelocityIncompressible, VelocityCompressible, VelocityAndPressureCompressible
    //!
    //! - no-slip: NoSlipBounceBack, NoSlipIncompressible, NoSlipCompressible, NoSlip3rdMomentsCompressible
    //!
    //! - slip: only use a slip boundary condition which sets the normals
    void setGeometryBoundaryCondition(std::variant<VelocityBC, NoSlipBC, SlipBC> boundaryConditionType);
    void setAdvectionDiffusionNoFluxBoundaryCondition(AdvectionDiffusionNoFluxBC boundaryConditionType);
    void setAdvectionDiffusionFluxBoundaryCondition(AdvectionDiffusionFluxBC boundaryConditionType);
    void setAdvectionDiffusionDirichletBoundaryCondition(AdvectionDiffusionDirichletBC boundaryConditionType);
    void setAdvectionDiffusionNeumannBoundaryCondition(AdvectionDiffusionNeumannBC boundaryConditionType);
    void setSurfaceLayerBoundaryCondition(StressBC momentumBoundaryConditionType, SurfaceLayerBC surfaceLayerBoundaryConditionType);
    // void setOutflowBoundaryCondition(...); // TODO:
    // https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16

    [[nodiscard]] virtual BoundaryConditionKernel getVelocityBoundaryConditionPost(bool isGeometryBC = false) const;
    [[nodiscard]] BoundaryConditionKernel getNoSlipBoundaryConditionPost(bool isGeometryBC = false) const;
    [[nodiscard]] BoundaryConditionKernel getSlipBoundaryConditionPost(bool isGeometryBC = false) const;
    [[nodiscard]] BoundaryConditionKernel getGeometryBoundaryConditionPost() const;
    [[nodiscard]] virtual std::variant<BoundaryConditionKernel, DirectionalBoundaryConditionKernel> getPressureBoundaryConditionPre() const;
    [[nodiscard]] BoundaryConditionKernel getStressBoundaryConditionPost() const;
    [[nodiscard]] PrecursorBoundaryConditionKernel getPrecursorBoundaryConditionPost() const;
    [[nodiscard]] AdvectionDiffusionNoFluxBoundaryConditionKernel getAdvectionDiffusionNoFluxBoundaryConditionPost() const;
    [[nodiscard]] AdvectionDiffusionFluxBoundaryConditionKernel getAdvectionDiffusionFluxBoundaryConditionPost() const;
    [[nodiscard]] AdvectionDiffusionDirichletBoundaryConditionKernel getAdvectionDiffusionDirichletBoundaryConditionPost() const;
    [[nodiscard]] AdvectionDiffusionNeumannBoundaryConditionKernel getAdvectionDiffusionNeumannBoundaryConditionPost() const;
    [[nodiscard]] BoundaryConditionKernel getSurfaceLayerBoundaryConditionPost() const;
    [[nodiscard]] virtual bool hasDirectionalPressureBoundaryCondition() const;

private:
    VelocityBC velocityBoundaryCondition = VelocityBC::NotSpecified;
    NoSlipBC noSlipBoundaryCondition = NoSlipBC::NoSlipDelayBounceBack;
    SlipBC slipBoundaryCondition = SlipBC::NotSpecified;
    PressureBC pressureBoundaryCondition = PressureBC::NotSpecified;
    std::variant<VelocityBC, NoSlipBC, SlipBC> geometryBoundaryCondition = NoSlipBC::NoSlipDelayBounceBack;
    StressBC stressBoundaryCondition = StressBC::NotSpecified;
    PrecursorBC precursorBoundaryCondition = PrecursorBC::NotSpecified;
    AdvectionDiffusionNoFluxBC advectionDiffusionNoFluxBoundaryCondition = AdvectionDiffusionNoFluxBC::NoFluxDelayedBounceBack;
    AdvectionDiffusionFluxBC advectionDiffusionFluxBoundaryCondition = AdvectionDiffusionFluxBC::NotSpecified;
    AdvectionDiffusionDirichletBC advectionDiffusionDirichletBoundaryCondition = AdvectionDiffusionDirichletBC::NotSpecified;
    AdvectionDiffusionNeumannBC advectionDiffusionNeumannBoundaryCondition = AdvectionDiffusionNeumannBC::NotSpecified;
    std::pair<StressBC, SurfaceLayerBC> surfaceLayerBoundaryCondition = {StressBC::NotSpecified, SurfaceLayerBC::NotSpecified};

    // OutflowBoundaryConditon outflowBC // TODO: https://git.rz.tu-bs.de/m.schoenherr/VirtualFluids_dev/-/issues/16
};

#endif

//! \}
