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
//! \addtogroup gpu_BoundaryConditions_tests BoundaryConditions
//! \ingroup gpu_core_tests
//! \{
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#include "Calculation/Calculation.h"
#include <gmock/gmock.h>

#include <typeindex>

#include <gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h>

#include <gpu/core/BoundaryConditions/BoundaryConditionFactory.h>
#include <gpu/core/BoundaryConditions/NoSlip/NoSlip.h>
#include <gpu/core/BoundaryConditions/Outflow/Outflow.h>
#include <gpu/core/BoundaryConditions/Pressure/Pressure.h>
#include <gpu/core/BoundaryConditions/Slip/Slip.h>
#include <gpu/core/BoundaryConditions/Stress/Stress.h>
#include <gpu/core/BoundaryConditions/Stress/SurfaceLayer.h>
#include <gpu/core/BoundaryConditions/Velocity/Velocity.h>
#include <gpu/core/BoundaryConditions/AdvectionDiffusion/AdvectionDiffusion.h>

using bcFunction = void (*)(LBMSimulationParameter *, QforBoundaryConditions *);
using bcFunctionDirectional = void (*)(LBMSimulationParameter *, QforDirectionalBoundaryCondition *);
using bcFunctionParameter = void (*)(Parameter *, QforBoundaryConditions *, const int level);
using adBCNoFluxFunction = void (*)(LBMSimulationParameter *, AdvectionDiffusionNoFluxBoundaryConditions);

// tests for default boundary conditions
TEST(BoundaryConditionFactoryTest, defaultVelocityBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getVelocityBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultNoSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getNoSlipBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw
}

TEST(BoundaryConditionFactoryTest, defaultSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getSlipBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultPressureBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = std::get<BoundaryConditionKernel>(bcFactory.getPressureBoundaryConditionPre());
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultGeometryBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getGeometryBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw
}

TEST(BoundaryConditionFactoryTest, defaultStressBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getStressBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultADNoFluxBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getAdvectionDiffusionNoFluxBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, AdvectionDiffusionNoFluxBoundaryConditions{})); // empty lambda function should not throw
}


TEST(BoundaryConditionFactoryTest, defaultADFluxBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getAdvectionDiffusionFluxBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, AdvectionDiffusionFluxBoundaryConditions{}), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultADDirichletBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getAdvectionDiffusionDirichletBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, AdvectionDiffusionDirichletBoundaryConditions{}), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultADNeumannBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getAdvectionDiffusionNeumannBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, AdvectionDiffusionNeumannBoundaryConditions{}), std::bad_function_call);
}

// tests for boundary conditions which are set by the user (tests both set and get functions)

bcFunction getVelocityBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getVelocityBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, velocityBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityBounceBack);
    EXPECT_TRUE(*(getVelocityBcTarget(bcFactory)) == VelocityBounceBack)
        << "The returned boundary condition is not the expected function VelocityBounceBack.";

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedIncompressible);
    EXPECT_TRUE(*(getVelocityBcTarget(bcFactory)) == VelocityInterpolatedIncompressible)
        << "The returned boundary condition is not the expected function VelocityInterpolatedIncompressible.";

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);
    EXPECT_TRUE(*(getVelocityBcTarget(bcFactory)) == VelocityInterpolatedCompressible)
        << "The returned boundary condition is not the expected function VelocityInterpolatedCompressible.";

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityWithPressureInterpolatedCompressible);
    EXPECT_TRUE(*(getVelocityBcTarget(bcFactory)) == VelocityWithPressureInterpolatedCompressible)
        << "The returned boundary condition is not the expected function VelocityWithPressureInterpolatedCompressible.";
}

bcFunction getNoSlipBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getNoSlipBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, noSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipDelayBounceBack);
    auto bc = bcFactory.getNoSlipBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipBounceBack);
    EXPECT_TRUE( *(getNoSlipBcTarget(bcFactory)) == NoSlipBounceBack)
        << "The returned boundary condition is not the expected function NoSlipBounceBack.";

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipInterpolatedIncompressible);
    EXPECT_TRUE( *(getNoSlipBcTarget(bcFactory)) == NoSlipInterpolatedIncompressible)
        << "The returned boundary condition is not the expected function NoSlipInterpolatedIncompressible.";

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipInterpolatedCompressible);
    EXPECT_TRUE( *(getNoSlipBcTarget(bcFactory)) == NoSlipInterpolatedCompressible)
        << "The returned boundary condition is not the expected function NoSlipInterpolatedCompressible.";
}

bcFunction getSlipBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getSlipBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, slipBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);
    EXPECT_TRUE( *(getSlipBcTarget(bcFactory)) == SlipCompressible)
        << "The returned boundary condition is not the expected function SlipCompressible.";

    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipTurbulentViscosityCompressible);
    EXPECT_TRUE( *(getSlipBcTarget(bcFactory)) == SlipTurbulentViscosityCompressible)
        << "The returned boundary condition is not the expected function QSlipDevCompTurbulentViscosity27.";
}


bcFunction getSurfaceLayerBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getSurfaceLayerBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, surfaceLayerBoundaryCondition)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setSurfaceLayerBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackCompressible, BoundaryConditionFactory::SurfaceLayerBC::SurfaceHeatFlux);
    auto bc = bcFactory.getSurfaceLayerBoundaryConditionPost();
    auto bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == SurfaceLayerBounceBackCompressibleHeatFlux)
        << "The returned boundary condition is not the expected function SurfaceLayerBounceBackCompressibleHeatFlux.";

    bcFactory.setSurfaceLayerBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackWithPressureCompressible, BoundaryConditionFactory::SurfaceLayerBC::SurfaceHeatFlux);
    bc = bcFactory.getSurfaceLayerBoundaryConditionPost();
    bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == SurfaceLayerBounceBackWithPressureCompressibleHeatFlux)
        << "The returned boundary condition is not the expected function SurfaceLayerBounceBackWithPressureCompressibleHeatFlux.";

    bcFactory.setSurfaceLayerBoundaryCondition(BoundaryConditionFactory::StressBC::StressInterpolatedCompressible, BoundaryConditionFactory::SurfaceLayerBC::SurfaceHeatFlux);
    bc = bcFactory.getSurfaceLayerBoundaryConditionPost();
    bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == SurfaceLayerInterpolatedCompressibleHeatFlux)
        << "The returned boundary condition is not the expected function SurfaceLayerInterpolatedCompressibleHeatFlux.";

    bcFactory.setSurfaceLayerBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackCompressible, BoundaryConditionFactory::SurfaceLayerBC::SurfaceTemperature);
    bc = bcFactory.getSurfaceLayerBoundaryConditionPost();
    bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == SurfaceLayerBounceBackCompressibleSurfaceTemperature)
        << "The returned boundary condition is not the expected function SurfaceLayerBounceBackCompressibleSurfaceTemperature.";

    bcFactory.setSurfaceLayerBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackWithPressureCompressible, BoundaryConditionFactory::SurfaceLayerBC::SurfaceTemperature);
    bc = bcFactory.getSurfaceLayerBoundaryConditionPost();
    bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == SurfaceLayerBounceBackWithPressureCompressibleSurfaceTemperature)
        << "The returned boundary condition is not the expected function SurfaceLayerBounceBackWithPressureCompressibleSurfaceTemperature.";

    bcFactory.setSurfaceLayerBoundaryCondition(BoundaryConditionFactory::StressBC::StressInterpolatedCompressible, BoundaryConditionFactory::SurfaceLayerBC::SurfaceTemperature);
    bc = bcFactory.getSurfaceLayerBoundaryConditionPost();
    bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == SurfaceLayerInterpolatedCompressibleSurfaceTemperature)
        << "The returned boundary condition is not the expected function SurfaceLayerInterpolatedCompressibleSurfaceTemperature.";
}


bcFunctionDirectional getDirectionalPressureBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = std::get<DirectionalBoundaryConditionKernel>(bcFactory.getPressureBoundaryConditionPre());
    void (*bcTarget)(LBMSimulationParameter*, QforDirectionalBoundaryCondition*) =
        (*bc.target<void (*)(LBMSimulationParameter*, QforDirectionalBoundaryCondition*)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, pressureBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumIncompressible);
    EXPECT_TRUE(*(getDirectionalPressureBcTarget(bcFactory)) == PressureNonEquilibriumIncompressible)
        << "The returned boundary condition is not the expected function PressureNonEquilibriumIncompressible.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
    EXPECT_TRUE(*(getDirectionalPressureBcTarget(bcFactory)) == PressureNonEquilibriumCompressible)
        << "The returned boundary condition is not the expected function PressureNonEquilibriumCompressible.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
    EXPECT_TRUE(*(getDirectionalPressureBcTarget(bcFactory)) == OutflowNonReflecting)
        << "The returned boundary condition is not the expected function OutflowNonReflecting.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflectivePressureCorrection);
    EXPECT_TRUE(*(getDirectionalPressureBcTarget(bcFactory)) == OutflowNonReflectingPressureCorrection)
        << "The returned boundary condition is not the expected function OutflowNonReflectingPressureCorrection.";
}

bcFunction getGeometryBcTarget(BoundaryConditionFactory& bcFactory)
{
    auto bc = bcFactory.getGeometryBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, geometryBC)
{
    auto bcFactory = BoundaryConditionFactory();

    // velocity
    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedIncompressible);
    EXPECT_TRUE( *(getGeometryBcTarget(bcFactory)) == VelocityInterpolatedIncompressible)
        << "The returned boundary condition is not the expected function VelocityInterpolatedIncompressible.";

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityInterpolatedCompressible);
    EXPECT_TRUE( *(getGeometryBcTarget(bcFactory)) == VelocityInterpolatedCompressible)
        << "The returned boundary condition is not the expected function VelocityInterpolatedCompressible.";

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityWithPressureInterpolatedCompressible);
    EXPECT_TRUE(*(getGeometryBcTarget(bcFactory)) == VelocityWithPressureInterpolatedCompressible)
        << "The returned boundary condition is not the expected function VelocityWithPressureInterpolatedCompressible.";

    // no slip
    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipDelayBounceBack);
    auto bc = bcFactory.getGeometryBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipInterpolatedIncompressible);
    EXPECT_TRUE( *(getGeometryBcTarget(bcFactory)) == NoSlipInterpolatedIncompressible)
        << "The returned boundary condition is not the expected function NoSlipInterpolatedIncompressible.";

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipInterpolatedCompressible);
    EXPECT_TRUE( *(getGeometryBcTarget(bcFactory)) == NoSlipInterpolatedCompressible)
        << "The returned boundary condition is not the expected function NoSlipInterpolatedCompressible.";
}

TEST(BoundaryConditionFactoryTest, stressBoundaryConditions)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackCompressible);
    auto bc = bcFactory.getStressBoundaryConditionPost();
    auto bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == StressBounceBackCompressible)
        << "The returned boundary condition is not the expected function StressBounceBackCompressible.";

    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBackWithPressureCompressible);
    bc = bcFactory.getStressBoundaryConditionPost();
    bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == StressBounceBackWithPressureCompressible)
        << "The returned boundary condition is not the expected function StressBounceBackWithPressureCompressible.";

    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressInterpolatedCompressible);
    bc = bcFactory.getStressBoundaryConditionPost();
    bcTarget = *bc.target<bcFunction>();
    EXPECT_TRUE(*bcTarget == StressInterpolatedCompressible)
        << "The returned boundary condition is not the expected function StressInterpolatedCompressible.";
}

TEST(BoundaryConditionFactoryTest, hasDirectionalPressureBoundaryCondition_whenDirectionalBC_returnsTrue)
{
    auto bcFactory = BoundaryConditionFactory();

    // check all directional pressure boundary conditions
    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
    EXPECT_TRUE(bcFactory.hasDirectionalPressureBoundaryCondition());

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumIncompressible);
    EXPECT_TRUE(bcFactory.hasDirectionalPressureBoundaryCondition());

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
    EXPECT_TRUE(bcFactory.hasDirectionalPressureBoundaryCondition());

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflectivePressureCorrection);
    EXPECT_TRUE(bcFactory.hasDirectionalPressureBoundaryCondition());
}

TEST(BoundaryConditionFactoryTest, hasDirectionalPressureBoundaryCondition_noPressureBC_returnsTrue)
{
    auto bcFactory=BoundaryConditionFactory();
    EXPECT_FALSE(bcFactory.hasDirectionalPressureBoundaryCondition());
}

adBCNoFluxFunction getADNoFluxBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getAdvectionDiffusionNoFluxBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, AdvectionDiffusionNoFluxBoundaryConditions) = (*bc.target<adBCNoFluxFunction>());
    return bcTarget;
}


TEST(BoundaryConditionFactoryTest, ADNoFluxBoundaryConditions)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setAdvectionDiffusionNoFluxBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionNoFluxBC::NoFluxDelayedBounceBack);
    auto bc = bcFactory.getAdvectionDiffusionNoFluxBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, AdvectionDiffusionNoFluxBoundaryConditions{})); // empty lambda function should not throw
    
    bcFactory.setAdvectionDiffusionNoFluxBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionNoFluxBC::NoFluxBounceBack);
    auto bcTarget = getADNoFluxBcTarget(bcFactory);
    EXPECT_TRUE(bcTarget == AdvectionDiffusionNoFluxBounceBack)
        << "The returned boundary condition is not the expected function AdvectionDiffusionNoFluxBounceBack.";
}

auto getADFluxBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getAdvectionDiffusionFluxBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, AdvectionDiffusionFluxBoundaryConditions) =
        (*bc.target<void (*)(LBMSimulationParameter *, AdvectionDiffusionFluxBoundaryConditions)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, ADFluxBoundaryConditions)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setAdvectionDiffusionFluxBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxBounceBack);
    auto bc = getADFluxBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionFluxBounceBack)
        << "The returned boundary condition is not the expected function AdvectionDiffusionFluxBounceBack.";

    bcFactory.setAdvectionDiffusionFluxBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxCompressible);
    bc = getADFluxBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionFluxCompressible)
        << "The returned boundary condition is not the expected function AdvectionDiffusionFluxCompressible.";

    bcFactory.setAdvectionDiffusionFluxBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionFluxBC::FluxTurbulentViscosityCompressible);
    bc = getADFluxBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionFluxTurbulentViscosityCompressible)
        << "The returned boundary condition is not the expected function AdvectionDiffusionFluxTurbulentViscosityCompressible.";
}

auto getADDirichletBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getAdvectionDiffusionDirichletBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, AdvectionDiffusionDirichletBoundaryConditions) =
        (*bc.target<void (*)(LBMSimulationParameter *, AdvectionDiffusionDirichletBoundaryConditions)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, ADDirichletBoundaryConditions)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setAdvectionDiffusionDirichletBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackNoSlip);
    auto bc = getADDirichletBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionDirichletAntiBounceBackNoSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionDirichletAntiBounceBackNoSlip.";

    bcFactory.setAdvectionDiffusionDirichletBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletAntiBounceBackSlip);
    bc = getADDirichletBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionDirichletAntiBounceBackSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionDirichletAntiBounceBackSlip.";

    bcFactory.setAdvectionDiffusionDirichletBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletInterpolatedNoSlip);
    bc = getADDirichletBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionDirichletInterpolatedNoSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionDirichletInterpolatedNoSlip.";

    bcFactory.setAdvectionDiffusionDirichletBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionDirichletBC::DirichletInterpolatedSlip);
    bc = getADDirichletBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionDirichletInterpolatedSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionDirichletInterpolatedSlip.";
}

auto getADNeumannBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getAdvectionDiffusionNeumannBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, AdvectionDiffusionNeumannBoundaryConditions) =
        (*bc.target<void (*)(LBMSimulationParameter *, AdvectionDiffusionNeumannBoundaryConditions)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, ADNeumannBoundaryConditions)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setAdvectionDiffusionNeumannBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackNoSlip);
    auto bc = getADNeumannBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionNeumannAntiBounceBackNoSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionNeumannAntiBounceBackNoSlip.";

    bcFactory.setAdvectionDiffusionNeumannBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannAntiBounceBackSlip);
    bc = getADNeumannBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionNeumannAntiBounceBackSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionNeumannAntiBounceBackSlip.";

    bcFactory.setAdvectionDiffusionNeumannBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannInterpolatedNoSlip);
    bc = getADNeumannBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionNeumannInterpolatedNoSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionNeumannInterpolatedNoSlip.";

    bcFactory.setAdvectionDiffusionNeumannBoundaryCondition(BoundaryConditionFactory::AdvectionDiffusionNeumannBC::NeumannInterpolatedSlip);
    bc = getADNeumannBcTarget(bcFactory);
    EXPECT_TRUE(bc == AdvectionDiffusionNeumannInterpolatedSlip)
        << "The returned boundary condition is not the expected function AdvectionDiffusionNeumannInterpolatedSlip.";
}

//! \}
