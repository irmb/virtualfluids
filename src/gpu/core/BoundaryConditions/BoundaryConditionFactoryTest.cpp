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
#include <gmock/gmock.h>
#include <typeindex>

#include "BoundaryConditionFactory.h"
#include "gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"

#include "BoundaryConditions/Outflow/Outflow.h"
#include "BoundaryConditions/Pressure/Pressure.h"
#include "BoundaryConditions/NoSlip/NoSlip.h"
#include "BoundaryConditions/Velocity/Velocity.h"
#include "BoundaryConditions/Slip/Slip.h"
#include "BoundaryConditions/Stress/Stress.h"


using bcFunction = void (*)(LBMSimulationParameter *, QforBoundaryConditions *);
using bcFunctionParamter = void (*)(Parameter *, QforBoundaryConditions *, const int level);

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
    auto bc = bcFactory.getPressureBoundaryConditionPre();
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
    EXPECT_THROW(bc(nullptr, nullptr, 0), std::bad_function_call);
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

bcFunction getPressureBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getPressureBoundaryConditionPre();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, pressureBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumIncompressible);
    EXPECT_TRUE(*(getPressureBcTarget(bcFactory)) == PressureNonEquilibriumIncompressible)
        << "The returned boundary condition is not the expected function PressureNonEquilibriumIncompressible.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
    EXPECT_TRUE(*(getPressureBcTarget(bcFactory)) == PressureNonEquilibriumCompressible)
        << "The returned boundary condition is not the expected function PressureNonEquilibriumCompressible.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
    EXPECT_TRUE(*(getPressureBcTarget(bcFactory)) == OutflowNonReflecting)
        << "The returned boundary condition is not the expected function OutflowNonReflecting_Device.";
}

bcFunction getGeometryBcTarget(BoundaryConditionFactory &bcFactory)
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
    auto bcTarget = *bc.target<bcFunctionParamter>();
    EXPECT_TRUE(*bcTarget == StressBounceBackCompressible)
        << "The returned boundary condition is not the expected function StressBounceBackCompressible.";

    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressCompressible);
    bc = bcFactory.getStressBoundaryConditionPost();
    bcTarget = *bc.target<bcFunctionParamter>();
    EXPECT_TRUE(*bcTarget == StressCompressible)
        << "The returned boundary condition is not the expected function StressCompressible.";
}

//! \}
