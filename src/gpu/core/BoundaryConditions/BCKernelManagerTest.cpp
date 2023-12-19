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

#include <gmock/gmock.h>
#include <stdexcept>

#include "BoundaryConditionKernelManager.h"
#include "BoundaryConditions/BoundaryConditionFactory.h"
#include "Parameter/Parameter.h"
#include <basics/PointerDefinitions.h>

class BoundaryConditionKernelManagerTest_BCsNotSpecified : public testing::Test
{
protected:
    BoundaryConditionFactory bcFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();

    void SetUp() override
    {
        para->initLBMSimulationParameter();
    }
};

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, velocityBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, velocityBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, noSlipBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->noSlipBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, noSlipBoundaryConditionPostNotSpecified_withBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->noSlipBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory)); // no throw, as a default is specified
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, slipBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->slipBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, slipBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->slipBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, pressureBoundaryConditionPreNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->pressureBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, pressureBoundaryConditionPreNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->pressureBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, geometryBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->geometryBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, geometryBoundaryConditionPostNotSpecified_withBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->geometryBC.numberOfBCnodes = 1;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory)); // no throw, as a default is specified
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, stressBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->stressBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, stressBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->stressBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, precursorBoundaryConditionPostNotSpecified_noBoundaryNodes_doesNotThrow)
{
    para->getParD(0)->precursorBC.numberOfBCnodes = 0;
    EXPECT_NO_THROW(BoundaryConditionKernelManager(para, &bcFactory));
}

TEST_F(BoundaryConditionKernelManagerTest_BCsNotSpecified, precursorBoundaryConditionPostNotSpecified_withBoundaryNodes_throws)
{
    para->getParD(0)->precursorBC.numberOfBCnodes = 1;
    EXPECT_THROW(BoundaryConditionKernelManager(para, &bcFactory), std::runtime_error);
}

class BoundaryConditionFactoryMock : public BoundaryConditionFactory
{
public:
    mutable uint numberOfCalls = 0;

    [[nodiscard]] boundaryCondition getVelocityBoundaryConditionPost(bool) const override
    {
        return [this](LBMSimulationParameter *, QforBoundaryConditions *) { numberOfCalls++; };
    }
};

class BoundaryConditionKernelManagerTest_runBCs : public testing::Test
{
protected:
    BoundaryConditionFactoryMock bcFactory;
    SPtr<Parameter> para = std::make_shared<Parameter>();
    UPtr<BoundaryConditionKernelManager> sut;

    void SetUp() override
    {
        para->initLBMSimulationParameter();
        sut = std::make_unique<BoundaryConditionKernelManager>(para, &bcFactory);
    }
};

TEST_F(BoundaryConditionKernelManagerTest_runBCs, runVelocityBCKernelPost)
{
    para->getParD(0)->velocityBC.numberOfBCnodes = 1;
    sut->runVelocityBCKernelPost(0);
    EXPECT_THAT(bcFactory.numberOfCalls, testing::Eq(1));

    bcFactory.numberOfCalls = 0;
    para->getParD(0)->velocityBC.numberOfBCnodes = 0;
    sut->runVelocityBCKernelPost(0);
    EXPECT_THAT(bcFactory.numberOfCalls, testing::Eq(0));
}

//! \}
