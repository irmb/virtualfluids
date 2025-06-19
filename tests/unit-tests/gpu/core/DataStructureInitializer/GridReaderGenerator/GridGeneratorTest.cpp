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
//! \addtogroup gpu_DataStructureInitializer_tests DataStructureInitializer
//! \ingroup gpu_core_tests core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include <gmock/gmock.h>

#include <gpu/core/DataStructureInitializer/GridReaderGenerator/GridGenerator.h>
#include <gpu/core/DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h>

#include <basics/DataTypes.h>
#include "Calculation/Calculation.h"
#include "Cuda/CudaMemoryManager.h"

#include "Parameter/Parameter.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridImp.h"
#include "gpu/GridGenerator/utilities/communication.h"
#include "parallel/Communicator.h"

#include <parallel/NullCommunicator.h>

namespace grid_generator_test
{

class LevelGridBuilderStub : public LevelGridBuilder
{
private:
    const SPtr<Grid> grid;
    LevelGridBuilderStub() = default;

public:
    uint numberOfSendIndices = 0;

    explicit LevelGridBuilderStub(SPtr<Grid> grid) : grid(std::move(grid)){};

    uint getCommunicationProcess(int direction) override
    {
        uint process = 0;
        if (direction != communication_directions::MX)
            process = (uint)INVALID_INDEX;
        return process;
    }

    uint getNumberOfGridLevels() const override
    {
        return 2;
    }

    uint getNumberOfSendIndices(int direction, uint level) override
    {
        return numberOfSendIndices;
    }

    uint getNumberOfReceiveIndices(int direction, uint level) override
    {
        return 0;
    }

    void getSendIndices(uint *sendIndices, int direction, int level) override
    {
    }

    void getReceiveIndices(uint *recvIndices, int direction, int level) override
    {
    }
};

class CudaMemoryManagerDouble : public CudaMemoryManager
{
public:
    explicit CudaMemoryManagerDouble(std::shared_ptr<Parameter> parameter) : CudaMemoryManager(parameter){};

    void cudaAllocProcessNeighbor(ProcessNeighbor27&, ProcessNeighbor27&, ProcessNeighbor27&, ProcessNeighbor27&) override{};
    void cudaCopyProcessNeighborIndex(ProcessNeighbor27&, ProcessNeighbor27&, ProcessNeighbor27&, ProcessNeighbor27&) override{};
};

class IndexRearrangementForStreamsDouble : public IndexRearrangementForStreams
{
public:
    IndexRearrangementForStreamsDouble(std::shared_ptr<Parameter> para, std::shared_ptr<GridBuilder> builder,
                                       vf::parallel::Communicator &communicator)
        : IndexRearrangementForStreams(para, builder, communicator){};

    void initCommunicationArraysForCommAfterFinetoCoarse(ProcessNeighbor27& /**/, ProcessNeighbor27& /**/, ProcessNeighbor27& /**/, ProcessNeighbor27& /**/, ProcessNeighbor27& /**/, ProcessNeighbor27& /**/, ProcessNeighbor27& /**/, ProcessNeighbor27& /**/, int/**/, int/**/) const override{};
};

} // namespace grid_generator_test

using namespace grid_generator_test;

class GridGeneratorTests_initalValuesDomainDecompostion : public testing::Test
{
public:
    void act() const
    {
        gridGenerator->initalValuesDomainDecompostion();
    }

protected:
    SPtr<Parameter> para;
    std::shared_ptr<LevelGridBuilderStub> builder;

    const uint level = 1;
    const uint direction = communication_directions::MX;

    SPtr<GridGenerator> gridGenerator;

private:
    void SetUp() override
    {
        para = std::make_shared<Parameter>();
        para->setMaxLevel(level + 1); // setMaxLevel resizes parH and parD
        for (uint i = 0; i <= level; i++) {
            para->parH[i] = std::make_shared<LBMSimulationParameter>();
            para->parD[i] = std::make_shared<LBMSimulationParameter>();
        }
        para->setNumprocs(2);

        builder = std::make_shared<LevelGridBuilderStub>(nullptr);
        auto communicator = vf::parallel::NullCommunicator::getInstance();

        gridGenerator = std::make_shared<GridGenerator>(builder, para, std::make_shared<CudaMemoryManagerDouble>(para), *communicator);
        gridGenerator->setIndexRearrangementForStreams(
            std::make_unique<IndexRearrangementForStreamsDouble>(para, builder, *communicator));
    }
};

TEST_F(GridGeneratorTests_initalValuesDomainDecompostion, whenNoCommunication_sendProcessNeighborShouldNotExist)
{
    act();
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborsX.size(), testing::Eq(0));
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborsY.size(), testing::Eq(0));
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborsZ.size(), testing::Eq(0));
}

TEST_F(GridGeneratorTests_initalValuesDomainDecompostion, whenCommunicationInX_sendProcessNeighborShouldExistInX)
{
    builder->numberOfSendIndices = 1;
    act();
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborsX.size(),
                testing::Eq(1)); // one entry for communication_directions::MX
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborsY.size(), testing::Eq(0));
    EXPECT_THAT(para->getParH(level)->sendProcessNeighborsZ.size(), testing::Eq(0));
}

//! \}
