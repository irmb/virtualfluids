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
//! \addtogroup cpu_Connectors_tests Connectors
//! \ingroup cpu_core_tests core
//! \{
#include <gmock/gmock.h>

#include <memory>

#include <basics/container/CbVector.h>

#include <parallel/transmitter/TbTransmitter.h>

#include "Block3D.h"
#include "CreateTransmittersHelper.h"
#include "FineToCoarseVectorConnector.h"

class FineToCoarseVectorConnectorTest : public testing::Test
{

    void SetUp() override
    {
        block = std::make_shared<Block3D>();
    }

public:
    CreateTransmittersHelper::TransmitterPtr senderFCevenEvenSW, receiverFCevenEvenSW;
    std::shared_ptr<Block3D> block;
};

TEST_F(FineToCoarseVectorConnectorTest, getLocalMinMax)
{
    using namespace vf::lbm::dir;

    int sendDir = dP00;
    block->setInterpolationFlagFC(sendDir);
    // FineToCoarseVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver,
    // int sendDir, InterpolationProcessorPtr iprocessor, CFconnectorType connType);
    InterpolationProcessorPtr iprocessor;
    auto sut = FineToCoarseVectorConnector<TbTransmitter<CbVector<real>>>(block, senderFCevenEvenSW, receiverFCevenEvenSW,
                                                                          sendDir, iprocessor, EvenOddNW);

    //(int &minX1, int &minX2, int &minX3, int &maxX1, int &maxX2, int &maxX3);
    // SPtr<DistributionArray3D> fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
    int maxX1 = 5; //(int)fFrom->getNX1();
    int maxX2 = 5; //(int)fFrom->getNX2();
    int maxX3 = 5; //(int)fFrom->getNX3();
    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;
    sut.getLocalMinMax(minX1, minX2, minX3, maxX1, maxX2, maxX3);

    int expectedMaxX1 = 2;
    EXPECT_THAT(maxX1, testing::Eq(expectedMaxX1));
}
//! \}
