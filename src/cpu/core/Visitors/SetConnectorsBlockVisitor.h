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
//! \addtogroup cpu_Visitors Visitors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef SETCONNECTORSBLOCKVISITOR_H
#define SETCONNECTORSBLOCKVISITOR_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "D3Q27System.h"
#include "Grid3D.h"
#include "CreateTransmittersHelper.h"
#include <parallel/Communicator.h>
#include "OneDistributionFullDirectConnector.h"
#include "OneDistributionFullVectorConnector.h"
#include "TwoDistributionsFullDirectConnector.h"
#include "TwoDistributionsFullVectorConnector.h"
#include "TwoDistributionsDoubleGhostLayerFullDirectConnector.h"
#include "TwoDistributionsDoubleGhostLayerFullVectorConnector.h"
#include "ThreeDistributionsFullDirectConnector.h"
#include "ThreeDistributionsFullVectorConnector.h"
#include "ThreeDistributionsDoubleGhostLayerFullDirectConnector.h"
#include "ThreeDistributionsDoubleGhostLayerFullVectorConnector.h"

#include <parallel/transmitter/TbTransmitterLocal.h>

//! \brief  A class sets connectors between blocks.
template <class T1, class T2>
class SetConnectorsBlockVisitor : public Block3DVisitor
{
public:
    using LocalConnector  = T1;
    using RemoteConnector = T2;
public:
    SetConnectorsBlockVisitor(std::shared_ptr<vf::parallel::Communicator> comm);
    ~SetConnectorsBlockVisitor() override;
    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
    //////////////////////////////////////////////////////////////////////////
protected:
    void setSameLevelConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block);
    void setRemoteConnectors(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir);
    std::shared_ptr<vf::parallel::Communicator> comm;
    int gridRank{0};
};

template <class T1, class T2>
SetConnectorsBlockVisitor<T1, T2>::SetConnectorsBlockVisitor(std::shared_ptr<vf::parallel::Communicator> comm)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), comm(comm)
{
}
//////////////////////////////////////////////////////////////////////////
template <class T1, class T2>
SetConnectorsBlockVisitor<T1, T2>::~SetConnectorsBlockVisitor(void)
{
}
//////////////////////////////////////////////////////////////////////////
template <class T1, class T2>
void SetConnectorsBlockVisitor<T1, T2>::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    if (!block)
        return;

    UBLOG(logDEBUG5, "SetConnectorsBlockVisitor::visit() - start");
    UBLOG(logDEBUG5, block->toString());

    gridRank = comm->getProcessID();
    grid->setRank(gridRank);

    setSameLevelConnectors(grid, block);

    UBLOG(logDEBUG5, "SetConnectorsBlockVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
template <class T1, class T2>
void SetConnectorsBlockVisitor<T1, T2>::setSameLevelConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;

    UBLOG(logDEBUG5, "SetConnectorsBlockVisitor::setSameLevelConnectors() - start");
    int blockRank = block->getRank();
    if (gridRank == blockRank && block->isActive()) {
        block->clearWeight();
        std::vector<SPtr<Block3D>> neighbors;
        int ix1   = block->getX1();
        int ix2   = block->getX2();
        int ix3   = block->getX3();
        int level = block->getLevel();

        for (int dir = D3Q27System::FSTARTDIR; dir <= D3Q27System::FENDDIR; dir++) { 
            SPtr<Block3D> neighBlock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);

            if (neighBlock) {
                int neighBlockRank = neighBlock->getRank();
                if (blockRank == neighBlockRank && neighBlock->isActive()) {
                    SPtr<Block3DConnector> connector;
                    connector = SPtr<Block3DConnector>(new LocalConnector(block, neighBlock, dir));
                    block->setConnector(connector);
                } else if (blockRank != neighBlockRank && neighBlock->isActive()) {
                    setRemoteConnectors(block, neighBlock, dir);

                    if (dir >= (int)dP00 && dir <= (int)d00M) {
                        int weight = block->getWeight(neighBlockRank);
                        weight++;
                        block->setWeight(neighBlockRank, weight);
                    }
                }
            }
        }

        int weight = block->getNumberOfLocalConnectorsForSurfaces();
        weight     = 6 - weight;
        block->addWeightForAll(weight);
    }
    UBLOG(logDEBUG5, "SetConnectorsBlockVisitor::setSameLevelConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
template <class T1, class T2>
void SetConnectorsBlockVisitor<T1, T2>::setRemoteConnectors(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir)
{
    UBLOG(logDEBUG5, "SetConnectorsBlockVisitor::setRemoteConnectors() - start");
    CreateTransmittersHelper helper;
    CreateTransmittersHelper::TransmitterPtr sender, receiver;
    helper.createTransmitters(sblock, tblock, dir, CreateTransmittersHelper::NONE, sender, receiver, comm,
                              CreateTransmittersHelper::MPI);

    SPtr<Block3DConnector> connector;
    connector = SPtr<Block3DConnector>(new RemoteConnector(sblock, sender, receiver, dir));
    sblock->setConnector(connector);
    UBLOG(logDEBUG5, "SetConnectorsBlockVisitor::setRemoteConnectors() - end");
}
//

using OneDistributionSetConnectorsBlockVisitor  = SetConnectorsBlockVisitor<OneDistributionFullDirectConnector, OneDistributionFullVectorConnector>;
using TwoDistributionsSetConnectorsBlockVisitor = SetConnectorsBlockVisitor<TwoDistributionsFullDirectConnector, TwoDistributionsFullVectorConnector>;
using TwoDistributionsDoubleGhostLayerSetConnectorsBlockVisitor = SetConnectorsBlockVisitor<TwoDistributionsDoubleGhostLayerFullDirectConnector, TwoDistributionsDoubleGhostLayerFullVectorConnector>;
using ThreeDistributionsSetConnectorsBlockVisitor = SetConnectorsBlockVisitor<ThreeDistributionsFullDirectConnector, ThreeDistributionsFullVectorConnector>;
using ThreeDistributionsDoubleGhostLayerSetConnectorsBlockVisitor = SetConnectorsBlockVisitor<ThreeDistributionsDoubleGhostLayerFullDirectConnector, ThreeDistributionsDoubleGhostLayerFullVectorConnector>;

#endif // SETCONNECTORSBLOCKVISITOR_H

//! \}
