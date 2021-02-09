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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file SetConnectorsBlockVisitor.h
//! \ingroup Visitors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef ConnectorBlockVisitor_H
#define ConnectorBlockVisitor_H

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "CreateTransmittersHelper.h"
#include "D3Q27System.h"

class Grid3D;
class Block3D;
class InterpolationProcessor;
class ConnectorFactory;

class ConnectorBlockVisitor : public Block3DVisitor
{
public:
    ConnectorBlockVisitor(SPtr<Communicator> comm, LBMReal nu, SPtr<InterpolationProcessor> iProcessor,
                          SPtr<ConnectorFactory> cFactory);
    ~ConnectorBlockVisitor() override;
    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;
    //////////////////////////////////////////////////////////////////////////
protected:
    void setSameLevelConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block);
    void setRemoteConnectors(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir);
    void setInterpolationConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block);
    void setInterpolationConnectors(SPtr<Block3D> fBlockSW, SPtr<Block3D> fBlockSE, SPtr<Block3D> fBlockNW,
                                    SPtr<Block3D> fBlockNE, SPtr<Block3D> cBlock, int dir);
    void createTransmitters(SPtr<Block3D> cBlock, SPtr<Block3D> fBlock, int dir, CreateTransmittersHelper::IBlock ib,
                            CreateTransmittersHelper::TransmitterPtr &senderCF,
                            CreateTransmittersHelper::TransmitterPtr &receiverCF,
                            CreateTransmittersHelper::TransmitterPtr &senderFC,
                            CreateTransmittersHelper::TransmitterPtr &receiverFC);
    SPtr<Communicator> comm;
    int gridRank;
    LBMReal nu;
    SPtr<InterpolationProcessor> iProcessor;
    SPtr<ConnectorFactory> cFactory;
};

#endif // ConnectorBlockVisitor_H
