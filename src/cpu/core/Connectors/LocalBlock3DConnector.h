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
//! \file LocalBlock3DConnector.h
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef LocalBlock3DConnector_H
#define LocalBlock3DConnector_H

#include "Block3D.h"
#include "Block3DConnector.h"

#include <basics/PointerDefinitions.h>
#include <lbm/constants/D3Q27.h>

//! A class provides an interface for connectors in shared memory
class LocalBlock3DConnector : public Block3DConnector
{
public:
    LocalBlock3DConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir)
        : Block3DConnector(sendDir), from(from), to(to)
    {
    }
    ~LocalBlock3DConnector() override = default;
    void sendTransmitterDataSize() override {}
    void receiveTransmitterDataSize() override {}
    void init() override = 0;
    void prepareForReceive() override {}
    void prepareForSend() override {}
    void fillSendVectors() override {}
    void sendVectors() override = 0;
    void receiveVectors() override {}

    void distributeReceiveVectors() override {}

    bool isLocalConnector() override { return true; }
    bool isRemoteConnector() override { return false; }
    bool isInterpolationConnectorCF() override { return false; }
    bool isInterpolationConnectorFC() override { return false; }

    real getSendRecieveTime();

    void prepareForSendX1() override {}
    void prepareForSendX2() override {}
    void prepareForSendX3() override {}

    void sendVectorsX1() override {}
    void sendVectorsX2() override {}
    void sendVectorsX3() override {}

    void prepareForReceiveX1() override {}
    void prepareForReceiveX2() override {}
    void prepareForReceiveX3() override {}

    void receiveVectorsX1() override {}
    void receiveVectorsX2() override {}
    void receiveVectorsX3() override {}

protected:
    WPtr<Block3D> from;
    WPtr<Block3D> to;
};

#endif // LocalBlock3DConnector_H
