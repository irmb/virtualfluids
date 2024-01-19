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
//! \addtogroup cpu_Connectors Connectors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef RemoteBlock3DConnector_H
#define RemoteBlock3DConnector_H

#include <vector>

#include "Block3D.h"
#include "Block3DConnector.h"
#include "D3Q27System.h"
#include "LBMKernel.h"
#include "TransmitterType.h"

//! \brief The same level connector between two distributed blocks.
//! \details Data are copied into a vector (this is located in the transmitter) the vector is transmitted via transmitter.
class RemoteBlock3DConnector : public Block3DConnector
{
public:
    RemoteBlock3DConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver,
                           int sendDir);

    bool isLocalConnector() override;
    bool isRemoteConnector() override;

    void init() override = 0;

    void sendTransmitterDataSize() override;
    void receiveTransmitterDataSize() override;

    void prepareForSend() override;
    void sendVectors() override;

    void prepareForReceive() override;
    void receiveVectors() override;

    void fillSendVectors() override          = 0;
    void distributeReceiveVectors() override = 0;

    bool isInterpolationConnectorCF() override { return false; }
    bool isInterpolationConnectorFC() override { return false; }

    real getSendRecieveTime() { return 0; }

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
    WPtr<Block3D> block;
    VectorTransmitterPtr sender;
    VectorTransmitterPtr receiver;
};

#endif // RemoteBlock3DConnector_H

//! \}
