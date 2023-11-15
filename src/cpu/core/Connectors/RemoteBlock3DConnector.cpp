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
//! \file RemoteBlock3DConnector.cpp
//! \ingroup Connectors
//! \author Konstantin Kutscher
//=======================================================================================
#include "RemoteBlock3DConnector.h"

//////////////////////////////////////////////////////////////////////////
RemoteBlock3DConnector::RemoteBlock3DConnector(SPtr<Block3D> block, VectorTransmitterPtr sender,
                                               VectorTransmitterPtr receiver, int sendDir)
    : Block3DConnector(sendDir), block(block), sender(sender), receiver(receiver)
{
    if (!block || !sender || !receiver)
        UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));
}
//////////////////////////////////////////////////////////////////////////

bool RemoteBlock3DConnector::isLocalConnector() { return !this->isRemoteConnector(); }
//////////////////////////////////////////////////////////////////////////

bool RemoteBlock3DConnector::isRemoteConnector()
{
    return ((sender && sender->isRemoteTransmitter()) || (receiver && receiver->isRemoteTransmitter()));
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::sendTransmitterDataSize()
{
    assert(sender != NULL);
    sender->sendDataSize();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::receiveTransmitterDataSize()
{
    assert(receiver != NULL);
    receiver->receiveDataSize();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::prepareForSend()
{
    assert(sender != NULL);
    sender->prepareForSend();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::sendVectors()
{
    assert(sender != NULL);
    sender->sendData();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::prepareForReceive()
{
    assert(receiver != NULL);
    receiver->prepareForReceive();
}
//////////////////////////////////////////////////////////////////////////

void RemoteBlock3DConnector::receiveVectors()
{
    assert(receiver != NULL);
    receiver->receiveData();
}
