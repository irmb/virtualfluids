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
//! \file TbTransmitter.h
//! \ingroup transmitter
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef TBTRANSMITTER_H
#define TBTRANSMITTER_H

#include <string>

/*================================================================================*/
/*  TbTransmitter                                                                 */
/*                                                                                */
/**
This Class provides the base for sending and receiving of data.
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 08.11.07
*/

/*
usage: ...
*/

//////////////////////////////////////////////////////////////////////////
// Transmitter
// macht nichts ausser daten senden und empfangen
template <typename T>
class TbTransmitter
{
public:
    using value_type = T;

public:
    TbTransmitter()          = default;
    virtual ~TbTransmitter() = default;

    virtual bool isLocalTransmitter() const  = 0;
    virtual bool isRemoteTransmitter() const = 0;

    // preprocess (e.g. synchronizing send-/receive-buffer)
    virtual void sendDataSize()    = 0;
    virtual void receiveDataSize() = 0;

    // calculation
    virtual void prepareForSend() {}
    virtual void sendData() = 0;
    virtual void prepareForReceive() {}
    virtual value_type &receiveData() = 0;

    // data-access
    inline value_type &getData() { return this->data; }
    inline const value_type &getData() const { return this->data; }

    // info-section (usable for remote transmitter)
    virtual int getSendToRank() const { return -1; }
    virtual int getSendToTag() const { return -1; }
    virtual int getRecvFromRank() const { return -1; }
    virtual int getRecvFromTag() const { return -1; }

    virtual std::string toString() const = 0;

protected:
    value_type data;
};

#endif // TBTRANSMITTER_H
