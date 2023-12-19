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
//! \addtogroup transmitter
//! \ingroup parallel
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef PARALLEL_TOTRANSMITTERLOCAL_H
#define PARALLEL_TOTRANSMITTERLOCAL_H

#include <basics/PointerDefinitions.h>
#include <basics/utilities/UbException.h>

#include "parallel/transmitter/TbTransmitter.h"

//////////////////////////////////////////////////////////////////////////
// LocalTransmitter lokalen Datenaustausch
// data = send- und zugleich receive-buffer
template <typename T>
class TbLocalTransmitter : public TbTransmitter<T>
{
public:
    using TbLocalTransmitterPtr = SPtr<TbLocalTransmitter<T>>;

    using value_type = T;

public:
    TbLocalTransmitter() : TbTransmitter<T>() {}

    bool isLocalTransmitter() const override { return true; }
    bool isRemoteTransmitter() const override { return !this->isLocalTransmitter(); }

    // send buffer wird autom resized
    void sendDataSize() override {}
    // reiceive braucht nichts machen, da send==receive buffer ;-)
    void receiveDataSize() override {}

    void sendData() override {}
    value_type &receiveData() override { return this->data; }

    std::string toString() const override { return "TbLocalTransmitter" + (std::string) typeid(T).name(); }
};

//////////////////////////////////////////////////////////////////////////
// TbVectorSender/ReceiverLocal lokalen Datenaustausch ueber ZWEI vektoren
template <typename T>
class TbVectorReceiverLocal : public TbTransmitter<T>
{
public:
    using value_type = T;

public:
    TbVectorReceiverLocal() : TbTransmitter<value_type>() {}
    // virtual ~TbVectorReceiverLocal() { std::cout<<typeid(*this).name()<<" tot"<<std::endl;   }

    bool isLocalTransmitter() const { return true; }
    bool isRemoteTransmitter() const { return !this->isLocalTransmitter(); }

    // send buffer wird autom resized
    void sendDataSize() { UB_THROW(UbException(UB_EXARGS, "empfaengt nur")); }
    // reiceive braucht nichts machen, das macht der sender :-)
    void receiveDataSize() {}

    void sendData() { UB_THROW(UbException(UB_EXARGS, "empfaengt nur")); }
    value_type &receiveData() { return this->data; }

    std::string toString() const { return "TbVectorReceiverLocal<" + (std::string) typeid(T).name() + ">"; }
};

template <typename T>
class TbVectorSenderLocal : public TbTransmitter<T>
{
public:
    using value_type = T;

public:
    TbVectorSenderLocal(SPtr<TbVectorReceiverLocal<value_type>> receiver)
        : TbTransmitter<value_type>(), receiver(receiver)
    {
    }
    // virtual ~TbVectorSenderLocal() { std::cout<<typeid(*this).name()<<" tot"<<std::endl;   }

    bool isLocalTransmitter() const { return true; }
    bool isRemoteTransmitter() const { return !this->isLocalTransmitter(); }

    // send buffer wird autom resized
    void sendDataSize()
    {
        assert(receiver != NULL);
        receiver->getData().resize(this->data.size());
    }
    // reiceive braucht nichts machen, da send==receive buffer ;-)
    void receiveDataSize() { UB_THROW(UbException(UB_EXARGS, "sendet nur")); }

    void sendData()
    {
        assert(this->data.size() == receiver->getData().size());
        receiver->getData() = this->data;
        //       for(int i=(int)this->data.size()-1; i>=0; --i)
        //          receiver->getData()[i]= this->data[i];
    }
    value_type &receiveData() { UB_THROW(UbException(UB_EXARGS, "sendet nur")); }

    std::string toString() const { return "TbVectorSenderLocal<" + (std::string) typeid(T).name() + ">"; }

protected:
    SPtr<TbVectorReceiverLocal<value_type>> receiver;
};

#endif // TOTRANSMITTERLOCAL_H

//! \}
