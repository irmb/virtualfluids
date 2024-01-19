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
//! \addtogroup parallel
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include "NullCommunicator.h"

#include <memory>

namespace vf::parallel
{

std::shared_ptr<Communicator> NullCommunicator::getInstance()
{
    std::lock_guard<std::mutex> myLock(instantiation_mutex);
    if (!instance) {
        instance = std::make_shared<NullCommunicator>();
    }
    return instance;
}
//////////////////////////////////////////////////////////////////////////
double NullCommunicator::Wtime()
{
    return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getBundleID() const
{
    return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getNumberOfBundles() const
{
    return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getProcessID() const
{
    return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getNumberOfProcesses() const
{
    return 1;
}
//////////////////////////////////////////////////////////////////////////
void *NullCommunicator::getNativeCommunicator()
{
    return NULL;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getRoot() const
{
    return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getBundleRoot() const
{
    return 0;
}
//////////////////////////////////////////////////////////////////////////
int NullCommunicator::getProcessRoot() const
{
    return 0;
}
//////////////////////////////////////////////////////////////////////////
std::vector<std::string> NullCommunicator::gather(const std::string & /*str*/)
{
    return {};
}
//////////////////////////////////////////////////////////////////////////

void NullCommunicator::sendSerializedObject(std::stringstream &stream, int target)
{
}
//////////////////////////////////////////////////////////////////////////
void NullCommunicator::receiveSerializedObject(std::stringstream &stream, int source)
{
}

int NullCommunicator::getProcessID(int bundle, int rank) const
{
    return 0;
}
bool NullCommunicator::isRoot() const
{
    return true;
}

int NullCommunicator::getNumberOfProcessesInBundle(int bundle) const
{
    return 0;
}
void NullCommunicator::barrier()
{
}
void NullCommunicator::abort(int errorcode)
{
}

std::vector<int> NullCommunicator::gather(std::vector<int> &values)
{
    return {};
}
std::vector<float> NullCommunicator::gather(std::vector<float> &values)
{
    return {};
}
std::vector<double> NullCommunicator::gather(std::vector<double> &values)
{
    return {};
}
std::vector<unsigned long long> NullCommunicator::gather(std::vector<unsigned long long> &values)
{
    return {};
}

void NullCommunicator::allGather(std::vector<int> &svalues, std::vector<int> &rvalues)
{
}
void NullCommunicator::allGather(std::vector<float> &svalues, std::vector<float> &rvalues)
{
}
void NullCommunicator::allGather(std::vector<double> &svalues, std::vector<double> &rvalues)
{
}
void NullCommunicator::allGather(std::vector<unsigned long long> &svalues, std::vector<unsigned long long> &rvalues)
{
}

void NullCommunicator::broadcast(int &value)
{
}
void NullCommunicator::broadcast(float &value)
{
}
void NullCommunicator::broadcast(double &value)
{
}
void NullCommunicator::broadcast(long int &value)
{
}
void NullCommunicator::broadcast(std::vector<int> &values)
{
}
void NullCommunicator::broadcast(std::vector<float> &values)
{
}
void NullCommunicator::broadcast(std::vector<double> &values)
{
}
void NullCommunicator::broadcast(std::vector<long int> &values)
{
}

void NullCommunicator::receiveSend(uint *buffer_receive, int size_buffer_recv, int neighbor_rank_recv, uint *buffer_send,
                                   int size_buffer_send, int neighbor_rank_send) const
{
}

void NullCommunicator::send(real *sbuf, int count_s, int nb_rank) const {};
double NullCommunicator::reduceSum(double /*quantityPerProcess*/) const
{
    return 0.0;
};
int NullCommunicator::mapCudaDevicesOnHosts(const std::vector<unsigned int> &devices, int numberOfDevices) const
{
    return 0;
}

void NullCommunicator::receiveSend(real *buffer_send, int size_buffer_send, real *buffer_receive, int size_buffer_recv,
                                   int neighbor_rank) const
{
}

void NullCommunicator::receiveNonBlocking(real *rbuf, int count_r, int sourceRank)
{
}
void NullCommunicator::sendNonBlocking(real *sbuf, int count_s, int destinationRank)
{
}

void NullCommunicator::send(real *sbuf, int count_s, int destinationRank)
{
}

void NullCommunicator::waitAll()
{
}

void NullCommunicator::resetRequests()
{
}
} // namespace vf::parallel

//! \}
