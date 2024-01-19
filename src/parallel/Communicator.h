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

#ifndef MPI_COMMUNICATOR_H
#define MPI_COMMUNICATOR_H

#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <basics/DataTypes.h>

namespace vf::parallel
{

//! \brief An abstract class for communication between processes in parallel computation
class Communicator
{
public:
    Communicator(const Communicator &) = delete;
    Communicator &operator=(const Communicator &rhs) = delete;
    static std::shared_ptr<Communicator> getInstance();

    virtual ~Communicator() = default;


    virtual double Wtime() = 0;
    virtual int getBundleID() const                      = 0;
    virtual int getNumberOfBundles() const               = 0;
    virtual int getProcessID() const                     = 0;
    virtual int getProcessID(int bundle, int rank) const = 0;
    virtual bool isRoot() const                          = 0;
    virtual void *getNativeCommunicator()                = 0;

    virtual void sendSerializedObject(std::stringstream &ss, int target)    = 0;
    virtual void receiveSerializedObject(std::stringstream &ss, int source) = 0;

    virtual int getRoot() const                                = 0;
    virtual int getBundleRoot() const                          = 0;
    virtual int getProcessRoot() const                         = 0;
    virtual int getNumberOfProcessesInBundle(int bundle) const = 0;
    virtual void barrier()                                     = 0;
    virtual void abort(int errorcode)                          = 0;

    virtual std::vector<std::string> gather(const std::string &str)                         = 0;
    virtual std::vector<int> gather(std::vector<int> &values)                               = 0;
    virtual std::vector<float> gather(std::vector<float> &values)                           = 0;
    virtual std::vector<double> gather(std::vector<double> &values)                         = 0;
    virtual std::vector<unsigned long long> gather(std::vector<unsigned long long> &values) = 0;

    virtual void allGather(std::vector<int> &svalues, std::vector<int> &rvalues)                               = 0;
    virtual void allGather(std::vector<float> &svalues, std::vector<float> &rvalues)                           = 0;
    virtual void allGather(std::vector<double> &svalues, std::vector<double> &rvalues)                         = 0;
    virtual void allGather(std::vector<unsigned long long> &svalues, std::vector<unsigned long long> &rvalues) = 0;

    virtual void broadcast(int &value)                    = 0;
    virtual void broadcast(float &value)                  = 0;
    virtual void broadcast(double &value)                 = 0;
    virtual void broadcast(long int &value)               = 0;
    virtual void broadcast(std::vector<int> &values)      = 0;
    virtual void broadcast(std::vector<float> &values)    = 0;
    virtual void broadcast(std::vector<double> &values)   = 0;
    virtual void broadcast(std::vector<long int> &values) = 0;

    virtual void receiveSend(uint *buffer_receive, int size_buffer_recv, int neighbor_rank_recv, uint *buffer_send,
                             int size_buffer_send, int neighbor_rank_send) const = 0;
    virtual int getNumberOfProcesses() const = 0;
    virtual void send(real *sbuf, int count_s, int nb_rank) const = 0;
    virtual double reduceSum(double quantityPerProcess) const = 0;
    virtual int mapCudaDevicesOnHosts(const std::vector<unsigned int> &devices, int numberOfDevices) const = 0;
    virtual void receiveSend(real *buffer_send, int size_buffer_send, real *buffer_receive, int size_buffer_recv,
                             int neighbor_rank) const = 0;
    virtual void receiveNonBlocking(real *rbuf, int count_r, int sourceRank) = 0;
    virtual void sendNonBlocking(real *sbuf, int count_s, int destinationRank) = 0;
    virtual void send(real *sbuf, int count_s, int destinationRank) = 0;
    virtual void waitAll() = 0;
    virtual void resetRequests() = 0;

protected:
    Communicator() = default;

    static std::mutex instantiation_mutex;

    static std::shared_ptr<Communicator> instance;
};

}

#endif

//! \}
