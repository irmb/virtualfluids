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
//! \file Communicator.h
//! \ingroup Parallel
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <string>
#include <vector>

#include <PointerDefinitions.h>

//! \brief An abstract class for communication between processes in parallel computation
class Communicator
{
public:
    virtual ~Communicator() = default;
    static SPtr<Communicator> getInstance();
    virtual int getBundleID()                      = 0;
    virtual int getNumberOfBundles()               = 0;
    virtual int getProcessID()                     = 0;
    virtual int getProcessID(int bundle, int rank) = 0;
    virtual int getNumberOfProcesses()             = 0;
    virtual bool isRoot()                          = 0;
    virtual void *getNativeCommunicator()          = 0;

    virtual void sendSerializedObject(std::stringstream &ss, int target)    = 0;
    virtual void receiveSerializedObject(std::stringstream &ss, int source) = 0;

    virtual int getRoot()                                = 0;
    virtual int getBundleRoot()                          = 0;
    virtual int getProcessRoot()                         = 0;
    virtual int getNumberOfProcessesInBundle(int bundle) = 0;
    virtual void barrier()                               = 0;
    virtual void abort(int errorcode)                    = 0;

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

protected:
    Communicator()                     = default;
    Communicator(const Communicator &) = default;
    static SPtr<Communicator> instance;
};

#endif
