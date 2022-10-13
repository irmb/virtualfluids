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
//! \file NullCommunicator.h
//! \ingroup Parallel
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef MPI_NullCommunicator_H
#define MPI_NullCommunicator_H

#include "Communicator.h"

namespace vf::mpi
{

//! \brief A class implements Communicator for shared memory.
//! \details NullCommunicator is only a place-holder. It is only one process in shared memory.
class NullCommunicator : public Communicator
{
public:
    static std::shared_ptr<Communicator> getInstance();

    int getBundleID();
    int getNumberOfBundles();
    int getProcessID();
    int getProcessID(int bundle, int rank);
    int getNumberOfProcesses();
    bool isRoot();
    void *getNativeCommunicator();

    void sendSerializedObject(std::stringstream &ss, int target);
    void receiveSerializedObject(std::stringstream &ss, int source);

    int getRoot();
    int getBundleRoot();
    int getProcessRoot();
    int getNumberOfProcessesInBundle(int bundle);
    void barrier();
    void abort(int errorcode);

    std::vector<std::string> gather(const std::string &str);
    std::vector<int> gather(std::vector<int> &values);
    std::vector<float> gather(std::vector<float> &values);
    std::vector<double> gather(std::vector<double> &values);
    std::vector<unsigned long long> gather(std::vector<unsigned long long> &values);

    void allGather(std::vector<int> &svalues, std::vector<int> &rvalues);
    void allGather(std::vector<float> &svalues, std::vector<float> &rvalues);
    void allGather(std::vector<double> &svalues, std::vector<double> &rvalues);
    void allGather(std::vector<unsigned long long> &svalues, std::vector<unsigned long long> &rvalues);

    void broadcast(int &value);
    void broadcast(float &value);
    void broadcast(double &value);
    void broadcast(long int &value);
    void broadcast(std::vector<int> &values);
    void broadcast(std::vector<float> &values);
    void broadcast(std::vector<double> &values);
    void broadcast(std::vector<long int> &values);
};

}

#endif
