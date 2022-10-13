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
//! \file NullCommunicator.cpp
//! \ingroup Parallel
//! \author Konstantin Kutscher
//=======================================================================================

#include "NullCommunicator.h"

namespace vf::mpi
{

    std::shared_ptr<Communicator> NullCommunicator::getInstance()
    {
        std::lock_guard<std::mutex> myLock(instantiation_mutex);
        if (!instance){
            instance = std::shared_ptr<NullCommunicator>(new NullCommunicator);
        }
        return instance;
    }
    //////////////////////////////////////////////////////////////////////////
    int NullCommunicator::getBundleID() { return 0; }
    //////////////////////////////////////////////////////////////////////////
    int NullCommunicator::getNumberOfBundles() { return 0; }
    //////////////////////////////////////////////////////////////////////////
    int NullCommunicator::getProcessID() { return 0; }
    //////////////////////////////////////////////////////////////////////////
    int NullCommunicator::getNumberOfProcesses() { return 0; }
    //////////////////////////////////////////////////////////////////////////
    void *NullCommunicator::getNativeCommunicator() { return NULL; }
    //////////////////////////////////////////////////////////////////////////
    int NullCommunicator::getRoot() { return 0; }
    //////////////////////////////////////////////////////////////////////////
    int NullCommunicator::getBundleRoot() { return 0; }
    //////////////////////////////////////////////////////////////////////////
    int NullCommunicator::getProcessRoot() { return 0; }
    //////////////////////////////////////////////////////////////////////////
    std::vector<std::string> NullCommunicator::gather(const std::string & /*str*/) { return std::vector<std::string>(); }
    //////////////////////////////////////////////////////////////////////////

    void NullCommunicator::sendSerializedObject(std::stringstream &ss, int target) {}
    //////////////////////////////////////////////////////////////////////////
    void NullCommunicator::receiveSerializedObject(std::stringstream &ss, int source) {}

    int NullCommunicator::getProcessID(int bundle, int rank) { return 0; }
    bool NullCommunicator::isRoot() {return true; }

    int NullCommunicator::getNumberOfProcessesInBundle(int bundle) {return 0;}
    void NullCommunicator::barrier() {}
    void NullCommunicator::abort(int errorcode) {}


    std::vector<int> NullCommunicator::gather(std::vector<int> &values){ return std::vector<int>(); }
    std::vector<float> NullCommunicator::gather(std::vector<float> &values){ return std::vector<float>(); }
    std::vector<double> NullCommunicator::gather(std::vector<double> &values){ return std::vector<double>(); }
    std::vector<unsigned long long> NullCommunicator::gather(std::vector<unsigned long long> &values){ return std::vector<unsigned long long>(); }

    void NullCommunicator::allGather(std::vector<int> &svalues, std::vector<int> &rvalues){ }
    void NullCommunicator::allGather(std::vector<float> &svalues, std::vector<float> &rvalues){ }
    void NullCommunicator::allGather(std::vector<double> &svalues, std::vector<double> &rvalues){ }
    void NullCommunicator::allGather(std::vector<unsigned long long> &svalues, std::vector<unsigned long long> &rvalues){ }

    void NullCommunicator::broadcast(int &value){ }
    void NullCommunicator::broadcast(float &value){ }
    void NullCommunicator::broadcast(double &value){ }
    void NullCommunicator::broadcast(long int &value){ }
    void NullCommunicator::broadcast(std::vector<int> &values){ }
    void NullCommunicator::broadcast(std::vector<float> &values){ }
    void NullCommunicator::broadcast(std::vector<double> &values){ }
    void NullCommunicator::broadcast(std::vector<long int> &values){ }
}
