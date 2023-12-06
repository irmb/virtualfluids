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
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef MPI_NullCommunicator_H
#define MPI_NullCommunicator_H

#include "Communicator.h"

namespace vf::parallel
{

//! \brief A class implements Communicator for shared memory.
//! \details NullCommunicator is only a place-holder. It is only one process in shared memory.
class NullCommunicator : public Communicator
{
public:
    static std::shared_ptr<Communicator> getInstance();

    double Wtime() override;
    int getBundleID() const override;
    int getNumberOfBundles() const override;
    int getProcessID() const override;
    int getProcessID(int bundle, int rank) const override;
    int getNumberOfProcesses() const override;
    bool isRoot() const override;
    void *getNativeCommunicator() override;

    void sendSerializedObject(std::stringstream &stream, int target) override;
    void receiveSerializedObject(std::stringstream &stream, int source) override;

    int getRoot() const override;
    int getBundleRoot() const override;
    int getProcessRoot() const override;
    int getNumberOfProcessesInBundle(int bundle) const override;
    void barrier() override;
    void abort(int errorcode) override;

    std::vector<std::string> gather(const std::string &str) override;
    std::vector<int> gather(std::vector<int> &values) override;
    std::vector<float> gather(std::vector<float> &values) override;
    std::vector<double> gather(std::vector<double> &values) override;
    std::vector<unsigned long long> gather(std::vector<unsigned long long> &values) override;

    void allGather(std::vector<int> &svalues, std::vector<int> &rvalues) override;
    void allGather(std::vector<float> &svalues, std::vector<float> &rvalues) override;
    void allGather(std::vector<double> &svalues, std::vector<double> &rvalues) override;
    void allGather(std::vector<unsigned long long> &svalues, std::vector<unsigned long long> &rvalues) override;

    void broadcast(int &value) override;
    void broadcast(float &value) override;
    void broadcast(double &value) override;
    void broadcast(long int &value) override;
    void broadcast(std::vector<int> &values) override;
    void broadcast(std::vector<float> &values) override;
    void broadcast(std::vector<double> &values) override;
    void broadcast(std::vector<long int> &values) override;

    void receiveSend(uint *buffer_receive, int size_buffer_recv, int neighbor_rank_recv, uint *buffer_send,
                     int size_buffer_send, int neighbor_rank_send) const override;

    void send(real *sbuf, int count_s, int nb_rank) const override;

    double reduceSum(double quantityPerProcess) const override;
    int mapCudaDevicesOnHosts(const std::vector<unsigned int> &devices, int numberOfDevices) const override;

    void receiveSend(real *buffer_send, int size_buffer_send, real *buffer_receive, int size_buffer_recv,
                     int neighbor_rank) const override;

    void receiveNonBlocking(real *rbuf, int count_r, int sourceRank) override;
    void sendNonBlocking(real *sbuf, int count_s, int destinationRank) override;
    void send(real *sbuf, int count_s, int destinationRank) override;
    void waitAll() override;
    void resetRequests() override;
};

} // namespace vf::parallel

#endif
