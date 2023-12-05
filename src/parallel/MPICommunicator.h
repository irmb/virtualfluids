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
#if defined VF_MPI

#ifndef MPI_MPICOMMUNICATOR_H
#define MPI_MPICOMMUNICATOR_H

#include "Communicator.h"
#include <basics/PointerDefinitions.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbLogger.h>
#include <mpi.h>
#include <string>
#include <vector>

//////////////////////////////////
#ifdef VF_DOUBLE_ACCURACY
#define VF_MPI_REAL MPI_DOUBLE
#else
#define VF_MPI_REAL MPI_FLOAT
#endif
//////////////////////////////////

namespace vf::parallel
{

//! \brief A class uses MPI library to communication.
//! \details Support MPI communication. Implements singleton pattern.
//! \author K. Kutscher
class MPICommunicator : public Communicator
{
public:
    MPICommunicator(MPICommunicator const&) = delete;
    MPICommunicator& operator=(MPICommunicator const&) = delete;

    ~MPICommunicator() override;
    static std::shared_ptr<Communicator> getInstance();
    double Wtime() override;
    int getBundleID() const override;
    int getNumberOfBundles() const override;
    int getProcessID() const override;
    int getProcessID(int bundle, int rank) const override;
    int getNumberOfProcesses() const override;
    void *getNativeCommunicator() override;
    int getRoot() const override;
    int getBundleRoot() const override;
    int getProcessRoot() const override;
    int getNumberOfProcessesInBundle(int bundle) const override;
    bool isRoot() const override;
    void abort(int errorcode) override;

    void sendSerializedObject(std::stringstream &ss, int target) override;
    void receiveSerializedObject(std::stringstream &ss, int source) override;

    void barrier() override;

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

    template <class T>
    std::vector<T> gather(std::vector<T> &values);

    template <class T>
    void allGather(std::vector<T> &svalues, std::vector<T> &rvalues);

    template <class T>
    void broadcast(std::vector<T> &values);

    template <class T>
    void broadcast(T &value);

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

private:
    MPICommunicator();

    int numprocs, PID;
    MPI_Comm comm;
    int root;

    std::vector<MPI_Request> requests;
};

//////////////////////////////////////////////////////////////////////////
template <class T>
std::vector<T> MPICommunicator::gather(std::vector<T> &values)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(unsigned long long).name())
        mpiDataType = MPI_UNSIGNED_LONG_LONG;
    else if ((std::string) typeid(T).name() == (std::string) typeid(char).name())
        mpiDataType = MPI_CHAR;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    int count = static_cast<int>(values.size());
    std::vector<T> rvalues(1);

    if (PID == root) {
        rvalues.resize(numprocs * count);
    }

    MPI_Gather(&values[0], count, mpiDataType, &rvalues[0], count, mpiDataType, root, comm);

    return rvalues;
}
//////////////////////////////////////////////////////////////////////////
template <class T>
void MPICommunicator::allGather(std::vector<T> &svalues, std::vector<T> &rvalues)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(unsigned long long).name())
        mpiDataType = MPI_UNSIGNED_LONG_LONG;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    int scount;
    std::vector<int> displs, rcounts;

    scount = (int)(svalues.size());

    rcounts.resize(numprocs);
    MPI_Allgather(&scount, 1, MPI_INT, &rcounts[0], 1, MPI_INT, comm);
    displs.resize(numprocs);

    displs[0] = 0;
    for (int i = 1; i < numprocs; ++i) {
        displs[i] = displs[i - 1] + rcounts[i - 1];
    }

    rvalues.resize(displs[numprocs - 1] + rcounts[numprocs - 1]);

    T* sval = NULL;
    T* rval = NULL;

    if (svalues.size() > 0) {
        //svalues.resize(1);
        //svalues[0] = 999;
        sval = &svalues[0];
    }

    if (rvalues.size() > 0) {
        //rvalues.resize(1);
        //rvalues[0] = 999;
        rval = &rvalues[0];
    }

    //MPI_Allgatherv(&svalues[0], scount, mpiDataType, &rvalues[0], &rcounts[0], &displs[0], mpiDataType, comm);
    MPI_Allgatherv(sval, scount, mpiDataType, rval, &rcounts[0], &displs[0], mpiDataType, comm);
}
//////////////////////////////////////////////////////////////////////////
template <class T>
void MPICommunicator::broadcast(std::vector<T> &values)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(long int).name())
        mpiDataType = MPI_LONG_INT;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    int rcount;
    if (this->PID == this->root) {
        rcount = (int)values.size();
    }

    MPI_Bcast(&rcount, 1, MPI_INT, this->root, comm);

    if (this->PID != this->root) {
        values.resize(rcount);
    }

    MPI_Bcast(&values[0], (int)values.size(), mpiDataType, this->root, comm);
}
//////////////////////////////////////////////////////////////////////////
template <class T>
void MPICommunicator::broadcast(T &value)
{
    MPI_Datatype mpiDataType;
    if ((std::string) typeid(T).name() == (std::string) typeid(double).name())
        mpiDataType = MPI_DOUBLE;
    else if ((std::string) typeid(T).name() == (std::string) typeid(float).name())
        mpiDataType = MPI_FLOAT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(int).name())
        mpiDataType = MPI_INT;
    else if ((std::string) typeid(T).name() == (std::string) typeid(long int).name())
        mpiDataType = MPI_LONG_INT;
    else
        throw UbException(UB_EXARGS, "no MpiDataType for T" + (std::string) typeid(T).name());

    MPI_Bcast(&value, 1, mpiDataType, this->root, comm);
}
//////////////////////////////////////////////////////////////////////////


#endif

}

#endif
