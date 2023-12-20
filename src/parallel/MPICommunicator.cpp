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
#if defined VF_MPI
#if defined (_WIN32) || defined (_WIN64)
   #include <Winsock2.h>
#elif defined (__unix__)
   #include <unistd.h>
#endif

#include "MPICommunicator.h"

#include <mpi.h>

#include <sstream>

#include <logger/Logger.h>

using namespace std;

namespace vf::parallel
{
//////////////////////////////////////////////////////////////////////////
double MPICommunicator::Wtime()
{
    return MPI_Wtime();
}
std::shared_ptr<Communicator> MPICommunicator::getInstance()
{
    std::lock_guard<std::mutex> myLock(instantiation_mutex);
    if (!instance) {
        instance = std::shared_ptr<MPICommunicator>(new MPICommunicator);
    }
    return instance;
}
//////////////////////////////////////////////////////////////////////////
MPICommunicator::MPICommunicator()
{
    // proof if MPI is initialized
    int mpiInitialized = 0; // false
    MPI_Initialized(&mpiInitialized);
    if (mpiInitialized == 0) {
        MPI_Init(NULL, NULL);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &PID);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    comm = MPI_COMM_WORLD;
    root = 0;
}
//////////////////////////////////////////////////////////////////////////
MPICommunicator::~MPICommunicator()
{
    // proof if MPI is finalized
    int _mpiFinalized = 0; // false
    MPI_Finalized(&_mpiFinalized);
    if (_mpiFinalized == 0) {
        MPI_Finalize();
    }
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::abort(int errorcode) { MPI_Abort(comm, errorcode); }
////////////////////////////////////////////////////////////////////////////
vector<string> MPICommunicator::gather(const string &str)
{
    vector<string> parts;
    vector<string> strings;
    int scount;
    vector<char> rbuf(1);
    vector<int> rcounts(1);
    MPI_Status status;

    if (PID == root) {
        rcounts.resize(numprocs - 1);
        strings.push_back(str);

        for (int i = 1; i < numprocs; i++) {
            MPI_Recv(&rcounts[i - 1], 1, MPI_INT, i, 0, comm, &status);
        }
        for (int i = 1; i < numprocs; i++) {
            rbuf.resize(rcounts[i - 1]);
            MPI_Recv(&rbuf[0], rcounts[i - 1], MPI_CHAR, i, 0, comm, &status);
            string s(&rbuf[0], rcounts[i - 1]);
            if (s != "")
                strings.push_back(s);
        }
    } else {
        scount = (int)str.length();
        MPI_Send(&scount, 1, MPI_INT, root, 0, comm);
        MPI_Send((char *)str.c_str(), scount, MPI_CHAR, root, 0, comm);
    }
    return strings;
}
//////////////////////////////////////////////////////////////////////////
vector<int> MPICommunicator::gather(vector<int> &values) { return gather<int>(values); }
//////////////////////////////////////////////////////////////////////////
vector<float> MPICommunicator::gather(vector<float> &values) { return gather<float>(values); }
//////////////////////////////////////////////////////////////////////////
vector<double> MPICommunicator::gather(vector<double> &values) { return gather<double>(values); }
//////////////////////////////////////////////////////////////////////////
std::vector<unsigned long long> MPICommunicator::gather(std::vector<unsigned long long> &values)
{
    return gather<unsigned long long>(values);
}
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getProcessID() const { return PID; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getProcessID(int /*bundle*/, int /*rank*/) const { return PID; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getNumberOfProcesses() const { return numprocs; }
//////////////////////////////////////////////////////////////////////////
void *MPICommunicator::getNativeCommunicator() { return &comm; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getBundleID() const { return 0; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getNumberOfBundles() const { return 1; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getRoot() const { return root; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getBundleRoot() const { return 0; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getProcessRoot() const { return 0; }
//////////////////////////////////////////////////////////////////////////
int MPICommunicator::getNumberOfProcessesInBundle(int /*bundle*/) const { return numprocs; }
//////////////////////////////////////////////////////////////////////////
bool MPICommunicator::isRoot() const { return PID == root; }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::sendSerializedObject(std::stringstream &ss, int target)
{
    string str = ss.str();
    int scount = static_cast<int>(str.length());
    MPI_Send(&scount, 1, MPI_INT, target, 0, comm);
    MPI_Send((char *)str.c_str(), scount, MPI_CHAR, target, 0, comm);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::receiveSerializedObject(std::stringstream &ss, int source)
{
    vector<char> rbuf;
    int rcount;
    MPI_Status status;
    MPI_Recv(&rcount, 1, MPI_INT, source, 0, comm, &status);
    rbuf.resize(rcount);
    MPI_Recv(&rbuf[0], rcount, MPI_CHAR, source, 0, comm, &status);
    ss.rdbuf()->pubsetbuf(&rbuf[0], rcount);
    string str(&rbuf[0]);
    ss.str(str);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::barrier() { MPI_Barrier(comm); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<int> &svalues, std::vector<int> &rvalues)
{
    allGather<int>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<float> &svalues, std::vector<float> &rvalues)
{
    allGather<float>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<double> &svalues, std::vector<double> &rvalues)
{
    allGather<double>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::allGather(std::vector<unsigned long long> &svalues, std::vector<unsigned long long> &rvalues)
{
    allGather<unsigned long long>(svalues, rvalues);
}
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<int> &values) { broadcast<int>(values); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<float> &values) { broadcast<float>(values); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<double> &values) { broadcast<double>(values); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(std::vector<long int> &values) { broadcast<long int>(values); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(int &value) { broadcast<int>(value); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(float &value) { broadcast<float>(value); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(double &value) { broadcast<double>(value); }
//////////////////////////////////////////////////////////////////////////
void MPICommunicator::broadcast(long int &value) { broadcast<long int>(value); }

void MPICommunicator::receiveSend(uint *buffer_receive, int size_buffer_recv,
                                  int neighbor_rank_recv, uint *buffer_send, int size_buffer_send,
                                  int neighbor_rank_send) const
{
    MPI_Request recv_request;
    MPI_Irecv(buffer_receive, size_buffer_recv, MPI_UNSIGNED, neighbor_rank_recv, 0, comm,
              &recv_request);
    // printf("receive_send PID: %i,   nbRev: nb_rank_recv: %i", this->getPID(), nb_rank_r);
    // fflush(stdout);
    MPI_Send(buffer_send, size_buffer_send, MPI_UNSIGNED, neighbor_rank_send, 0, comm);
    // printf("receive_send PID: %i,   sendUintGPU: nb_rank_send: %i", this->getPID(), nb_rank_s);
    // fflush(stdout);
    MPI_Wait(&recv_request, MPI_STATUSES_IGNORE); // TODO: Do we have a benefit here or could we simply do a blocking receiv.
}

void MPICommunicator::receiveSend(real *buffer_send, int size_buffer_send, real *buffer_receive, int size_buffer_recv,
                     int neighbor_rank) const
{
    MPI_Send(buffer_send, size_buffer_send, VF_MPI_REAL, neighbor_rank, 0, comm);
    MPI_Recv(buffer_receive, size_buffer_recv, VF_MPI_REAL, neighbor_rank, 0, comm, MPI_STATUS_IGNORE);
}

void MPICommunicator::send(real *sbuf, int count_s, int nb_rank) const
{
    MPI_Send(sbuf, count_s, VF_MPI_REAL, nb_rank, 0, comm);
}

double MPICommunicator::reduceSum(double quantityPerProcess) const
{
    double *buffer_send = &quantityPerProcess;
    double *buffer_recv = (double *)malloc(sizeof(double));

    MPI_Reduce(buffer_send, buffer_recv, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    return *buffer_recv;
}

int MPICommunicator::mapCudaDevicesOnHosts(const std::vector<unsigned int> &devices, int numberOfDevices) const
{
    int device        = -1;
    char *host        = (char *)malloc(sizeof(char) * getNumberOfProcesses() * 255);
    unsigned int *map = (unsigned int *)malloc(sizeof(unsigned int) * getNumberOfProcesses());

    char hostname[255];
    gethostname(hostname, 254);
    hostname[254] = 0;

    MPI_Gather(hostname, 255, MPI_BYTE, host, 255, MPI_BYTE, 0, MPI_COMM_WORLD);

    int i, j;
    if (isRoot()) {
        for (i = 0; i < getNumberOfProcesses(); i++) {
            int counter = 0;
            for (j = 0; j < i; j++) {
                if (strcmp(&host[i * 255], &host[j * 255]) == 0)
                    counter++;
            }
            if (counter >= numberOfDevices) {
                VF_LOG_CRITICAL("More processes than GPUs!");
                exit(1);
            }
            map[i] = devices[counter];
        }
    }

    MPI_Scatter(map, 1, MPI_UNSIGNED, &device, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    VF_LOG_INFO("Rank: {} runs on host: {} with GPU: {}", getProcessID(), hostname, device);

    free(map);
    free(host);
    return device;
}

void MPICommunicator::receiveNonBlocking(real *rbuf, int count_r, int sourceRank)
{
    // printf("\n Start Recv Rank: %d, neighbor Rank: %d, request = %d \n", PID, nb_rank, (int)requestGPU.size());
    // fflush(stdout);

    MPI_Request request;
    MPI_Irecv(rbuf, count_r, VF_MPI_REAL, sourceRank, 0, comm, &request);
    requests.push_back(request);

    // printf("\n End Recv - Rank: %d , neighbor Rank: %d \n", PID, nb_rank);
    // fflush(stdout);
}

void MPICommunicator::sendNonBlocking(real *sbuf, int count_s, int destinationRank)
{
    // printf("\n Start Send Rank: %d, neighbor Rank: %d, request = %d \n", PID, nb_rank, (int)requestGPU.size());
    // fflush(stdout);

    MPI_Request request;
    MPI_Isend(sbuf, count_s, VF_MPI_REAL, destinationRank, 0, comm, &request);
    requests.push_back(request);
    // printf("\n End Send - Rank: %d , neighbor Rank: %d \n", PID, nb_rank);
    // fflush(stdout);
}

void MPICommunicator::send(real *sbuf, int count_s, int destinationRank)
{
    MPI_Send(sbuf, count_s, VF_MPI_REAL, destinationRank, 0, comm);
}

void MPICommunicator::waitAll()
{
    MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

void MPICommunicator::resetRequests()
{
    requests.clear();
}

}

#endif

//! \}
