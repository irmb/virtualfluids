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
//! \addtogroup gpu_Communication Communication
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//======================================================================================

#include "ExchangeData27.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <parallel/Communicator.h>

#include "Cuda/CudaStreamManager.h"
#include "ExchangeData27_Device.cuh"

using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition: functions used by all directions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collectNodesInSendBufferGPU(Parameter* para, int level, CudaStreamIndex streamIndex, real* distributions,
                                 std::vector<ProcessNeighbor27>& sendProcessNeighbor)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);

    for (auto& neighbor : sendProcessNeighbor) {
        GetSendFsPostDev27(distributions,
                           neighbor.f[0],
                           neighbor.index,
                           neighbor.numberOfNodes,
                           para->getParD(level)->neighborX,
                           para->getParD(level)->neighborY,
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes,
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads, 
                           stream);
    }
}

void scatterNodesFromRecvBufferGPU(Parameter* para, int level, CudaStreamIndex streamIndex, real* distributions,
                                   std::vector<ProcessNeighbor27>& recvProcessNeighborDev)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    for (auto& neighbor : recvProcessNeighborDev) {
        SetRecvFsPostDev27(distributions, 
                           neighbor.f[0], 
                           neighbor.index, 
                           neighbor.numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads, 
                           stream);
    }
}

void startBlockingMpiSend(vf::parallel::Communicator& comm, std::vector<ProcessNeighbor27>& sendProcessNeighborHost)
{
    for (auto& neighbor : sendProcessNeighborHost)
        comm.send(neighbor.f[0], neighbor.numberOfFs, neighbor.rankNeighbor);
}

void startNonBlockingMpiReceive(vf::parallel::Communicator& comm, std::vector<ProcessNeighbor27>& recvProcessNeighborHost)
{
    for (auto& neighbor : recvProcessNeighborHost)
        comm.receiveNonBlocking(neighbor.f[0], neighbor.numberOfFs, neighbor.rankNeighbor);
}

void copyEdgeNodes(std::vector<LBMSimulationParameter::EdgeNodePositions>& edgeNodes,
                   std::vector<ProcessNeighbor27>& recvProcessNeighborHost,
                   std::vector<ProcessNeighbor27>& sendProcessNeighborHost)
{
    int indexInSubdomainRecv = 0;
    int indexInSubdomainSend = 0;
    int numNodesInBufferRecv = 0;
    int numNodesInBufferSend = 0;

#pragma omp parallel for
    for (size_t i = 0; i < edgeNodes.size(); i++) {
        indexInSubdomainRecv = edgeNodes[i].indexOfProcessNeighborRecv;
        indexInSubdomainSend = edgeNodes[i].indexOfProcessNeighborSend;
        numNodesInBufferRecv = recvProcessNeighborHost[indexInSubdomainRecv].numberOfNodes;
        numNodesInBufferSend = sendProcessNeighborHost[indexInSubdomainSend].numberOfNodes;
        if (edgeNodes[i].indexInSendBuffer >= numNodesInBufferSend) {
            // for reduced communication after fine to coarse: only copy send nodes which are not part of the reduced comm
            continue;
        }

        // copy fs for all directions
        for (size_t direction = 0; direction <= ENDDIR; direction++) {
            (sendProcessNeighborHost[indexInSubdomainSend].f[0] +
             (direction * numNodesInBufferSend))[edgeNodes[i].indexInSendBuffer] =
                (recvProcessNeighborHost[indexInSubdomainRecv].f[0] +
                 (direction * numNodesInBufferRecv))[edgeNodes[i].indexInRecvBuffer];
        }
    }
}

void exchangeCollDataGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager,
                           CudaStreamIndex streamIndex, std::vector<ProcessNeighbor27>& sendProcessNeighborsDev,
                           std::vector<ProcessNeighbor27>& recvProcessNeighborsDev,
                           std::vector<ProcessNeighbor27>& sendProcessNeighborsHost,
                           std::vector<ProcessNeighbor27>& recvProcessNeighborsHost)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    const size_t numberOfProcessNeighbors = sendProcessNeighborsHost.size();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! \details steps:
    //! 1. copy data from device to host
    for (size_t i = 0; i < numberOfProcessNeighbors; i++)
        cudaMemoryManager->cudaCopyProcessNeighborFsDtoH(&sendProcessNeighborsHost[i], &sendProcessNeighborsDev[i]);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 2. start non-blocking receive (MPI)
    startNonBlockingMpiReceive(comm, recvProcessNeighborsHost);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 3. before sending data, wait for memcopy (from device to host) to finish
    if (para->getUseStreams())
        cudaStreamSynchronize(stream);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 4. send data to neighboring process (MPI)
    startBlockingMpiSend(comm, sendProcessNeighborsHost);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 5. wait until data is received
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 6. reset the request array, which was used for the mpi communication
    if (0 < numberOfProcessNeighbors)
        comm.resetRequests();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 7. copy received data from host to device
    for (size_t i = 0; i < numberOfProcessNeighbors; i++)
        cudaMemoryManager->cudaCopyProcessNeighborFsHtoD(&recvProcessNeighborsHost[i], &recvProcessNeighborsDev[i]);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataXGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.sendProcessNeighborsX);
}

void prepareExchangeCollDataXGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.sendProcessNeighborsAfterFtoCX);
}

void exchangeCollDataXGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);
    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsX, parD.recvProcessNeighborsX,
                          parH.sendProcessNeighborsX, parH.recvProcessNeighborsX);
}

void exchangeCollDataXGPU27AfterFtoC(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);
    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex,
                          parD.sendProcessNeighborsAfterFtoCX, parD.recvProcessNeighborsAfterFtoCX, 
                          parH.sendProcessNeighborsAfterFtoCX, parH.recvProcessNeighborsAfterFtoCX);
}

void scatterNodesFromRecvBufferXGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.recvProcessNeighborsX);
}

void scatterNodesFromRecvBufferXGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.recvProcessNeighborsAfterFtoCX);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataYGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.sendProcessNeighborsY);
}

void prepareExchangeCollDataYGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.sendProcessNeighborsAfterFtoCY);
}

void exchangeCollDataYGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);
    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex,
                          parD.sendProcessNeighborsY, parD.recvProcessNeighborsY,
                          parH.sendProcessNeighborsY, parH.recvProcessNeighborsY);
}

void exchangeCollDataYGPU27AfterFtoC(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);
    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsAfterFtoCY, parD.recvProcessNeighborsAfterFtoCY,
                          parH.sendProcessNeighborsAfterFtoCY, parH.recvProcessNeighborsAfterFtoCY);
}

void scatterNodesFromRecvBufferYGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.recvProcessNeighborsY);
}

void scatterNodesFromRecvBufferYGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.recvProcessNeighborsAfterFtoCY);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataZGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.sendProcessNeighborsZ);
}

void prepareExchangeCollDataZGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.sendProcessNeighborsAfterFtoCZ);
}

void exchangeCollDataZGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);
    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsZ, parD.recvProcessNeighborsZ,
                          parH.sendProcessNeighborsZ, parH.recvProcessNeighborsZ);
}
void exchangeCollDataZGPU27AfterFtoC(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);
    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex,
                          parD.sendProcessNeighborsAfterFtoCZ, parD.recvProcessNeighborsAfterFtoCZ, 
                          parH.sendProcessNeighborsAfterFtoCZ, parH.recvProcessNeighborsAfterFtoCZ);
}

void scatterNodesFromRecvBufferZGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.recvProcessNeighborsZ);
}

void scatterNodesFromRecvBufferZGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributions.f[0], parD.recvProcessNeighborsAfterFtoCZ);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition advection diffusion
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataXADGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributionsAD.f[0], parD.sendProcessNeighborsADX);
}
void exchangeCollDataXADGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm,
                                      CudaMemoryManager* cudaMemoryManager, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);

    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsADX, parD.recvProcessNeighborsADX, 
                          parH.sendProcessNeighborsADX, parH.recvProcessNeighborsADX);
}
void scatterNodesFromRecvBufferXADGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributionsAD.f[0], parD.recvProcessNeighborsADX);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataYADGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributionsAD.f[0], parD.sendProcessNeighborsADY);
}
void exchangeCollDataYADGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm,
                                      CudaMemoryManager* cudaMemoryManager, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);

    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsADY, parD.recvProcessNeighborsADY, 
                          parH.sendProcessNeighborsADY, parH.recvProcessNeighborsADY);
}
void scatterNodesFromRecvBufferYADGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributionsAD.f[0], parD.recvProcessNeighborsADY);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataZADGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.distributionsAD.f[0], parD.sendProcessNeighborsADZ);
}
void exchangeCollDataZADGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm,
                                      CudaMemoryManager* cudaMemoryManager, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);
    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsADZ, parD.recvProcessNeighborsADZ, 
                          parH.sendProcessNeighborsADZ, parH.recvProcessNeighborsADZ);
}
void scatterNodesFromRecvBufferZADGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.distributionsAD.f[0], parD.recvProcessNeighborsADZ);
}

//! \}
