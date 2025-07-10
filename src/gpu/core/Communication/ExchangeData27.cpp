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

#include "Calculation/Calculation.h"
#include "Cuda/CudaStreamManager.h"
#include "ExchangeData27_Device.cuh"
#include "Utilities/KernelUtilities.h"

using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition: functions used by all directions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collectNodesInSendBufferGPU(Parameter* para, int level, CudaStreamIndex streamIndex,
                                 std::vector<ProcessNeighbor27>& sendProcessNeighborsDevice)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);

    for (auto& neighbor : sendProcessNeighborsDevice) {
        GetSendFsPostDev27(para->getParD(level)->distributions.f[0],
                           neighbor.populations[0],
                           neighbor.index,
                           neighbor.numberOfNodes,
                           para->getParD(level)->neighborX,
                           para->getParD(level)->neighborY,
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes,
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads, 
                           stream);
        if(para->getDiffOn())
            GetSendFsPostDev27(para->getParD(level)->distributionsAD.f[0],
                           neighbor.populationsAD[0],
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

void scatterNodesFromRecvBufferGPU(Parameter* para, int level, CudaStreamIndex streamIndex,
                                   std::vector<ProcessNeighbor27>& recvProcessNeighborsDevice)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    for (auto& neighbor : recvProcessNeighborsDevice) {
        SetRecvFsPostDev27(para->getParD(level)->distributions.f[0], 
                           neighbor.populations[0], 
                           neighbor.index, 
                           neighbor.numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads, 
                           stream);
        if(para->getDiffOn())
            SetRecvFsPostDev27(para->getParD(level)->distributionsAD.f[0], 
                            neighbor.populationsAD[0], 
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

void startBlockingMpiSend(vf::parallel::Communicator& comm, std::vector<ProcessNeighbor27>& sendProcessNeighborHost,
                          bool diffOn)
{
    for (auto& neighbor : sendProcessNeighborHost) {
        comm.send(neighbor.populations[0], neighbor.numberOfFs, neighbor.rankNeighbor);
        if (diffOn)
            comm.send(neighbor.populationsAD[0], neighbor.numberOfFs, neighbor.rankNeighbor);
    }
}

void startNonBlockingMpiReceive(vf::parallel::Communicator& comm, std::vector<ProcessNeighbor27>& recvProcessNeighborHost,
                                bool diffOn)
{
    for (auto& neighbor : recvProcessNeighborHost) {
        comm.receiveNonBlocking(neighbor.populations[0], neighbor.numberOfFs, neighbor.rankNeighbor);
        if (diffOn)
            comm.receiveNonBlocking(neighbor.populationsAD[0], neighbor.numberOfFs, neighbor.rankNeighbor);
    }
}

void copyEdgeNodes(std::vector<LBMSimulationParameter::EdgeNodePositions>& edgeNodes,
                   std::vector<ProcessNeighbor27>& recvProcessNeighborsHost,
                   std::vector<ProcessNeighbor27>& sendProcessNeighborsHost, bool diffOn)
{
#pragma omp parallel for
    for (int i=0; i< int(edgeNodes.size()); i++) {
        const uint indexInSubdomainRecv = edgeNodes[i].indexOfProcessNeighborRecv;
        const uint indexInSubdomainSend = edgeNodes[i].indexOfProcessNeighborSend;
        const uint numNodesInBufferRecv = recvProcessNeighborsHost[indexInSubdomainRecv].numberOfNodes;
        const uint numNodesInBufferSend = sendProcessNeighborsHost[indexInSubdomainSend].numberOfNodes;
        const uint indexInSend = edgeNodes[i].indexInSendBuffer;
        const uint indexInRecv = edgeNodes[i].indexInRecvBuffer;
        if (indexInSend >= numNodesInBufferSend)
            // for reduced communication after fine to coarse: only copy send nodes which are not part of the reduced comm
            continue;

        Distributions27 populationsSend = vf::gpu::getDistributionReferences27(
            sendProcessNeighborsHost[indexInSubdomainSend].populations[0], numNodesInBufferSend, true);
        Distributions27 populationsRecv = vf::gpu::getDistributionReferences27(
            recvProcessNeighborsHost[indexInSubdomainRecv].populations[0], numNodesInBufferRecv, true);
        forEachDirection([&](auto direction) {
            (populationsSend.f[direction])[indexInSend] = (populationsRecv.f[direction])[indexInRecv];
        });
        if (diffOn) {
            Distributions27 populationsADSend = vf::gpu::getDistributionReferences27(
                sendProcessNeighborsHost[indexInSubdomainSend].populationsAD[0], numNodesInBufferSend, true);
            Distributions27 populationsADRecv = vf::gpu::getDistributionReferences27(
                recvProcessNeighborsHost[indexInSubdomainRecv].populationsAD[0], numNodesInBufferRecv, true);

            forEachDirection([&](auto direction) {
                (populationsADSend.f[direction])[indexInSend] = (populationsADRecv.f[direction])[indexInRecv];
            });
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
    startNonBlockingMpiReceive(comm, recvProcessNeighborsHost, para->getDiffOn());
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 3. before sending data, wait for memcopy (from device to host) to finish
    if (para->getUseStreams())
        cudaStreamSynchronize(stream);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 4. send data to neighboring process (MPI)
    startBlockingMpiSend(comm, sendProcessNeighborsHost, para->getDiffOn());
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
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.sendProcessNeighborsX);
}

void prepareExchangeCollDataXGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.sendProcessNeighborsAfterFtoCX);
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
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.recvProcessNeighborsX);
}

void scatterNodesFromRecvBufferXGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.recvProcessNeighborsAfterFtoCX);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataYGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.sendProcessNeighborsY);
}

void prepareExchangeCollDataYGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.sendProcessNeighborsAfterFtoCY);
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
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.recvProcessNeighborsY);
}

void scatterNodesFromRecvBufferYGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.recvProcessNeighborsAfterFtoCY);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataZGPU27AllNodes(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.sendProcessNeighborsZ);
}

void prepareExchangeCollDataZGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    collectNodesInSendBufferGPU(para, level, streamIndex, parD.sendProcessNeighborsAfterFtoCZ);
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
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.recvProcessNeighborsZ);
}

void scatterNodesFromRecvBufferZGPU27AfterFtoC(Parameter* para, int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, parD.recvProcessNeighborsAfterFtoCZ);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \}
