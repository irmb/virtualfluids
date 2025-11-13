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
#include <vector>

#include "Calculation/Calculation.h"
#include "Cuda/CudaStreamManager.h"
#include "ExchangeData27_Device.cuh"
#include "Parameter/Parameter.h"
#include "Utilities/KernelUtilities.h"

using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void collectNodesInSendBufferGPU(Parameter* para, int level, CudaStreamIndex streamIndex,
                                 const std::vector<ProcessNeighbor27>& sendProcessNeighborsDevice)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);

    for (const auto& neighbor : sendProcessNeighborsDevice) {
        GetSendFsPostDev27(para->getParD(level)->distributions.f[0],
                           neighbor.populations[0],
                           para->getParD(level)->distributionsAD.f[0],
                           neighbor.populationsAD[0],
                           neighbor.index,
                           neighbor.numberOfNodes,
                           para->getParD(level)->neighborX,
                           para->getParD(level)->neighborY,
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes,
                           para->getParD(level)->isEvenTimestep, 
                           para->getDiffOn(),
                           para->getParD(level)->numberofthreads, 
                           stream);
    }
}

void scatterNodesFromRecvBufferGPU(Parameter* para, int level, CudaStreamIndex streamIndex,
                                   const std::vector<ProcessNeighbor27>& recvProcessNeighborsDevice)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    for (const auto& neighbor : recvProcessNeighborsDevice) {
        SetRecvFsPostDev27(para->getParD(level)->distributions.f[0], 
                           neighbor.populations[0], 
                           para->getParD(level)->distributionsAD.f[0], 
                           neighbor.populationsAD[0], 
                           neighbor.index, 
                           neighbor.numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getDiffOn(), 
                           para->getParD(level)->numberofthreads,
                           stream);
    }
}

void startNonBlockingMpiSend(vf::parallel::Communicator& comm, const std::vector<ProcessNeighbor27>& sendProcessNeighborsHost,
                             const bool diffOn)
{
    for (const auto& neighbor : sendProcessNeighborsHost) {
        comm.sendNonBlocking(neighbor.populations[0], neighbor.numberOfFs, neighbor.rankNeighbor);
        if (diffOn)
            comm.sendNonBlocking(neighbor.populationsAD[0], neighbor.numberOfFs, neighbor.rankNeighbor);
    }
}

void startNonBlockingMpiReceive(vf::parallel::Communicator& comm, const std::vector<ProcessNeighbor27>& recvProcessNeighborsHost,
                                const bool diffOn)
{
    for (const auto& neighbor : recvProcessNeighborsHost) {
        comm.receiveNonBlocking(neighbor.populations[0], neighbor.numberOfFs, neighbor.rankNeighbor);
        if (diffOn)
            comm.receiveNonBlocking(neighbor.populationsAD[0], neighbor.numberOfFs, neighbor.rankNeighbor);
    }
}

void copyEdgeNodes(const std::vector<LBMSimulationParameter::EdgeNodePositions>& edgeNodes,
                   const std::vector<ProcessNeighbor27>& recvProcessNeighborsHost,
                   const std::vector<ProcessNeighbor27>& sendProcessNeighborsHost, const bool diffOn)
{
#pragma omp parallel for
    for (int i = 0; i < int(edgeNodes.size()); i++) {
        const auto& edgeNode = edgeNodes[i];
        const auto& sendNeighbor = sendProcessNeighborsHost[edgeNode.indexOfProcessNeighborSend];
        const auto& recvNeighbor = recvProcessNeighborsHost[edgeNode.indexOfProcessNeighborRecv];

        if (edgeNode.indexInSendBuffer >= sendNeighbor.numberOfNodes)
            // for reduced communication after fine to coarse: only copy send nodes which are not part of the reduced comm
            continue;

        const Distributions27 populationsSend =
            vf::gpu::getDistributionReferences27(sendNeighbor.populations[0], sendNeighbor.numberOfNodes, true);
        const Distributions27 populationsRecv =
            vf::gpu::getDistributionReferences27(recvNeighbor.populations[0], recvNeighbor.numberOfNodes, true);
        forEachDirection([&](auto direction) { (populationsSend.f[direction])[edgeNode.indexInSendBuffer] = (populationsRecv.f[direction])[edgeNode.indexInRecvBuffer]; });

        if (diffOn) {
            const Distributions27 populationsADSend =
                vf::gpu::getDistributionReferences27(sendNeighbor.populationsAD[0], sendNeighbor.numberOfNodes, true);
            const Distributions27 populationsADRecv =
                vf::gpu::getDistributionReferences27(recvNeighbor.populationsAD[0], recvNeighbor.numberOfNodes, true);

            forEachDirection([&](auto direction) { (populationsADSend.f[direction])[edgeNode.indexInSendBuffer] = (populationsADRecv.f[direction])[edgeNode.indexInRecvBuffer]; });
        }
    }
}

void exchangeCollDataGPU27(Parameter* para, vf::parallel::Communicator& comm, const CudaMemoryManager* cudaMemoryManager,
                            const CudaStreamIndex streamIndex, 
                            const std::vector<ProcessNeighbor27>& sendProcessNeighborsDevice,
                            const std::vector<ProcessNeighbor27>& recvProcessNeighborsDevice,
                            const std::vector<ProcessNeighbor27>& sendProcessNeighborsHost,
                            const std::vector<ProcessNeighbor27>& recvProcessNeighborsHost,
                            const std::optional<std::vector<ProcessNeighbor27>>& recvProcessNeighborsHostX,
                            const std::optional<std::vector<LBMSimulationParameter::EdgeNodePositions>>& edgeNodesX,
                            const std::optional<std::vector<ProcessNeighbor27>>& recvProcessNeighborsHostY,
                            const std::optional<std::vector<LBMSimulationParameter::EdgeNodePositions>>& edgeNodesY)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    const size_t numberOfProcessNeighbors = sendProcessNeighborsHost.size();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! \details steps:
    //! 1. copy data from device to host
    for (size_t i = 0; i < numberOfProcessNeighbors; i++)
        cudaMemoryManager->cudaCopyProcessNeighborFsDtoH(sendProcessNeighborsHost[i], sendProcessNeighborsDevice[i]);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 2. start non-blocking receive (MPI)
    startNonBlockingMpiReceive(comm, recvProcessNeighborsHost, para->getDiffOn());
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 3. before sending data, wait for memcopy (from device to host) to finish
    if (para->getUseStreams())
        cudaStreamSynchronize(stream);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 4. copy received edge node values from x and y
    if (para->getUseStreams() && recvProcessNeighborsHostX && !recvProcessNeighborsHostX->empty() && !sendProcessNeighborsHost.empty())
        copyEdgeNodes(edgeNodesX.value(), recvProcessNeighborsHostX.value(), sendProcessNeighborsHost, para->getDiffOn());
    if (para->getUseStreams() && recvProcessNeighborsHostY && !recvProcessNeighborsHostY->empty() && !sendProcessNeighborsHost.empty())
        copyEdgeNodes(edgeNodesY.value(), recvProcessNeighborsHostY.value(), sendProcessNeighborsHost, para->getDiffOn());
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 5. send data to neighboring process (MPI) and wait
    startNonBlockingMpiSend(comm, sendProcessNeighborsHost, para->getDiffOn());
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 6. reset the request array, which was used for the mpi communication
    if (0 < numberOfProcessNeighbors)
        comm.resetRequests();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 7. copy received data from host to device
    for (size_t i = 0; i < numberOfProcessNeighbors; i++)
        cudaMemoryManager->cudaCopyProcessNeighborFsHtoD(recvProcessNeighborsHost[i], recvProcessNeighborsDevice[i]);
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

void exchangeCollDataXGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm, const CudaMemoryManager* cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);

    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsX, parD.recvProcessNeighborsX,
                          parH.sendProcessNeighborsX, parH.recvProcessNeighborsX);
}

void exchangeCollDataXGPU27AfterFtoC(Parameter* para, vf::parallel::Communicator& comm, const CudaMemoryManager* cudaMemoryManager,
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

void exchangeCollDataYGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm, const CudaMemoryManager* cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);

    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex,
                          parD.sendProcessNeighborsY, parD.recvProcessNeighborsY,
                          parH.sendProcessNeighborsY, parH.recvProcessNeighborsY,
                          parH.recvProcessNeighborsX, parH.edgeNodesXtoY);
}


void exchangeCollDataYGPU27AfterFtoC(Parameter* para, vf::parallel::Communicator& comm, const CudaMemoryManager* cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);

    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsAfterFtoCY, parD.recvProcessNeighborsAfterFtoCY,
                          parH.sendProcessNeighborsAfterFtoCY, parH.recvProcessNeighborsAfterFtoCY, 
                          parH.recvProcessNeighborsAfterFtoCX, parH.edgeNodesXtoY);
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

void exchangeCollDataZGPU27AllNodes(Parameter* para, vf::parallel::Communicator& comm, const CudaMemoryManager* cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);

    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex, 
                          parD.sendProcessNeighborsZ, parD.recvProcessNeighborsZ,
                          parH.sendProcessNeighborsZ, parH.recvProcessNeighborsZ,
                          parH.recvProcessNeighborsX, parH.edgeNodesXtoZ,
                          parH.recvProcessNeighborsY, parH.edgeNodesYtoZ);
}
void exchangeCollDataZGPU27AfterFtoC(Parameter* para, vf::parallel::Communicator& comm, const CudaMemoryManager* cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    auto& parD = para->getParDeviceAsReference(level);
    auto& parH = para->getParHostAsReference(level);

    exchangeCollDataGPU27(para, comm, cudaMemoryManager, streamIndex,
                          parD.sendProcessNeighborsAfterFtoCZ, parD.recvProcessNeighborsAfterFtoCZ, 
                          parH.sendProcessNeighborsAfterFtoCZ, parH.recvProcessNeighborsAfterFtoCZ, 
                          parH.recvProcessNeighborsAfterFtoCX, parH.edgeNodesXtoZ,
                          parH.recvProcessNeighborsAfterFtoCY, parH.edgeNodesYtoZ);
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

//! \}
