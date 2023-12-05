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
//! \author Martin Schoenherr
//======================================================================================

#include "ExchangeData27.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <parallel/Communicator.h>

#include "ExchangeData27_Device.cuh"
#include "Parameter/CudaStreamManager.h"

using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D domain decomposition: functions used by all directions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collectNodesInSendBufferGPU(Parameter *para, int level, CudaStreamIndex streamIndex,
                                 std::vector<ProcessNeighbor27> *sendProcessNeighbor,
                                 unsigned int numberOfSendProcessNeighbors)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    
    for (unsigned int i = 0; i < numberOfSendProcessNeighbors; i++) {
        GetSendFsPostDev27(para->getParD(level)->distributions.f[0], 
                           (*sendProcessNeighbor)[i].f[0],
                           (*sendProcessNeighbor)[i].index, 
                           (*sendProcessNeighbor)[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY,
                           para->getParD(level)->neighborZ, 
                           para->getParD(level)->numberOfNodes,
                           para->getParD(level)->isEvenTimestep, 
                           para->getParD(level)->numberofthreads, 
                           stream);
    }
}

void scatterNodesFromRecvBufferGPU(Parameter *para, int level, CudaStreamIndex streamIndex,
                                   std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                   unsigned int numberOfRecvProcessNeighbors)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    for (unsigned int i = 0; i < numberOfRecvProcessNeighbors; i++) {
        SetRecvFsPostDev27(para->getParD(level)->distributions.f[0], 
                           (*recvProcessNeighborDev)[i].f[0],
                           (*recvProcessNeighborDev)[i].index, 
                           (*recvProcessNeighborDev)[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY,
                           para->getParD(level)->neighborZ, 
                           para->getParD(level)->numberOfNodes,
                           para->getParD(level)->isEvenTimestep, 
                           para->getParD(level)->numberofthreads, 
                           stream);
    }
}

void startBlockingMpiSend(unsigned int numberOfSendProcessNeighbors, vf::parallel::Communicator &comm,
                          std::vector<ProcessNeighbor27> *sendProcessNeighborHost)
{
    for (unsigned int i = 0; i < numberOfSendProcessNeighbors; i++) {
        comm.send((*sendProcessNeighborHost)[i].f[0], (*sendProcessNeighborHost)[i].numberOfFs,
                  (*sendProcessNeighborHost)[i].rankNeighbor);
    }
}

void startNonBlockingMpiReceive(unsigned int numberOfSendProcessNeighbors, vf::parallel::Communicator &comm,
                                std::vector<ProcessNeighbor27> *recvProcessNeighborHost)
{
    for (unsigned int i = 0; i < numberOfSendProcessNeighbors; i++) {
        comm.receiveNonBlocking((*recvProcessNeighborHost)[i].f[0], (*recvProcessNeighborHost)[i].numberOfFs,
                                (*recvProcessNeighborHost)[i].rankNeighbor);
    }
}

void copyEdgeNodes(std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes,
                   std::vector<ProcessNeighbor27> &recvProcessNeighborHost,
                   std::vector<ProcessNeighbor27> &sendProcessNeighborHost)
{
    int indexInSubdomainRecv = 0;
    int indexInSubdomainSend = 0;
    int numNodesInBufferRecv = 0;
    int numNodesInBufferSend = 0;

#pragma omp parallel for
    for (int i = 0; i < (int)edgeNodes.size(); i++) {
        indexInSubdomainRecv = edgeNodes[i].indexOfProcessNeighborRecv;
        indexInSubdomainSend = edgeNodes[i].indexOfProcessNeighborSend;
        numNodesInBufferRecv = recvProcessNeighborHost[indexInSubdomainRecv].numberOfNodes;
        numNodesInBufferSend = sendProcessNeighborHost[indexInSubdomainSend].numberOfNodes;
        if(edgeNodes[i].indexInSendBuffer >= numNodesInBufferSend){
            // for reduced communication after fine to coarse: only copy send nodes which are not part of the reduced comm
            continue;
        }

        // copy fs for all directions
        for (int direction = 0; direction <= (int)ENDDIR; direction++) {
            (sendProcessNeighborHost[indexInSubdomainSend].f[0] +
             (direction * numNodesInBufferSend))[edgeNodes[i].indexInSendBuffer] =
                (recvProcessNeighborHost[indexInSubdomainRecv].f[0] +
                 (direction * numNodesInBufferRecv))[edgeNodes[i].indexInRecvBuffer];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataXGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    collectNodesInSendBufferGPU(para, level, streamIndex, &para->getParD(level)->sendProcessNeighborX,
                                (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")));
}

void prepareExchangeCollDataXGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    collectNodesInSendBufferGPU(para, level, streamIndex, &para->getParD(level)->sendProcessNeighborsAfterFtoCX,
                                (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")));
}

void exchangeCollDataXGPU27AllNodes(Parameter *para, vf::parallel::Communicator &comm, CudaMemoryManager *cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    exchangeCollDataXGPU27(para, comm, cudaMemoryManager, level, streamIndex,
                           &para->getParD(level)->sendProcessNeighborX,
                           &para->getParD(level)->recvProcessNeighborX,
                           &para->getParH(level)->sendProcessNeighborX,
                           &para->getParH(level)->recvProcessNeighborX);
}

void exchangeCollDataXGPU27AfterFtoC(Parameter *para, vf::parallel::Communicator &comm, CudaMemoryManager *cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    exchangeCollDataXGPU27(para, comm, cudaMemoryManager, level, streamIndex,
                           &para->getParD(level)->sendProcessNeighborsAfterFtoCX,
                           &para->getParD(level)->recvProcessNeighborsAfterFtoCX,
                           &para->getParH(level)->sendProcessNeighborsAfterFtoCX,
                           &para->getParH(level)->recvProcessNeighborsAfterFtoCX);
}

void scatterNodesFromRecvBufferXGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    scatterNodesFromRecvBufferGPU(para, level, streamIndex,&para->getParD(level)->recvProcessNeighborX,
                                  (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")));
}

void scatterNodesFromRecvBufferXGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    scatterNodesFromRecvBufferGPU(para, level, streamIndex,&para->getParD(level)->recvProcessNeighborsAfterFtoCX,
                                  (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")));
}

void exchangeCollDataXGPU27(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager, 
                            int level, CudaStreamIndex streamIndex,
                            std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                            std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                            std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                            std::vector<ProcessNeighbor27> *recvProcessNeighborHost)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! \details steps:
    //! 1. copy data from device to host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
        cudaMemoryManager->cudaCopyProcessNeighborXFsDH(level, i, (*sendProcessNeighborDev)[i].memsizeFs);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 2. start non-blocking receive (MPI)
    startNonBlockingMpiReceive((unsigned int)(*sendProcessNeighborHost).size(), comm, recvProcessNeighborHost);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 3. before sending data, wait for memcopy (from device to host) to finish
    if (para->getUseStreams()) cudaStreamSynchronize(stream); 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 4. send data to neighboring process (MPI)
    startBlockingMpiSend((unsigned int)(*sendProcessNeighborHost).size(), comm, sendProcessNeighborHost);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 5. wait for until data is received
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 6. reset the request array, which was used for the mpi communication
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send"))) comm.resetRequests();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! 7. copy received data from host to device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
        cudaMemoryManager->cudaCopyProcessNeighborXFsHD(level, i, (*recvProcessNeighborDev)[i].memsizeFs);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataYGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    collectNodesInSendBufferGPU(para, level, streamIndex, &para->getParD(level)->sendProcessNeighborY,
                                (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")));
}

void prepareExchangeCollDataYGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    collectNodesInSendBufferGPU(para, level, streamIndex, &para->getParD(level)->sendProcessNeighborsAfterFtoCY,
                                (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")));
}

void exchangeCollDataYGPU27AllNodes(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    exchangeCollDataYGPU27(para, comm, cudaMemoryManager, level, streamIndex,
                           &para->getParD(level)->sendProcessNeighborY,
                           &para->getParD(level)->recvProcessNeighborY, 
                           &para->getParH(level)->sendProcessNeighborY,
                           &para->getParH(level)->recvProcessNeighborY);
}

void exchangeCollDataYGPU27AfterFtoC(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    exchangeCollDataYGPU27(para, comm, cudaMemoryManager, level, streamIndex,
                           &para->getParD(level)->sendProcessNeighborsAfterFtoCY,
                           &para->getParD(level)->recvProcessNeighborsAfterFtoCY, 
                           &para->getParH(level)->sendProcessNeighborsAfterFtoCY,
                           &para->getParH(level)->recvProcessNeighborsAfterFtoCY);
}

void scatterNodesFromRecvBufferYGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, &para->getParD(level)->recvProcessNeighborY,
                                  (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")));
}

void scatterNodesFromRecvBufferYGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, &para->getParD(level)->recvProcessNeighborsAfterFtoCY,
                                  (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")));
}

void exchangeCollDataYGPU27(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager, int level,
                            CudaStreamIndex streamIndex,
                            std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                            std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                            std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                            std::vector<ProcessNeighbor27> *recvProcessNeighborHost)
{
    cudaStream_t stream = para->getStreamManager()->getStream(CudaStreamIndex::SubDomainBorder);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
        cudaMemoryManager->cudaCopyProcessNeighborYFsDH(level, i, (*sendProcessNeighborDev)[i].memsizeFs);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    startNonBlockingMpiReceive((unsigned int)(*sendProcessNeighborHost).size(), comm, recvProcessNeighborHost);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // wait for memcopy device to host to finish before sending data
    if (para->getUseStreams())
        cudaStreamSynchronize(stream);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // edge nodes: copy received node values from x
    if (para->getUseStreams() && para->getNumberOfProcessNeighborsX(level, "recv") > 0 &&
        para->getParH(level)->sendProcessNeighborY.size() != 0) {
        if (para->getParH(level)->sendProcessNeighborY[0].numberOfNodes ==
            (*sendProcessNeighborHost)[0].numberOfNodes) {
            // check if in communication of all nodes (as opposed to reduced communication after fine to coarse)
            copyEdgeNodes(para->getParH(level)->edgeNodesXtoY, para->getParH(level)->recvProcessNeighborX,
                          *sendProcessNeighborHost);
        } else {
            copyEdgeNodes(para->getParH(level)->edgeNodesXtoY, para->getParH(level)->recvProcessNeighborsAfterFtoCX,
                          *sendProcessNeighborHost);
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    startBlockingMpiSend((unsigned int)(*sendProcessNeighborHost).size(), comm, sendProcessNeighborHost);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send"))) comm.resetRequests();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++) {
        cudaMemoryManager->cudaCopyProcessNeighborYFsHD(level, i, (*recvProcessNeighborDev)[i].memsizeFs);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataZGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    collectNodesInSendBufferGPU(para, level, streamIndex, &para->getParD(level)->sendProcessNeighborZ,
                                (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")));
}

void prepareExchangeCollDataZGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    collectNodesInSendBufferGPU(para, level, streamIndex, &para->getParD(level)->sendProcessNeighborsAfterFtoCZ,
                                (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")));
}

void exchangeCollDataZGPU27AllNodes(Parameter *para, vf::parallel::Communicator &comm, CudaMemoryManager *cudaMemoryManager,
                                    int level, CudaStreamIndex streamIndex)
{
    exchangeCollDataZGPU27(para, comm, cudaMemoryManager, level, streamIndex,
                           &para->getParD(level)->sendProcessNeighborZ,
                           &para->getParD(level)->recvProcessNeighborZ, 
                           &para->getParH(level)->sendProcessNeighborZ,
                           &para->getParH(level)->recvProcessNeighborZ);
}
void exchangeCollDataZGPU27AfterFtoC(Parameter *para, vf::parallel::Communicator &comm, CudaMemoryManager *cudaMemoryManager,
                                     int level, CudaStreamIndex streamIndex)
{
    exchangeCollDataZGPU27(para, comm, cudaMemoryManager, level, streamIndex,
                           &para->getParD(level)->sendProcessNeighborsAfterFtoCZ,
                           &para->getParD(level)->recvProcessNeighborsAfterFtoCZ, 
                           &para->getParH(level)->sendProcessNeighborsAfterFtoCZ,
                           &para->getParH(level)->recvProcessNeighborsAfterFtoCZ);
}

void scatterNodesFromRecvBufferZGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, &para->getParD(level)->recvProcessNeighborZ,
                                  (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")));
}

void scatterNodesFromRecvBufferZGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex)
{
    scatterNodesFromRecvBufferGPU(para, level, streamIndex, &para->getParD(level)->recvProcessNeighborsAfterFtoCZ,
                                  (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataZGPU27(Parameter *para, vf::parallel::Communicator &comm, CudaMemoryManager *cudaMemoryManager, int level, 
                            CudaStreamIndex streamIndex,
                            std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                            std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                            std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                            std::vector<ProcessNeighbor27> *recvProcessNeighborHost)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
        cudaMemoryManager->cudaCopyProcessNeighborZFsDH(level, i, (*sendProcessNeighborDev)[i].memsizeFs);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    startNonBlockingMpiReceive((unsigned int)(*sendProcessNeighborHost).size(), comm, recvProcessNeighborHost);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // wait for memcopy device to host to finish before sending data
    if (para->getUseStreams())
        cudaStreamSynchronize(stream);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // edge nodes: copy received node values from x
    if (para->getUseStreams() && para->getNumberOfProcessNeighborsX(level, "recv") > 0 &&
        para->getParH(level)->sendProcessNeighborZ.size() != 0) {
        if (para->getParH(level)->sendProcessNeighborZ[0].numberOfNodes ==
            (*sendProcessNeighborHost)[0].numberOfNodes) {
            // check if in communication of all nodes (as opposed to reduced communication after fine to coarse)
            copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, para->getParH(level)->recvProcessNeighborX,
                          *sendProcessNeighborHost);
        } else {
            copyEdgeNodes(para->getParH(level)->edgeNodesXtoZ, para->getParH(level)->recvProcessNeighborsAfterFtoCX,
                          *sendProcessNeighborHost);
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // edge nodes: copy received node values from y
    if (para->getUseStreams() && para->getNumberOfProcessNeighborsY(level, "recv") > 0 &&
        para->getParH(level)->sendProcessNeighborZ.size() != 0) {
        if (para->getParH(level)->sendProcessNeighborZ[0].numberOfNodes ==
            (*sendProcessNeighborHost)[0].numberOfNodes) {
            // check if in communication of all nodes (as opposed to reduced communication after fine to coarse)
            copyEdgeNodes(para->getParH(level)->edgeNodesYtoZ, para->getParH(level)->recvProcessNeighborY,
                          *sendProcessNeighborHost);
        } else {
            copyEdgeNodes(para->getParH(level)->edgeNodesYtoZ, para->getParH(level)->recvProcessNeighborsAfterFtoCY,
                          *sendProcessNeighborHost);
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    startBlockingMpiSend((unsigned int)(*sendProcessNeighborHost).size(), comm, sendProcessNeighborHost);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send"))) comm.resetRequests();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborZFsHD(level, i, (*recvProcessNeighborDev)[i].memsizeFs);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//1D domain decomposition
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighbors(level, "send")); i++)
    {
        //////////////////////////////////////////////////////////////////////////
        GetSendFsPreDev27(para->getParD(level)->distributions.f[0],
                          para->getParD(level)->sendProcessNeighbor[i].f[0],
                          para->getParD(level)->sendProcessNeighbor[i].index,
                          para->getParD(level)->sendProcessNeighbor[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborFsDH(level, i);
        //////////////////////////////////////////////////////////////////////////
        comm.receiveSend(para->getParH(level)->sendProcessNeighbor[i].f[0], 
                            para->getParH(level)->sendProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighbor[i].f[0],
                            para->getParH(level)->recvProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->sendProcessNeighbor[i].rankNeighbor);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->distributions.f[0],
                          para->getParD(level)->recvProcessNeighbor[i].f[0],
                          para->getParD(level)->recvProcessNeighbor[i].index,
                          para->getParD(level)->recvProcessNeighbor[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighbors(level, "send")); i++)
    {
        //////////////////////////////////////////////////////////////////////////
        GetSendFsPostDev27(para->getParD(level)->distributions.f[0],
                           para->getParD(level)->sendProcessNeighbor[i].f[0],
                           para->getParD(level)->sendProcessNeighbor[i].index,
                           para->getParD(level)->sendProcessNeighbor[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborFsDH(level, i);
        //////////////////////////////////////////////////////////////////////////
        comm.receiveSend(para->getParH(level)->sendProcessNeighbor[i].f[0], 
                            para->getParH(level)->sendProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighbor[i].f[0],
                            para->getParH(level)->recvProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->sendProcessNeighbor[i].rankNeighbor);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->distributions.f[0],
                           para->getParD(level)->recvProcessNeighbor[i].f[0],
                           para->getParD(level)->recvProcessNeighbor[i].index,
                           para->getParD(level)->recvProcessNeighbor[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePostCollDataXGPU27(Parameter* para, vf::parallel::Communicator& comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPostDev27( para->getParD(level)->distributions.f[0],
//                            para->getParD(level)->sendProcessNeighborX[i].f[0],
//                            para->getParD(level)->sendProcessNeighborX[i].index,
//                            para->getParD(level)->sendProcessNeighborX[i].numberOfNodes,
//                            para->getParD(level)->neighborX, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ,
//                            para->getParD(level)->numberOfNodes, 
//                            para->getParD(level)->isEvenTimestep,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborXFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm.receiveSend(para->getParH(level)->sendProcessNeighborX[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborX[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborX[i].f[0],
//                            para->getParH(level)->recvProcessNeighborX[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborX[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborXFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPostDev27( para->getParD(level)->distributions.f[0],
//                            para->getParD(level)->recvProcessNeighborX[i].f[0],
//                            para->getParD(level)->recvProcessNeighborX[i].index,
//                            para->getParD(level)->recvProcessNeighborX[i].numberOfNodes,
//                            para->getParD(level)->neighborX, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ,
//                            para->getParD(level)->numberOfNodes, 
//                            para->getParD(level)->isEvenTimestep,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Y
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePreCollDataYGPU27(Parameter* para, vf::parallel::Communicator& comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPreDev27(para->getParD(level)->distributions.f[0],
//                          para->getParD(level)->sendProcessNeighborY[i].f[0],
//                          para->getParD(level)->sendProcessNeighborY[i].index,
//                          para->getParD(level)->sendProcessNeighborY[i].numberOfNodes,
//                          para->getParD(level)->neighborX, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ,
//                          para->getParD(level)->numberOfNodes, 
//                          para->getParD(level)->isEvenTimestep,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm.receiveSend(para->getParH(level)->sendProcessNeighborY[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborY[i].f[0],
//                            para->getParH(level)->recvProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborY[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPreDev27(para->getParD(level)->distributions.f[0],
//                          para->getParD(level)->recvProcessNeighborY[i].f[0],
//                          para->getParD(level)->recvProcessNeighborY[i].index,
//                          para->getParD(level)->recvProcessNeighborY[i].numberOfNodes,
//                          para->getParD(level)->neighborX, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ,
//                          para->getParD(level)->numberOfNodes, 
//                          para->getParD(level)->isEvenTimestep,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePostCollDataYGPU27(Parameter* para, vf::parallel::Communicator& comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPostDev27( para->getParD(level)->distributions.f[0],
//                            para->getParD(level)->sendProcessNeighborY[i].f[0],
//                            para->getParD(level)->sendProcessNeighborY[i].index,
//                            para->getParD(level)->sendProcessNeighborY[i].numberOfNodes,
//                            para->getParD(level)->neighborX, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ,
//                            para->getParD(level)->numberOfNodes, 
//                            para->getParD(level)->isEvenTimestep,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm.receiveSend(para->getParH(level)->sendProcessNeighborY[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborY[i].f[0],
//                            para->getParH(level)->recvProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborY[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPostDev27( para->getParD(level)->distributions.f[0],
//                            para->getParD(level)->recvProcessNeighborY[i].f[0],
//                            para->getParD(level)->recvProcessNeighborY[i].index,
//                            para->getParD(level)->recvProcessNeighborY[i].numberOfNodes,
//                            para->getParD(level)->neighborX, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ,
//                            para->getParD(level)->numberOfNodes, 
//                            para->getParD(level)->isEvenTimestep,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Z
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePreCollDataZGPU27(Parameter* para, vf::parallel::Communicator& comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPreDev27(para->getParD(level)->distributions.f[0],
//                          para->getParD(level)->sendProcessNeighborZ[i].f[0],
//                          para->getParD(level)->sendProcessNeighborZ[i].index,
//                          para->getParD(level)->sendProcessNeighborZ[i].numberOfNodes,
//                          para->getParD(level)->neighborX, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ,
//                          para->getParD(level)->numberOfNodes, 
//                          para->getParD(level)->isEvenTimestep,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm.receiveSend(para->getParH(level)->sendProcessNeighborZ[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborZ[i].f[0],
//                            para->getParH(level)->recvProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborZ[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPreDev27(para->getParD(level)->distributions.f[0],
//                          para->getParD(level)->recvProcessNeighborZ[i].f[0],
//                          para->getParD(level)->recvProcessNeighborZ[i].index,
//                          para->getParD(level)->recvProcessNeighborZ[i].numberOfNodes,
//                          para->getParD(level)->neighborX, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ,
//                          para->getParD(level)->numberOfNodes, 
//                          para->getParD(level)->isEvenTimestep,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePostCollDataZGPU27(Parameter* para, vf::parallel::Communicator& comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPostDev27( para->getParD(level)->distributions.f[0],
//                            para->getParD(level)->sendProcessNeighborZ[i].f[0],
//                            para->getParD(level)->sendProcessNeighborZ[i].index,
//                            para->getParD(level)->sendProcessNeighborZ[i].numberOfNodes,
//                            para->getParD(level)->neighborX, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ,
//                            para->getParD(level)->numberOfNodes, 
//                            para->getParD(level)->isEvenTimestep,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm.receiveSend(para->getParH(level)->sendProcessNeighborZ[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborZ[i].f[0],
//                            para->getParH(level)->recvProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborZ[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPostDev27( para->getParD(level)->distributions.f[0],
//                            para->getParD(level)->recvProcessNeighborZ[i].f[0],
//                            para->getParD(level)->recvProcessNeighborZ[i].index,
//                            para->getParD(level)->recvProcessNeighborZ[i].numberOfNodes,
//                            para->getParD(level)->neighborX, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ,
//                            para->getParD(level)->numberOfNodes, 
//                            para->getParD(level)->isEvenTimestep,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


















































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3D domain decomposition convection diffusion
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataADXGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        GetSendFsPreDev27(para->getParD(level)->distributionsAD.f[0],
                          para->getParD(level)->sendProcessNeighborADX[i].f[0],
                          para->getParD(level)->sendProcessNeighborADX[i].index,
                          para->getParD(level)->sendProcessNeighborADX[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborADXFsDH(level, i);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm.receiveNonBlocking(para->getParH(level)->recvProcessNeighborADX[i].f[0],
                            para->getParH(level)->recvProcessNeighborADX[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADX[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    //{
    //    comm.nbSendDataGPU(para->getParH(level)->sendProcessNeighborADX[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    //{
    //    comm.waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm.send(para->getParH(level)->sendProcessNeighborADX[i].f[0],
                          para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborADXFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->distributionsAD.f[0],
                          para->getParD(level)->recvProcessNeighborADX[i].f[0],
                          para->getParD(level)->recvProcessNeighborADX[i].index,
                          para->getParD(level)->recvProcessNeighborADX[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataADXGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        GetSendFsPostDev27(para->getParD(level)->distributionsAD.f[0],
                           para->getParD(level)->sendProcessNeighborADX[i].f[0],
                           para->getParD(level)->sendProcessNeighborADX[i].index,
                           para->getParD(level)->sendProcessNeighborADX[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborADXFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm.receiveNonBlocking(para->getParH(level)->recvProcessNeighborADX[i].f[0],
                            para->getParH(level)->recvProcessNeighborADX[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADX[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    //{
    //    comm.nbSendDataGPU(para->getParH(level)->sendProcessNeighborADX[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    //{
    //    comm.waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm.send(para->getParH(level)->sendProcessNeighborADX[i].f[0],
                          para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborADXFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->distributionsAD.f[0],
                           para->getParD(level)->recvProcessNeighborADX[i].f[0],
                           para->getParD(level)->recvProcessNeighborADX[i].index,
                           para->getParD(level)->recvProcessNeighborADX[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataADYGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        GetSendFsPreDev27(para->getParD(level)->distributionsAD.f[0],
                          para->getParD(level)->sendProcessNeighborADY[i].f[0],
                          para->getParD(level)->sendProcessNeighborADY[i].index,
                          para->getParD(level)->sendProcessNeighborADY[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborADYFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm.receiveNonBlocking(para->getParH(level)->recvProcessNeighborADY[i].f[0],
                            para->getParH(level)->recvProcessNeighborADY[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADY[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    //{
    //    comm.nbSendDataGPU(para->getParH(level)->sendProcessNeighborADY[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    //{
    //    comm.waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm.send(para->getParH(level)->sendProcessNeighborADY[i].f[0],
                          para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborADYFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->distributionsAD.f[0],
                          para->getParD(level)->recvProcessNeighborADY[i].f[0],
                          para->getParD(level)->recvProcessNeighborADY[i].index,
                          para->getParD(level)->recvProcessNeighborADY[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataADYGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        GetSendFsPostDev27(para->getParD(level)->distributionsAD.f[0],
                           para->getParD(level)->sendProcessNeighborADY[i].f[0],
                           para->getParD(level)->sendProcessNeighborADY[i].index,
                           para->getParD(level)->sendProcessNeighborADY[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborADYFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm.receiveNonBlocking(para->getParH(level)->recvProcessNeighborADY[i].f[0],
                            para->getParH(level)->recvProcessNeighborADY[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADY[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    //{
    //    comm.nbSendDataGPU(para->getParH(level)->sendProcessNeighborADY[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    //{
    //    comm.waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm.send(para->getParH(level)->sendProcessNeighborADY[i].f[0],
                          para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborADYFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->distributionsAD.f[0],
                           para->getParD(level)->recvProcessNeighborADY[i].f[0],
                           para->getParD(level)->recvProcessNeighborADY[i].index,
                           para->getParD(level)->recvProcessNeighborADY[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataADZGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        GetSendFsPreDev27(para->getParD(level)->distributionsAD.f[0],
                          para->getParD(level)->sendProcessNeighborADZ[i].f[0],
                          para->getParD(level)->sendProcessNeighborADZ[i].index,
                          para->getParD(level)->sendProcessNeighborADZ[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborADZFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm.receiveNonBlocking(para->getParH(level)->recvProcessNeighborADZ[i].f[0],
                            para->getParH(level)->recvProcessNeighborADZ[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADZ[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    //{
    //    comm.nbSendDataGPU(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    //{
    //    comm.waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm.send(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
                          para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborADZFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->distributionsAD.f[0],
                          para->getParD(level)->recvProcessNeighborADZ[i].f[0],
                          para->getParD(level)->recvProcessNeighborADZ[i].index,
                          para->getParD(level)->recvProcessNeighborADZ[i].numberOfNodes,
                          para->getParD(level)->neighborX, 
                          para->getParD(level)->neighborY, 
                          para->getParD(level)->neighborZ,
                          para->getParD(level)->numberOfNodes, 
                          para->getParD(level)->isEvenTimestep,
                          para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataADZGPU27(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        GetSendFsPostDev27(para->getParD(level)->distributionsAD.f[0],
                           para->getParD(level)->sendProcessNeighborADZ[i].f[0],
                           para->getParD(level)->sendProcessNeighborADZ[i].index,
                           para->getParD(level)->sendProcessNeighborADZ[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborADZFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm.receiveNonBlocking(para->getParH(level)->recvProcessNeighborADZ[i].f[0],
                            para->getParH(level)->recvProcessNeighborADZ[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADZ[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    //{
    //    comm.nbSendDataGPU(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    //{
    //    comm.waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm.send(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
                          para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborADZFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->distributionsAD.f[0],
                           para->getParD(level)->recvProcessNeighborADZ[i].f[0],
                           para->getParD(level)->recvProcessNeighborADZ[i].index,
                           para->getParD(level)->recvProcessNeighborADZ[i].numberOfNodes,
                           para->getParD(level)->neighborX, 
                           para->getParD(level)->neighborY, 
                           para->getParD(level)->neighborZ,
                           para->getParD(level)->numberOfNodes, 
                           para->getParD(level)->isEvenTimestep,
                           para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

















































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3D domain decomposition F3 - K18/K20
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataF3XGPU(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        getSendGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->sendProcessNeighborF3X[i].g[0],
            para->getParD(level)->sendProcessNeighborF3X[i].index,
            para->getParD(level)->sendProcessNeighborF3X[i].numberOfNodes,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep,
            para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborF3XFsDH(level, i);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm.receiveNonBlocking(
            para->getParH(level)->recvProcessNeighborF3X[i].g[0],
            para->getParH(level)->recvProcessNeighborF3X[i].numberOfGs,
            para->getParH(level)->recvProcessNeighborF3X[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm.send(
            para->getParH(level)->sendProcessNeighborF3X[i].g[0],
            para->getParH(level)->sendProcessNeighborF3X[i].numberOfGs,
            para->getParH(level)->sendProcessNeighborF3X[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborF3XFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        setRecvGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->recvProcessNeighborF3X[i].g[0],
            para->getParD(level)->recvProcessNeighborF3X[i].index,
            para->getParD(level)->recvProcessNeighborF3X[i].numberOfNodes,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep,
            para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataF3YGPU(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        getSendGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->sendProcessNeighborF3Y[i].g[0],
            para->getParD(level)->sendProcessNeighborF3Y[i].index,
            para->getParD(level)->sendProcessNeighborF3Y[i].numberOfNodes,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep,
            para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborF3YFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm.receiveNonBlocking(
            para->getParH(level)->recvProcessNeighborF3Y[i].g[0],
            para->getParH(level)->recvProcessNeighborF3Y[i].numberOfGs,
            para->getParH(level)->recvProcessNeighborF3Y[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm.send(
            para->getParH(level)->sendProcessNeighborF3Y[i].g[0],
            para->getParH(level)->sendProcessNeighborF3Y[i].numberOfGs,
            para->getParH(level)->sendProcessNeighborF3Y[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborF3YFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        setRecvGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->recvProcessNeighborF3Y[i].g[0],
            para->getParD(level)->recvProcessNeighborF3Y[i].index,
            para->getParD(level)->recvProcessNeighborF3Y[i].numberOfNodes,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep,
            para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataF3ZGPU(Parameter* para, vf::parallel::Communicator& comm, CudaMemoryManager* cudaMemoryManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        getSendGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->sendProcessNeighborF3Z[i].g[0],
            para->getParD(level)->sendProcessNeighborF3Z[i].index,
            para->getParD(level)->sendProcessNeighborF3Z[i].numberOfNodes,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep,
            para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaMemoryManager->cudaCopyProcessNeighborF3ZFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm.receiveNonBlocking(
            para->getParH(level)->recvProcessNeighborF3Z[i].g[0],
            para->getParH(level)->recvProcessNeighborF3Z[i].numberOfGs,
            para->getParH(level)->recvProcessNeighborF3Z[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm.send(
            para->getParH(level)->sendProcessNeighborF3Z[i].g[0],
            para->getParH(level)->sendProcessNeighborF3Z[i].numberOfGs,
            para->getParH(level)->sendProcessNeighborF3Z[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    comm.waitAll();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    {
        comm.resetRequests();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaMemoryManager->cudaCopyProcessNeighborF3ZFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        setRecvGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->recvProcessNeighborF3Z[i].g[0],
            para->getParD(level)->recvProcessNeighborF3Z[i].index,
            para->getParD(level)->recvProcessNeighborF3Z[i].numberOfNodes,
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->isEvenTimestep,
            para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////