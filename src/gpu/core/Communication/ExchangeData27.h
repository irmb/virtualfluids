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

#ifndef EXCHANGEDATA27_H
#define EXCHANGEDATA27_H

#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include "LBM/LB.h"
#include "Parameter/CudaStreamManager.h"
#include "Parameter/Parameter.h"

namespace vf::parallel
{
class Communicator;
}

//! \file ExchangeData27.h
//! \ingroup GPU
//! \author Martin Schoenherr, Anna Wellmann
//! \brief Routines for data exchange when running simulations on multiple GPUs

//////////////////////////////////////////////////////////////////////////
// 1D domain decomposition
void exchangePreCollDataGPU27(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager, 
                                         int level);
void exchangePostCollDataGPU27(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager, 
                                          int level);
//////////////////////////////////////////////////////////////////////////
// 3D domain decomposition

// functions used for all directions

//! \brief Collect the send nodes in a buffer on the gpu
void collectNodesInSendBufferGPU(Parameter *para, int level, CudaStreamIndex streamIndex,
                                 std::vector<ProcessNeighbor27> *sendProcessNeighbor,
                                 unsigned int numberOfSendProcessNeighbors);
//! \brief Distribute the receive nodes from the buffer on the gpu
void scatterNodesFromRecvBufferGPU(Parameter *para, int level, CudaStreamIndex streamIndex,
                                   std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                   unsigned int numberOfRecvProcessNeighbors);
//! \brief Copy nodes which are part of the communication in multiple directions
//! \details The nodes are copied from the receive buffer in one direction to the send buffer in another direction. The
//! copy operation is conducted on the cpu. 
//! See [master thesis of Anna Wellmann (p. 56f: "Communication Hiding bei
//! der Verwendung eines uniformen Simulationsgitters")]
//! \param edgeNodes determines from where to where the nodes are
//! copied 
//! \param recvProcessNeighborHost is a reference to the receive buffer on the host, nodes are copied from here
//! \param sendProcessNeighborHost is a reference to the send buffer on the host, nodes are copied to here
void copyEdgeNodes(std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes,
                              std::vector<ProcessNeighbor27> &recvProcessNeighborHost,
                              std::vector<ProcessNeighbor27> &sendProcessNeighborHost);

//////////////////////////////////////////////////////////////////////////
// x

//! \brief Collect the send nodes for communication in the x direction in a buffer on the gpu
//! \details Needed to exchange all nodes, used in the communication after collision step
void prepareExchangeCollDataXGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex);
//! \brief Collect the send nodes for communication in the x direction in a buffer on the gpu
//! \details Only exchange nodes which are part of the interpolation process on refined grids. This function is used in
//! the exchange which takes place after the interpolation fine to coarse and before the interpolation coarse to fine.
//! See [master thesis of Anna Wellmann]
void prepareExchangeCollDataXGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex);
//! \brief Exchange routine in x direction for simulations on multiple gpus
//! \details Send and receive the nodes from the communication buffers on the gpus.
//! \param Communicator is needed for the communication between the processes with mpi
//! \param CudaMemoryManager is needed for moving the data between host and device
//! \param sendProcessNeighborDev, recvProcessNeighborDev, sendProcessNeighborHost, recvProcessNeighborHost are pointers
//! to the send and receive arrays, both on the device and the host
void exchangeCollDataXGPU27(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                       int level, CudaStreamIndex streamIndex,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborHost);
//! \brief Calls exchangeCollDataXGPU27() for exchanging all nodes
//! \details Used in the communication after collision step
void exchangeCollDataXGPU27AllNodes(Parameter *para, vf::parallel::Communicator& comm,
                                               CudaMemoryManager *cudaMemoryManager, int level, CudaStreamIndex streamIndex);
//! \brief Calls exchangeCollDataXGPU27() for exchanging the nodes, which are part of the communication between the two
//! interpolation processes on refined grids 
//! \details Only exchange nodes which are part of the interpolation process on
//! refined grids. This function is used in the exchange which takes place after the interpolation fine to coarse and
//! before the interpolation coarse to fine. See [master thesis of Anna Wellmann]
void exchangeCollDataXGPU27AfterFtoC(Parameter *para, vf::parallel::Communicator& comm,
                                                CudaMemoryManager *cudaMemoryManager, int level, CudaStreamIndex streamIndex);
//! \brief Distribute the receive nodes (x direction) from the buffer on the gpu
//! \details Needed to exchange all nodes, used in the communication after collision step
void scatterNodesFromRecvBufferXGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex);
//! \brief Distribute the receive nodes (x direction) from the buffer on the gpu
//! \details Only exchange nodes which are part of the interpolation process on refined grids. This function is used in
//! the exchange which takes place after the interpolation fine to coarse and before the interpolation coarse to fine.
//! See [master thesis of Anna Wellmann]
void scatterNodesFromRecvBufferXGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex);

//////////////////////////////////////////////////////////////////////////
// y

void prepareExchangeCollDataYGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex);
void prepareExchangeCollDataYGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex);

void exchangeCollDataYGPU27(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                       int level,CudaStreamIndex streamIndex,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborHos);
void exchangeCollDataYGPU27AllNodes(Parameter *para, vf::parallel::Communicator& comm,
                                               CudaMemoryManager *cudaMemoryManager, int level, CudaStreamIndex streamIndex);
void exchangeCollDataYGPU27AfterFtoC(Parameter *para, vf::parallel::Communicator& comm,
                                                CudaMemoryManager *cudaMemoryManager, int level, CudaStreamIndex streamIndex);
void scatterNodesFromRecvBufferYGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex);
void scatterNodesFromRecvBufferYGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex);

// z
void prepareExchangeCollDataZGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex);
void prepareExchangeCollDataZGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex);

void exchangeCollDataZGPU27(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                       int level, CudaStreamIndex streamIndex,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborHost);
void exchangeCollDataZGPU27AllNodes(Parameter *para, vf::parallel::Communicator& comm,
                                               CudaMemoryManager *cudaMemoryManager, int level, CudaStreamIndex streamIndex);
void exchangeCollDataZGPU27AfterFtoC(Parameter *para, vf::parallel::Communicator& comm,
                                                CudaMemoryManager *cudaMemoryManager, int level, CudaStreamIndex streamIndex);

void scatterNodesFromRecvBufferZGPU27AllNodes(Parameter *para, int level, CudaStreamIndex streamIndex);
void scatterNodesFromRecvBufferZGPU27AfterFtoC(Parameter *para, int level, CudaStreamIndex streamIndex);

//////////////////////////////////////////////////////////////////////////
// 3D domain decomposition convection diffusion
void exchangePreCollDataADXGPU27(Parameter *para, vf::parallel::Communicator& comm,
                                            CudaMemoryManager *cudaMemoryManager, int level);
void exchangePreCollDataADYGPU27(Parameter *para, vf::parallel::Communicator& comm,
                                            CudaMemoryManager *cudaMemoryManager, int level);
void exchangePreCollDataADZGPU27(Parameter *para, vf::parallel::Communicator& comm,
                                            CudaMemoryManager *cudaMemoryManager, int level);
void exchangePostCollDataADXGPU27(Parameter *para, vf::parallel::Communicator& comm,
                                             CudaMemoryManager *cudaMemoryManager, int level);
void exchangePostCollDataADYGPU27(Parameter *para, vf::parallel::Communicator& comm,
                                             CudaMemoryManager *cudaMemoryManager, int level);
void exchangePostCollDataADZGPU27(Parameter *para, vf::parallel::Communicator& comm,
                                             CudaMemoryManager *cudaMemoryManager, int level);
//////////////////////////////////////////////////////////////////////////
// 3D domain decomposition F3 - K18/K20
void exchangeCollDataF3XGPU(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                       int level);
void exchangeCollDataF3YGPU(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                       int level);
void exchangeCollDataF3ZGPU(Parameter *para, vf::parallel::Communicator& comm, CudaMemoryManager *cudaMemoryManager,
                                       int level);

#endif
