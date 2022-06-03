#ifndef EXCHANGEDATA27_H
#define EXCHANGEDATA27_H

#include "Communication/Communicator.h"
#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include "LBM/LB.h"
#include "Parameter/Parameter.h"

//! \file ExchangeData27.h
//! \ingroup GPU
//! \author Martin Schoenherr, Anna Wellmann
//! \brief routines for data exchange when running simulations on multiple GPUs

//////////////////////////////////////////////////////////////////////////
// 1D domain decomposition
extern "C" void exchangePreCollDataGPU27(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                         int level);
extern "C" void exchangePostCollDataGPU27(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                          int level);
//////////////////////////////////////////////////////////////////////////
// 3D domain decomposition

// functions used for all directions

//! \brief collect the send nodes in a buffer on the gpu
extern "C" void collectNodesInSendBufferGPU(Parameter *para, int level, int streamIndex,
                                            std::vector<ProcessNeighbor27> *sendProcessNeighbor,
                                            unsigned int numberOfSendProcessNeighbors);
//! \brief distribute the receive nodes from the buffer on the gpu
extern "C" void scatterNodesFromRecvBufferGPU(Parameter *para, int level, int streamIndex,
                                              std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                              unsigned int numberOfRecvProcessNeighbors);
//! \brief copy nodes which are part of the communication in multiple directions
//! \details The nodes are copied from the receive buffer in one direction to the send buffer in another direction. The
//! copy operation is conducted on the cpu. 
//! \ref see master thesis of Anna Wellmann (p. 56f: "Communication Hiding bei
//! der Verwendung eines uniformen Simulationsgitters") 
//! \param edgeNodes determines from where to where the nodes are
//! copied 
//! \param recvProcessNeighborHost is a reference to the receive buffer on the host, nodes are copied from here
//! \param sendProcessNeighborHost is a reference to the send buffer on the host, nodes are copied to here
extern "C" void copyEdgeNodes(std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes,
                              std::vector<ProcessNeighbor27> &recvProcessNeighborHost,
                              std::vector<ProcessNeighbor27> &sendProcessNeighborHost);

//////////////////////////////////////////////////////////////////////////
// x

//! \brief collect the send nodes for communication in the x direction in a buffer on the gpu
//! \details needed to exchange all nodes, used in the communication after collision step
extern "C" void prepareExchangeCollDataXGPU27AllNodes(Parameter *para, int level, int streamIndex);
//! \brief collect the send nodes for communication in the x direction in a buffer on the gpu
//! \details Only exchange nodes which are part of the interpolation process on refined grids. This function is used in
//! the exchange which takes place after the interpolation fine to coarse and before the interpolation coarse to fine.
//! \ref see master thesis of Anna Wellmann
extern "C" void prepareExchangeCollDataXGPU27AfterFtoC(Parameter *para, int level, int streamIndex);
//! \brief exchange routine in x direction for simulations on multiple gpus
//! \details send and receive the nodes from the communication buffers on the gpus
//! \param Communicator is needed for the communication between the processes with mpi
//! \param CudaMemoryManager is needed for moving the data between host and device
//! \param streamIndex is the index of a CUDA Stream, which is needed for communication hiding
//! \param sendProcessNeighborDev, recvProcessNeighborDev, sendProcessNeighborHost, recvProcessNeighborHost are pointers
//! to the send and receive arrays, both on the device and the host
extern "C" void exchangeCollDataXGPU27(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                       int level, int streamIndex,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborHost);
//! \brief calls exchangeCollDataXGPU27() for exchanging all nodes
//! \details used in the communication after collision step
extern "C" void exchangeCollDataXGPU27AllNodes(Parameter *para, vf::gpu::Communicator &comm,
                                               CudaMemoryManager *cudaManager, int level, int streamIndex);
//! \brief calls exchangeCollDataXGPU27() for exchanging the nodes, which are part of the communication between the two
//! interpolation processes on refined grids \details Only exchange nodes which are part of the interpolation process on
//! refined grids. This function is used in the exchange which takes place after the interpolation fine to coarse and
//! before the interpolation coarse to fine. \ref see master thesis of Anna Wellmann
extern "C" void exchangeCollDataXGPU27AfterFtoC(Parameter *para, vf::gpu::Communicator &comm,
                                                CudaMemoryManager *cudaManager, int level, int streamIndex);
//! \brief distribute the receive nodes (x direction) from the buffer on the gpu
//! \details needed to exchange all nodes, used in the communication after collision step
extern "C" void scatterNodesFromRecvBufferXGPU27AllNodes(Parameter *para, int level, int streamIndex);
//! \brief distribute the receive nodes (x direction) from the buffer on the gpu
//! \details Only exchange nodes which are part of the interpolation process on refined grids. This function is used in
//! the exchange which takes place after the interpolation fine to coarse and before the interpolation coarse to fine.
//! \ref see master thesis of Anna Wellmann
extern "C" void scatterNodesFromRecvBufferXGPU27AfterFtoC(Parameter *para, int level, int streamIndex);

//////////////////////////////////////////////////////////////////////////
// y

extern "C" void prepareExchangeCollDataYGPU27AllNodes(Parameter *para, int level, int streamIndex);
extern "C" void prepareExchangeCollDataYGPU27AfterFtoC(Parameter *para, int level, int streamIndex);

extern "C" void exchangeCollDataYGPU27(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                       int level, int streamIndex,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborHos);
extern "C" void exchangeCollDataYGPU27AllNodes(Parameter *para, vf::gpu::Communicator &comm,
                                               CudaMemoryManager *cudaManager, int level, int streamIndex);
extern "C" void exchangeCollDataYGPU27AfterFtoC(Parameter *para, vf::gpu::Communicator &comm,
                                                CudaMemoryManager *cudaManager, int level, int streamIndex);
extern "C" void scatterNodesFromRecvBufferYGPU27AllNodes(Parameter *para, int level, int streamIndex);
extern "C" void scatterNodesFromRecvBufferYGPU27AfterFtoC(Parameter *para, int level, int streamIndex);

// z
extern "C" void prepareExchangeCollDataZGPU27AllNodes(Parameter *para, int level, int streamIndex);
extern "C" void prepareExchangeCollDataZGPU27AfterFtoC(Parameter *para, int level, int streamIndex);

extern "C" void exchangeCollDataZGPU27(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                       int level, int streamIndex,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborDev,
                                       std::vector<ProcessNeighbor27> *sendProcessNeighborHost,
                                       std::vector<ProcessNeighbor27> *recvProcessNeighborHost);
extern "C" void exchangeCollDataZGPU27AllNodes(Parameter *para, vf::gpu::Communicator &comm,
                                               CudaMemoryManager *cudaManager, int level, int streamIndex);
extern "C" void exchangeCollDataZGPU27AfterFtoC(Parameter *para, vf::gpu::Communicator &comm,
                                                CudaMemoryManager *cudaManager, int level, int streamIndex);

extern "C" void scatterNodesFromRecvBufferZGPU27AllNodes(Parameter *para, int level, int streamIndex);
extern "C" void scatterNodesFromRecvBufferZGPU27AfterFtoC(Parameter *para, int level, int streamIndex);

//////////////////////////////////////////////////////////////////////////
// 3D domain decomposition convection diffusion
extern "C" void exchangePreCollDataADXGPU27(Parameter *para, vf::gpu::Communicator &comm,
                                            CudaMemoryManager *cudaManager, int level);
extern "C" void exchangePreCollDataADYGPU27(Parameter *para, vf::gpu::Communicator &comm,
                                            CudaMemoryManager *cudaManager, int level);
extern "C" void exchangePreCollDataADZGPU27(Parameter *para, vf::gpu::Communicator &comm,
                                            CudaMemoryManager *cudaManager, int level);
extern "C" void exchangePostCollDataADXGPU27(Parameter *para, vf::gpu::Communicator &comm,
                                             CudaMemoryManager *cudaManager, int level);
extern "C" void exchangePostCollDataADYGPU27(Parameter *para, vf::gpu::Communicator &comm,
                                             CudaMemoryManager *cudaManager, int level);
extern "C" void exchangePostCollDataADZGPU27(Parameter *para, vf::gpu::Communicator &comm,
                                             CudaMemoryManager *cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
// 3D domain decomposition F3 - K18/K20
extern "C" void exchangeCollDataF3XGPU(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                       int level);
extern "C" void exchangeCollDataF3YGPU(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                       int level);
extern "C" void exchangeCollDataF3ZGPU(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager,
                                       int level);
//////////////////////////////////////////////////////////////////////////
extern "C" void barrierGPU(vf::gpu::Communicator &comm);
//////////////////////////////////////////////////////////////////////////

#endif
