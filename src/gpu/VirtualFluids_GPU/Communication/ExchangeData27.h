#ifndef EXCHANGEDATA27_H
#define EXCHANGEDATA27_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"
#include "GPU/CudaMemoryManager.h"

//////////////////////////////////////////////////////////////////////////
//1D domain decomposition
extern "C" void exchangePreCollDataGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
//3D domain decomposition

extern "C" void prepareExchangeCollDataXGPU27(Parameter *para, int level, int streamIndex);
extern "C" void exchangeCollDataXGPU27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                                       int level, int streamIndex);
extern "C" void prepareExchangeCollDataYGPU27(Parameter *para, int level, int streamIndex);
extern "C" void exchangeCollDataYGPU27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                                       int level, int streamIndex);
extern "C" void prepareExchangeCollDataZGPU27(Parameter *para, int level, int streamIndex);
extern "C" void exchangeCollDataZGPU27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                                       int level, int streamIndex);
//////////////////////////////////////////////////////////////////////////
//3D domain decomposition convection diffusion
extern "C" void exchangePreCollDataADXGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePreCollDataADYGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePreCollDataADZGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataADXGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataADYGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataADZGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
//3D domain decomposition F3 - K18/K20
extern "C" void exchangeCollDataF3XGPU( Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangeCollDataF3YGPU( Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangeCollDataF3ZGPU( Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
extern "C" void barrierGPU(vf::gpu::Communicator* comm);
//////////////////////////////////////////////////////////////////////////

#endif
