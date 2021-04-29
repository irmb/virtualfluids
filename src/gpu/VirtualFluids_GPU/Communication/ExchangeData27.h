#ifndef EXCHANGEDATA27_H
#define EXCHANGEDATA27_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"
#include "GPU/CudaMemoryManager.h"

//////////////////////////////////////////////////////////////////////////
//1D domain decomposition
extern "C" void exchangePreCollDataGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
//3D domain decomposition
extern "C" void exchangePreCollDataXGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePreCollDataYGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePreCollDataZGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataXGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataYGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataZGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
//3D domain decomposition convection diffusion
extern "C" void exchangePreCollDataADXGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePreCollDataADYGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePreCollDataADZGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataADXGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataADYGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangePostCollDataADZGPU27(Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
//3D domain decomposition F3 - K18/K20
extern "C" void exchangeCollDataF3XGPU( Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangeCollDataF3YGPU( Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
extern "C" void exchangeCollDataF3ZGPU( Parameter* para, VF::GPU::Communicator* comm, CudaMemoryManager* cudaManager, int level);
//////////////////////////////////////////////////////////////////////////
extern "C" void barrierGPU(VF::GPU::Communicator* comm);
//////////////////////////////////////////////////////////////////////////

#endif
