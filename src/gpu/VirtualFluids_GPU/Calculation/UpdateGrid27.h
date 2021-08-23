#ifndef UPDATEGRID27_H
#define UPDATEGRID27_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"
#include "Communication/Communicator.h"
#include "Calculation/PorousMedia.h"

class Kernel;

extern "C" void updateGrid27(Parameter* para, 
                             vf::gpu::Communicator* comm, 
                             CudaMemoryManager* cudaManager, 
                             std::vector<std::shared_ptr<PorousMedia>>& pm, 
                             int level,
                             unsigned int t, 
                             std::vector < SPtr< Kernel>>& kernels);

extern "C" void collision(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,  std::vector<SPtr<Kernel>> &kernels);

extern "C" void collisionUsingIndex(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,  std::vector<SPtr<Kernel>> &kernels, uint *fluidNodeIndices = nullptr, uint numberOfFluidNodes = 0, int stream = -1);

extern "C" void collisionPorousMedia(Parameter* para, std::vector<std::shared_ptr<PorousMedia>>& pm, int level);

extern "C" void collisionAdvectionDiffusion(Parameter* para, int level);

extern "C" void prepareExchangeMultiGPU(Parameter *para, int level, int streamIndex);

extern "C" void exchangeMultiGPU(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                                 int level, int streamIndex);

extern "C" void postCollisionBC(Parameter* para, int level, unsigned int t);

extern "C" void swapBetweenEvenAndOddTimestep(Parameter* para, int level);

extern "C" void calcMacroscopicQuantities(Parameter* para, int level);

extern "C" void preCollisionBC(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t);

extern "C" void fineToCoarse(Parameter* para, int level);

extern "C" void coarseToFine(Parameter* para, int level);

#endif
