#ifndef UPDATEGRID27_H
#define UPDATEGRID27_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"
#include "Communication/Communicator.h"
#include "Calculation/PorousMedia.h"

class Kernel;

class UpdateGrid27
{
public:
    UpdateGrid27(Parameter *para);
    ~UpdateGrid27();
    void updateGrid27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                      std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,
                      std::vector<SPtr<Kernel>> &kernels);

private:
    UpdateGrid27();
    std::function<void(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,
                       std::vector<SPtr<Kernel>> &kernels, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager)>
        collisionAndExchange = nullptr;
    std::function<void(Parameter *para, int level, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager)>
        refinementAndExchange = nullptr;

    void chooseFunctionForCollisionAndExchange(Parameter *para);
    void chooseFunctionForRefinementAndExchange(Parameter *para);


};

extern "C" void collisionAndExchange_noStreams_indexKernel(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm,
                                                       int level, unsigned int t, std::vector<SPtr<Kernel>> &kernels,
                                                       vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager);

extern "C" void collisionAndExchange_noStreams_oldKernel(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm,
                                                       int level, unsigned int t, std::vector<SPtr<Kernel>> &kernels,
                                                       vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager);

extern "C" void collisionAndExchange_streams(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level,
                                     unsigned int t, std::vector<SPtr<Kernel>> &kernels, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager);

extern "C" void collision(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,  std::vector<SPtr<Kernel>> &kernels);

extern "C" void collisionUsingIndex(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,  std::vector<SPtr<Kernel>> &kernels, uint *fluidNodeIndices = nullptr, uint numberOfFluidNodes = 0, int stream = -1);

extern "C" void collisionPorousMedia(Parameter* para, std::vector<std::shared_ptr<PorousMedia>>& pm, int level);

extern "C" void collisionAdvectionDiffusion(Parameter* para, int level);

extern "C" void prepareExchangeMultiGPU(Parameter *para, int level, int streamIndex);
extern "C" void prepareExchangeMultiGPUAfterFtoC(Parameter *para, int level, int streamIndex);

extern "C" void exchangeMultiGPU(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                                 int level, int streamIndex);
extern "C" void exchangeMultiGPUAfterFtoC(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                                 int level, int streamIndex);

extern "C" void postCollisionBC(Parameter* para, int level, unsigned int t);

extern "C" void swapBetweenEvenAndOddTimestep(Parameter* para, int level);

extern "C" void calcMacroscopicQuantities(Parameter* para, int level);

extern "C" void preCollisionBC(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t);

extern "C" void fineToCoarse(Parameter* para, int level);
extern "C" void fineToCoarseWithStream(Parameter *para, int level, uint *iCellFCC, uint *iCellFCF, uint k_FC, int streamIndex);

extern "C" void coarseToFine(Parameter* para, int level);

extern "C" void refinementAndExchange_streams(Parameter *para, int level, vf::gpu::Communicator *comm,
                                              CudaMemoryManager *cudaManager);
extern "C" void refinementAndExchange_noStreams_onlyExchangeInterface(Parameter *para, int level,
                                                                      vf::gpu::Communicator *comm,
                                                                      CudaMemoryManager *cudaManager);
extern "C" void refinementAndExchange_noStreams_completeExchange(Parameter *para, int level,
                                                                 vf::gpu::Communicator *comm,
                                                                 CudaMemoryManager *cudaManager);
extern "C" void refinementAndExchange_noExchange(Parameter *para, int level, vf::gpu::Communicator *comm,
                                                  CudaMemoryManager *cudaManager);

#endif
