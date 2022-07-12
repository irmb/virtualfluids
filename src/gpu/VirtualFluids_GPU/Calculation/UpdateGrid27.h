#ifndef UPDATEGRID27_H
#define UPDATEGRID27_H

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"
#include "Communication/Communicator.h"
#include "Calculation/PorousMedia.h"

class BCKernelManager;
class ADKernelManager;
class GridScalingKernelManager;
class Kernel;
class BoundaryConditionFactory;

class UpdateGrid27
{
public:
    UpdateGrid27(SPtr<Parameter> para, vf::gpu::Communicator &comm, SPtr<CudaMemoryManager> cudaMemoryManager,
                 std::vector<std::shared_ptr<PorousMedia>> &pm, std::vector<SPtr<Kernel>> &kernels, BoundaryConditionFactory* bcFactory);
    void updateGrid(int level, unsigned int t);
    void exchangeData(int level);

private:
    void postCollisionBC(int level);
    void preCollisionBC(int level, unsigned int t);

    void collision(int level, unsigned int t);
    void collisionUsingIndex(int level, unsigned int t, uint *fluidNodeIndices = nullptr, uint numberOfFluidNodes = 0, int stream = -1);
    void collisionAdvectionDiffusion(int level);
    void collisionPorousMedia(int level);

    void fineToCoarse(int level, uint *iCellFCC, uint *iCellFCF, uint k_FC, int streamIndex);
    void coarseToFine(int level, uint *iCellCFC, uint *iCellCFF, uint k_CF, OffCF &offCF, int streamIndex);

    void prepareExchangeMultiGPU(int level, int streamIndex);
    void prepareExchangeMultiGPUAfterFtoC(int level, int streamIndex);

    void exchangeMultiGPU(int level, int streamIndex);
    void exchangeMultiGPUAfterFtoC(int level, int streamIndex);
    void exchangeMultiGPU_noStreams_withPrepare(int level, bool useReducedComm);

    void swapBetweenEvenAndOddTimestep(int level);

    void calcMacroscopicQuantities(int level);
    void calcTurbulentViscosity(int level);
    void interactWithActuators(int level, unsigned int t);
    void interactWithProbes(int level, unsigned int t);

private:
    typedef void (UpdateGrid27::*collisionAndExchangeFun)(int level, unsigned int t);
    typedef void (UpdateGrid27::*refinementAndExchangeFun)(int level);
    collisionAndExchangeFun collisionAndExchange   = nullptr;
    refinementAndExchangeFun refinementAndExchange  = nullptr;

    void chooseFunctionForCollisionAndExchange();
    void chooseFunctionForRefinementAndExchange();

    // functions for collision and exchange
    void collisionAndExchange_noStreams_indexKernel(int level, unsigned int t);
    void collisionAndExchange_noStreams_oldKernel(int level, unsigned int t);
    void collisionAndExchange_streams(int level, unsigned int t);

    // functions for refinement and exchange
    void refinementAndExchange_streams_onlyExchangeInterface(int level);
    void refinementAndExchange_streams_completeExchange(int level);
    void refinementAndExchange_noStreams_onlyExchangeInterface(int level);
    void refinementAndExchange_noStreams_completeExchange(int level);
    void refinementAndExchange_noRefinementAndExchange(int level);
    void refinementAndExchange_noExchange(int level);

private:
    SPtr<Parameter> para;
    vf::gpu::Communicator& comm;
    SPtr<CudaMemoryManager> cudaMemoryManager;
    std::vector<std::shared_ptr<PorousMedia>> pm;
    std::vector<SPtr<Kernel>> kernels;
    //! \property lbKernelManager is a shared pointer to an object of LBKernelManager
    std::shared_ptr<BCKernelManager> bcKernelManager;
    //! \property adKernelManager is a shared pointer to an object of ADKernelManager
    std::shared_ptr<ADKernelManager> adKernelManager;
    //! \property gridScalingKernelManager is a shared pointer to an object of GridScalingKernelManager
    std::shared_ptr<GridScalingKernelManager> gridScalingKernelManager;
};

#endif
