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
class TurbulenceModelManager;
class GridScalingKernelManager;
class Kernel;
class BoundaryConditionFactory;
class TurbulenceModelFactory;

class UpdateGrid27;
using CollisionStrategy = std::function<void (UpdateGrid27* updateGrid, Parameter* para, int level, unsigned int t)>;
using RefinementStrategy = std::function<void (UpdateGrid27* updateGrid, Parameter* para, int level)>;


class UpdateGrid27
{
public:
    UpdateGrid27(SPtr<Parameter> para, vf::gpu::Communicator &comm, SPtr<CudaMemoryManager> cudaMemoryManager,
                 std::vector<std::shared_ptr<PorousMedia>> &pm, std::vector<SPtr<Kernel>> &kernels, BoundaryConditionFactory* bcFactory, SPtr<TurbulenceModelFactory> tmFactory);
    void updateGrid(int level, unsigned int t);
    void exchangeData(int level);

private:
    void collisionAllNodes(int level, unsigned int t);
    void collisionUsingIndices(int level, unsigned int t, uint *fluidNodeIndices = nullptr, uint numberOfFluidNodes = 0, int stream = -1);
    void collisionAdvectionDiffusion(int level);

    void postCollisionBC(int level);
    void preCollisionBC(int level, unsigned int t);
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
    CollisionStrategy collision;
    friend class CollisionAndExchange_noStreams_indexKernel;
    friend class CollisionAndExchange_noStreams_oldKernel;
    friend class CollisionAndExchange_streams;

    RefinementStrategy refinement;
    friend class RefinementAndExchange_streams_exchangeInterface;
    friend class RefinementAndExchange_streams_exchangeAllNodes;
    friend class RefinementAndExchange_noStreams_exchangeInterface;
    friend class RefinementAndExchange_noStreams_exchangeAllNodes;
    friend class Refinement_noExchange;
    friend class NoRefinement;

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
    //! \property tmFactory is a shared pointer to an object of TurbulenceModelFactory
    std::shared_ptr<TurbulenceModelFactory> tmFactory;
};

#endif
