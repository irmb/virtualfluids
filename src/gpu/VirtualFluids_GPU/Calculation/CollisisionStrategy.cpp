#include "CollisionStrategy.h"
#include "Parameter/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include "logger/Logger.h"

std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t)>
getFunctionForCollisionAndExchange(const bool useStreams, const int numberOfMpiProcesses,
                                   const bool kernelNeedsFluidNodeIndicesToRun)
{
    VF_LOG_INFO("Function used for collisionAndExchange: ");

    if (useStreams && numberOfMpiProcesses > 1 && kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("collisionAndExchange_streams()");
        return CollisionAndExchange_streams();

    } else if (useStreams && !kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("Cuda Streams can only be used with kernels which run using fluidNodesIndices.");

    } else if (useStreams && numberOfMpiProcesses <= 1) {
        VF_LOG_INFO("Cuda Streams can only be used with multiple MPI processes.");

    } else if (!useStreams && kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("collisionAndExchange_noStreams_indexKernel()");
        return CollisionAndExchange_noStreams_indexKernel();

    } else if (!useStreams && !kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("collisionAndExchange_noStreams_oldKernel()");
        return CollisionAndExchange_noStreams_oldKernel();
    }

    throw std::runtime_error("Invalid Configuration for collision and exchange");
    return nullptr;
}

void CollisionAndExchange_noStreams_indexKernel::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level,
                                                            unsigned int t)
{
    updateGrid->collisionUsingIndex(level, t, para->getParD(level)->fluidNodeIndices,
                                    para->getParD(level)->numberOfFluidNodes, -1);
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, false);
}

void CollisionAndExchange_noStreams_oldKernel::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level,
                                                          unsigned int t)
{
    updateGrid->collisionAllNodes(level, t);
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, false);
}

void CollisionAndExchange_streams::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t)
{
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex = para->getStreamManager()->getBulkStreamIndex();

    // launch border kernel
    updateGrid->collisionUsingIndex(level, t, para->getParD(level)->fluidNodeIndicesBorder,
                                    para->getParD(level)->numberOfFluidNodesBorder, borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
    updateGrid->prepareExchangeMultiGPU(level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernel
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    updateGrid->collisionUsingIndex(level, t, para->getParD(level)->fluidNodeIndices,
                                    para->getParD(level)->numberOfFluidNodes, bulkStreamIndex);

    updateGrid->exchangeMultiGPU(level, borderStreamIndex);
}
