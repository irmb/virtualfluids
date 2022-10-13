#include "RefinementStrategy.h"
#include "Parameter/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include "logger/Logger.h"

std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level)>
    getFunctionForRefinementAndExchange(const bool useStreams, const int numberOfMpiProcesses, const int maxLevel,
                                        const bool useReducedCommunicationAfterFtoC) noexcept
{
    VF_LOG_INFO("Function used for refinementAndExchange: ");
    if (maxLevel == 0) {
        VF_LOG_INFO("only one level - no function needed.");
        return NoRefinement();

    } else if (numberOfMpiProcesses == 1) {
        VF_LOG_INFO("only one process - no exchange needed: Refinement_noExchange()");
        return Refinement_noExchange();

    } else if (numberOfMpiProcesses > 1 && useStreams && useReducedCommunicationAfterFtoC) {
        VF_LOG_INFO("RefinementAndExchange_streams_exchangeInterface()");
        return RefinementAndExchange_streams_exchangeInterface();

    } else if(numberOfMpiProcesses > 1 && useStreams && !useReducedCommunicationAfterFtoC){
        VF_LOG_INFO("refinementAndExchange_streams_completeExchange()");
        return RefinementAndExchange_streams_exchangeAllNodes();

    } else if (numberOfMpiProcesses > 1 && !useStreams && useReducedCommunicationAfterFtoC) {
        VF_LOG_INFO("RefinementAndExchange_noStreams_exchangeInterface()");
        return RefinementAndExchange_noStreams_exchangeInterface();

    } else {
        VF_LOG_INFO("RefinementAndExchange_noStreams_exchangeAllNodes()");
        return RefinementAndExchange_noStreams_exchangeAllNodes();
    }
}

void NoRefinement::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level){}

void RefinementAndExchange_streams_exchangeInterface::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. Interpolation fine to coarse for nodes which are at the border of the gpus/processes
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->intFCBorder, para->getParD(level)->offFC, CudaStreamIndex::Border);

    //! 2. prepare the exchange between gpus (collect the send nodes for communication in a buffer on the gpu) and trigger bulk kernel execution when finished
    //!
    updateGrid->prepareExchangeMultiGPUAfterFtoC(level, CudaStreamIndex::Border);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(CudaStreamIndex::Border);

    //! 3. launch the bulk kernels for both interpolation processes (fine to coarse and coarse to fine)
    //!
    para->getStreamManager()->waitOnStartBulkKernelEvent(CudaStreamIndex::Bulk);
    updateGrid->fineToCoarse(level, &para->getParD(level)->intFCBulk, para->getParD(level)->offFCBulk, CudaStreamIndex::Border);
    updateGrid->coarseToFine(level, &para->getParD(level)->intCFBulk, para->getParD(level)->offCFBulk, CudaStreamIndex::Border);

    //! 4. exchange information between GPUs (only nodes which are part of the interpolation)
    //!
    updateGrid->exchangeMultiGPUAfterFtoC(level, CudaStreamIndex::Border);

    // 5. interpolation fine to coarse for nodes which are at the border of the gpus/processes
    //!
    updateGrid->coarseToFine(level, &para->getParD(level)->intCFBorder, para->getParD(level)->offCF, CudaStreamIndex::Border);

    cudaDeviceSynchronize();
}

void RefinementAndExchange_streams_exchangeAllNodes::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse for nodes which are at the border of the gpus/processes
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->intFCBorder, para->getParD(level)->offFC, CudaStreamIndex::Border);

    //! 2. prepare the exchange between gpus (collect the send nodes for communication in a buffer on the gpu) and trigger bulk kernel execution when finished
    //!
    updateGrid->prepareExchangeMultiGPU(level, CudaStreamIndex::Border);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(CudaStreamIndex::Border);

    //! 3. launch the bulk kernels for both interpolation processes (fine to coarse and coarse to fine)
    //!
    para->getStreamManager()->waitOnStartBulkKernelEvent(CudaStreamIndex::Bulk);
    updateGrid->fineToCoarse(level, &para->getParD(level)->intFCBulk, para->getParD(level)->offFCBulk, CudaStreamIndex::Border);
    updateGrid->coarseToFine(level, &para->getParD(level)->intCFBulk, para->getParD(level)->offCFBulk, CudaStreamIndex::Border);

    //! 4. exchange information between GPUs (all nodes)
    //!
    updateGrid->exchangeMultiGPU(level, CudaStreamIndex::Border);

    // 5. interpolation fine to coarse for nodes which are at the border of the gpus/processes
    //!
    updateGrid->coarseToFine(level, &para->getParD(level)->intCFBorder, para->getParD(level)->offCF, CudaStreamIndex::Border);

    cudaDeviceSynchronize();
}

void RefinementAndExchange_noStreams_exchangeInterface::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->intFC, para->getParD(level)->offFC, CudaStreamIndex::Legacy);

    //! 2. exchange information between GPUs (only nodes which are part of the interpolation)
    //!
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, true);

    //! 3. interpolation coarse to fine
    updateGrid->coarseToFine(level, &para->getParD(level)->intCF, para->getParD(level)->offCF, CudaStreamIndex::Legacy);
}

void RefinementAndExchange_noStreams_exchangeAllNodes::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->intFC, para->getParD(level)->offFC, CudaStreamIndex::Legacy);

    //! 2. exchange information between GPUs (all nodes)
    //!
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, false);

    //! 3. interpolation coarse to fine
    updateGrid->coarseToFine(level, &para->getParD(level)->intCF, para->getParD(level)->offCF, CudaStreamIndex::Legacy);
}

void Refinement_noExchange::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->intFC, para->getParD(level)->offFC, CudaStreamIndex::Legacy);
    //! 2. interpolation coarse to fine
    updateGrid->coarseToFine(level, &para->getParD(level)->intCF, para->getParD(level)->offCF, CudaStreamIndex::Legacy);
}
