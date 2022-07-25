#include "RefinementStrategy.h"
#include "Parameter/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include "logger/Logger.h"

std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level)>
    getFunctionForRefinementAndExchange(const bool useStreams, const int numberOfMpiProcesses, const int maxLevel,
                                        const bool useReducedCommunicationAfterFtoC)
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
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex = para->getStreamManager()->getBulkStreamIndex();

    // fine to coarse border
    updateGrid->fineToCoarse(level, para->getParD(level)->intFCBorder.ICellFCC, para->getParD(level)->intFCBorder.ICellFCF,
                 para->getParD(level)->intFCBorder.kFC, borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
     updateGrid->prepareExchangeMultiGPUAfterFtoC(level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernels (fine to coarse and coarse to fine)
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    updateGrid->fineToCoarse(level, para->getParD(level)->intFCBulk.ICellFCC, para->getParD(level)->intFCBulk.ICellFCF,
                 para->getParD(level)->intFCBulk.kFC, bulkStreamIndex);
    updateGrid->coarseToFine(level, para->getParD(level)->intCFBulk.ICellCFC, para->getParD(level)->intCFBulk.ICellCFF,
                 para->getParD(level)->intCFBulk.kCF, para->getParD(level)->offCFBulk, bulkStreamIndex);

    // exchange
    updateGrid->exchangeMultiGPUAfterFtoC(level, borderStreamIndex);

    // coarse to fine border
    updateGrid->coarseToFine(level, para->getParD(level)->intCFBorder.ICellCFC, para->getParD(level)->intCFBorder.ICellCFF,
                 para->getParD(level)->intCFBorder.kCF, para->getParD(level)->offCF, borderStreamIndex);
    cudaDeviceSynchronize();
}

void RefinementAndExchange_streams_exchangeAllNodes::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level){
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex = para->getStreamManager()->getBulkStreamIndex();

    // fine to coarse border
    updateGrid->fineToCoarse(level, para->getParD(level)->intFCBorder.ICellFCC, para->getParD(level)->intFCBorder.ICellFCF,
                 para->getParD(level)->intFCBorder.kFC, borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
    updateGrid->prepareExchangeMultiGPU(level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernels (fine to coarse and coarse to fine)
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    updateGrid->fineToCoarse(level, para->getParD(level)->intFCBulk.ICellFCC, para->getParD(level)->intFCBulk.ICellFCF,
                 para->getParD(level)->intFCBulk.kFC, bulkStreamIndex);
    updateGrid->coarseToFine(level, para->getParD(level)->intCFBulk.ICellCFC, para->getParD(level)->intCFBulk.ICellCFF,
                 para->getParD(level)->intCFBulk.kCF, para->getParD(level)->offCFBulk, bulkStreamIndex);

    // exchange
    updateGrid->exchangeMultiGPU(level, borderStreamIndex);

    // coarse to fine border
    updateGrid->coarseToFine(level, para->getParD(level)->intCFBorder.ICellCFC, para->getParD(level)->intCFBorder.ICellCFF,
                 para->getParD(level)->intCFBorder.kCF, para->getParD(level)->offCF, borderStreamIndex);
    cudaDeviceSynchronize();
}

void RefinementAndExchange_noStreams_exchangeInterface::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    updateGrid->fineToCoarse(level, para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, para->getParD(level)->K_FC, -1);

    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, true);

    updateGrid->coarseToFine(level, para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, para->getParD(level)->K_CF,
                             para->getParD(level)->offCF, -1);
}

void RefinementAndExchange_noStreams_exchangeAllNodes::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. Interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, para->getParD(level)->K_FC, -1);
    
    //! 2. Exchange informations between GPUs
    //!
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, false);

    //! 3. Interpolation coarse to fine
    updateGrid->coarseToFine(level, para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, para->getParD(level)->K_CF,
                             para->getParD(level)->offCF, -1);
}

void Refinement_noExchange::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. Interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, para->getParD(level)->K_FC, -1);
    //! 2. Interpolation coarse to fine
    updateGrid->coarseToFine(level, para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, para->getParD(level)->K_CF,
                 para->getParD(level)->offCF, -1);
}
