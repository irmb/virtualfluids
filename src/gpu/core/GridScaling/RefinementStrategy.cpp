//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_GridScaling GridScaling
//! \ingroup gpu_core core
//! \{
//! \author Anna Wellmann, Martin Schönherr
//=======================================================================================
#include "RefinementStrategy.h"

#include <logger/Logger.h>

#include "Calculation/UpdateGrid27.h"
#include "Cuda/CudaStreamManager.h"
#include "Parameter/Parameter.h"

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
    updateGrid->fineToCoarse(level, &para->getParD(level)->fineToCoarseBorder, para->getParD(level)->neighborFineToCoarse, CudaStreamIndex::SubDomainBorder);

    //! 2. prepare the exchange between gpus (collect the send nodes for communication in a buffer on the gpu) and trigger bulk kernel execution when finished
    //!
    updateGrid->prepareExchangeMultiGPUAfterFtoC(level, CudaStreamIndex::SubDomainBorder);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(CudaStreamIndex::SubDomainBorder);

    //! 3. launch the bulk kernels for both interpolation processes (fine to coarse and coarse to fine)
    //!
    para->getStreamManager()->waitOnStartBulkKernelEvent(CudaStreamIndex::Bulk);
    updateGrid->fineToCoarse(level, &para->getParD(level)->fineToCoarseBulk, para->getParD(level)->neighborFineToCoarseBulk, CudaStreamIndex::SubDomainBorder);
    updateGrid->coarseToFine(level, &para->getParD(level)->coarseToFineBulk, para->getParD(level)->neighborCoarseToFineBulk, CudaStreamIndex::SubDomainBorder);

    //! 4. exchange information between GPUs (only nodes which are part of the interpolation)
    //!
    updateGrid->exchangeMultiGPUAfterFtoC(level, CudaStreamIndex::SubDomainBorder);

    // 5. interpolation fine to coarse for nodes which are at the border of the gpus/processes
    //!
    updateGrid->coarseToFine(level, &para->getParD(level)->coarseToFineBorder, para->getParD(level)->neighborCoarseToFine, CudaStreamIndex::SubDomainBorder);

    cudaDeviceSynchronize();
}

void RefinementAndExchange_streams_exchangeAllNodes::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse for nodes which are at the border of the gpus/processes
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->fineToCoarseBorder, para->getParD(level)->neighborFineToCoarse, CudaStreamIndex::SubDomainBorder);

    //! 2. prepare the exchange between gpus (collect the send nodes for communication in a buffer on the gpu) and trigger bulk kernel execution when finished
    //!
    updateGrid->prepareExchangeMultiGPU(level, CudaStreamIndex::SubDomainBorder);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(CudaStreamIndex::SubDomainBorder);

    //! 3. launch the bulk kernels for both interpolation processes (fine to coarse and coarse to fine)
    //!
    para->getStreamManager()->waitOnStartBulkKernelEvent(CudaStreamIndex::Bulk);
    updateGrid->fineToCoarse(level, &para->getParD(level)->fineToCoarseBulk, para->getParD(level)->neighborFineToCoarseBulk, CudaStreamIndex::SubDomainBorder);
    updateGrid->coarseToFine(level, &para->getParD(level)->coarseToFineBulk, para->getParD(level)->neighborCoarseToFineBulk, CudaStreamIndex::SubDomainBorder);

    //! 4. exchange information between GPUs (all nodes)
    //!
    updateGrid->exchangeMultiGPU(level, CudaStreamIndex::SubDomainBorder);

    // 5. interpolation fine to coarse for nodes which are at the border of the gpus/processes
    //!
    updateGrid->coarseToFine(level, &para->getParD(level)->coarseToFineBorder, para->getParD(level)->neighborCoarseToFine, CudaStreamIndex::SubDomainBorder);

    cudaDeviceSynchronize();
}

void RefinementAndExchange_noStreams_exchangeInterface::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->fineToCoarse, para->getParD(level)->neighborFineToCoarse, CudaStreamIndex::Legacy);

    //! 2. exchange information between GPUs (only nodes which are part of the interpolation)
    //!
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, true);

    //! 3. interpolation coarse to fine
    updateGrid->coarseToFine(level, &para->getParD(level)->coarseToFine, para->getParD(level)->neighborCoarseToFine, CudaStreamIndex::Legacy);
}

void RefinementAndExchange_noStreams_exchangeAllNodes::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->fineToCoarse, para->getParD(level)->neighborFineToCoarse, CudaStreamIndex::Legacy);

    //! 2. exchange information between GPUs (all nodes)
    //!
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, false);

    //! 3. interpolation coarse to fine
    updateGrid->coarseToFine(level, &para->getParD(level)->coarseToFine, para->getParD(level)->neighborCoarseToFine, CudaStreamIndex::Legacy);
}

void Refinement_noExchange::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level)
{
    //! \details steps:
    //!
    //! 1. interpolation fine to coarse
    //!
    updateGrid->fineToCoarse(level, &para->getParD(level)->fineToCoarse, para->getParD(level)->neighborFineToCoarse, CudaStreamIndex::Legacy);
    //! 2. interpolation coarse to fine
    updateGrid->coarseToFine(level, &para->getParD(level)->coarseToFine, para->getParD(level)->neighborCoarseToFine, CudaStreamIndex::Legacy);
}

//! \}
