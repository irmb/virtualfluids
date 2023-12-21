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
//! \addtogroup gpu_Calculation Calculation
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Stephan Lenz, Anna Wellmann
//=======================================================================================
#include "UpdateGrid27.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <logger/Logger.h>

#include <parallel/Communicator.h>

#include "BoundaryConditions/BoundaryConditionKernelManager.h"
#include "CollisionStrategy.h"
#include "Communication/ExchangeData27.h"
#include "Cuda/CudaStreamManager.h"
#include "GridScaling/GridScalingKernelManager.h"
#include "GridScaling/RefinementStrategy.h"
#include "Kernel/ADKernelManager.h"
#include "Kernel/Kernel.h"
#include "PostProcessor/MacroscopicQuantities.cuh"
#include "TurbulenceModels/TurbulenceModelFactory.h"

void UpdateGrid27::updateGrid(int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////

    if (level != para->getFine()) {
        updateGrid(level + 1, t);
        updateGrid(level + 1, t);
    }

    //////////////////////////////////////////////////////////////////////////
    
    interactWithProbes(level, t);

    //////////////////////////////////////////////////////////////////////////

    collision(this, para.get(), level, t);

    //////////////////////////////////////////////////////////////////////////

    postCollisionBC(level, t);

    //////////////////////////////////////////////////////////////////////////

    swapBetweenEvenAndOddTimestep(level);

    //////////////////////////////////////////////////////////////////////////

    calcTurbulentViscosity(level);

    //////////////////////////////////////////////////////////////////////////

    this->preCollisionBC(level, t);

    //////////////////////////////////////////////////////////////////////////
    if( level != para->getFine() )
    {   
        refinement(this, para.get(), level);
    }

    //////////////////////////////////////////////////////////////////////////
    
    interactWithActuators(level, t);

}

void UpdateGrid27::collisionAllNodes(int level, unsigned int t)
{
    kernels.at(level)->run();

    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
        collisionAdvectionDiffusion(level);
}

void UpdateGrid27::collisionUsingIndices(int level, unsigned int t, uint *taggedFluidNodeIndices, uint numberOfTaggedFluidNodes, CollisionTemplate collisionTemplate, CudaStreamIndex stream)
{
    if (taggedFluidNodeIndices != nullptr && numberOfTaggedFluidNodes != 0)
        kernels.at(level)->runOnIndices(taggedFluidNodeIndices, numberOfTaggedFluidNodes, collisionTemplate, stream);
    else
        VF_LOG_CRITICAL("In collision: fluidNodeIndices or numberOfFluidNodes not defined (level = {})", level);

    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
        collisionAdvectionDiffusion(level);
}

void UpdateGrid27::collisionAdvectionDiffusion(int level)
{
    this->adKernelManager->runADcollisionKernel(level);
}

void UpdateGrid27::prepareExchangeMultiGPU(int level, CudaStreamIndex streamIndex)
{
    prepareExchangeCollDataXGPU27AllNodes(para.get(), level, streamIndex);
    prepareExchangeCollDataYGPU27AllNodes(para.get(), level, streamIndex);
    prepareExchangeCollDataZGPU27AllNodes(para.get(), level, streamIndex);
}

void UpdateGrid27::prepareExchangeMultiGPUAfterFtoC(int level, CudaStreamIndex streamIndex)
{
    prepareExchangeCollDataXGPU27AfterFtoC(para.get(), level, streamIndex);
    prepareExchangeCollDataYGPU27AfterFtoC(para.get(), level, streamIndex);
    prepareExchangeCollDataZGPU27AfterFtoC(para.get(), level, streamIndex);
}

void UpdateGrid27::exchangeMultiGPU(int level, CudaStreamIndex streamIndex)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataYGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataZGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);

    scatterNodesFromRecvBufferXGPU27AllNodes(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferYGPU27AllNodes(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferZGPU27AllNodes(para.get(), level, streamIndex);

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            VF_LOG_WARNING("Warning: Cuda streams not yet implemented for convection diffusion");
        exchangePostCollDataADXGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADYGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADZGPU27(para.get(), comm, cudaMemoryManager.get(), level);
    }

}
void UpdateGrid27::exchangeMultiGPU_noStreams_withPrepare(int level, bool useReducedComm)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    if (useReducedComm) {
        // X
        prepareExchangeCollDataXGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataXGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferXGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        // Y
        prepareExchangeCollDataYGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataYGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferYGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        // Z
        prepareExchangeCollDataZGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataZGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferZGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
    } else {
        // X
        prepareExchangeCollDataXGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataXGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferXGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        // Y
        prepareExchangeCollDataYGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataYGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferYGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        // Z
        prepareExchangeCollDataZGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataZGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferZGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
    }

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            VF_LOG_WARNING("Warning: Cuda streams not yet implemented for convection diffusion");
        exchangePostCollDataADXGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADYGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADZGPU27(para.get(), comm, cudaMemoryManager.get(), level);
    }
}
void UpdateGrid27::exchangeMultiGPUAfterFtoC(int level, CudaStreamIndex streamIndex)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataYGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataZGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);

    scatterNodesFromRecvBufferXGPU27AfterFtoC(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferYGPU27AfterFtoC(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferZGPU27AfterFtoC(para.get(), level, streamIndex);

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            VF_LOG_WARNING("Warning: Cuda streams not yet implemented for convection diffusion");
        exchangePostCollDataADXGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADYGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADZGPU27(para.get(), comm, cudaMemoryManager.get(), level);
    }
}

void UpdateGrid27::postCollisionBC(int level, uint t)
{
    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    // V E L O C I T Y (I N F L O W)
    this->bcKernelManager->runVelocityBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // N O - S L I P
    this->bcKernelManager->runNoSlipBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // S L I P
    this->bcKernelManager->runSlipBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // S T R E S S (wall model)
    this->bcKernelManager->runStressWallModelKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    this->bcKernelManager->runGeoBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    this->bcKernelManager->runOutflowBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////
    // P R E C U R S O R
    this->bcKernelManager->runPrecursorBCKernelPost(level, t, cudaMemoryManager.get());

    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    if (para->getDiffOn())
    {
        this->adKernelManager->runADgeometryBCKernel(level);
        this->adKernelManager->runADDirichletBCKernel(level);
        this->adKernelManager->runADslipBCKernel(level);
    }
}

void UpdateGrid27::swapBetweenEvenAndOddTimestep(int level)
{
    if (para->getParD(level)->isEvenTimestep==true)  para->getParD(level)->isEvenTimestep=false;
    else                                        para->getParD(level)->isEvenTimestep=true;
}

void UpdateGrid27::calcMacroscopicQuantities(int level)
{
    calculateMacroscopicQuantitiesCompressible(para->getParD(level)->velocityX,
                    para->getParD(level)->velocityY,
                    para->getParD(level)->velocityZ,
                    para->getParD(level)->rho,
                    para->getParD(level)->pressure,
                    para->getParD(level)->typeOfGridNode,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->numberofthreads,
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->isEvenTimestep);
    getLastCudaError("calculateMacroscopicQuantities execution failed");
}

void UpdateGrid27::preCollisionBC(int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////
    // V E L O C I T Y (I N F L O W)
    this->bcKernelManager->runVelocityBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    this->bcKernelManager->runGeoBCKernelPre(level, t, cudaMemoryManager.get());

    //////////////////////////////////////////////////////////////////////////
    // P R E S S U R E
    this->bcKernelManager->runPressureBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    this->bcKernelManager->runOutflowBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////////////
    ////only for a round off error test
    //para->cudaCopyTestREtoHost(0,para->getParH(0)->pressureBC.numberOfBCnodes);
    //printRE(para, t);
    //////////////////////////////////////////////////////////////////////////////////
}

void UpdateGrid27::fineToCoarse(int level, InterpolationCells* fineToCoarse, ICellNeigh &neighborFineToCoarse, CudaStreamIndex streamIndex)
{
    gridScalingKernelManager->runFineToCoarseKernelLB(level, fineToCoarse, neighborFineToCoarse, streamIndex);
}

void UpdateGrid27::coarseToFine(int level, InterpolationCells* coarseToFine, ICellNeigh &neighborCoarseToFine, CudaStreamIndex streamIndex)
{
    this->gridScalingKernelManager->runCoarseToFineKernelLB(level, coarseToFine, neighborCoarseToFine, streamIndex);
}

void UpdateGrid27::interactWithActuators(int level, unsigned int t)
{
    for( SPtr<PreCollisionInteractor> actuator: para->getActuators() )
    {
        actuator->interact(level, t);
    }
}

void  UpdateGrid27::interactWithProbes(int level, unsigned int t)
{
    for( SPtr<PreCollisionInteractor> probe: para->getProbes() )
    {
        probe->interact(level, t);
    }
}

void  UpdateGrid27::calcTurbulentViscosity(int level)
{
    this->tmFactory->runTurbulenceModelKernel(level);
}

void UpdateGrid27::exchangeData(int level)
{
    exchangeMultiGPU_noStreams_withPrepare(level, false);
}

UpdateGrid27::UpdateGrid27(SPtr<Parameter> para, vf::parallel::Communicator &comm, SPtr<CudaMemoryManager> cudaMemoryManager,
                           std::vector<SPtr<Kernel>>& kernels,
                           std::vector<SPtr<AdvectionDiffusionKernel>>& adkernels, BoundaryConditionFactory* bcFactory,
                           SPtr<TurbulenceModelFactory> tmFactory, GridScalingFactory* scalingFactory)
    : para(para), comm(comm), cudaMemoryManager(cudaMemoryManager), kernels(kernels), tmFactory(tmFactory)
{
    this->collision = getFunctionForCollisionAndExchange(para->getUseStreams(), para->getNumprocs(), para->getKernelNeedsFluidNodeIndicesToRun());
    this->refinement = getFunctionForRefinementAndExchange(para->getUseStreams(), para->getNumprocs(), para->getMaxLevel(), para->useReducedCommunicationAfterFtoC);

    this->bcKernelManager = std::make_shared<BoundaryConditionKernelManager>(para, bcFactory);
    this->adKernelManager = std::make_shared<ADKernelManager>(para, adkernels);
    this->gridScalingKernelManager = std::make_shared<GridScalingKernelManager>(para, scalingFactory);
}

//! \}
