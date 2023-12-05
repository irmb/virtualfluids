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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr
//=======================================================================================
#include "CollisionStrategy.h"
#include "Cuda/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include <logger/Logger.h>

std::function<void(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t)>
getFunctionForCollisionAndExchange(const bool useStreams, const int numberOfMpiProcesses,
                                   const bool kernelNeedsFluidNodeIndicesToRun)
{
    VF_LOG_INFO("Function used for collisionAndExchange: ");

    if (useStreams && numberOfMpiProcesses > 1 && kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("CollisionAndExchange_streams()");
        return CollisionAndExchange_streams();

    } else if (useStreams && !kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("Cuda Streams can only be used with kernels which run using fluidNodesIndices.");

    } else if (useStreams && numberOfMpiProcesses <= 1) {
        VF_LOG_INFO("Cuda Streams can only be used with multiple MPI processes.");

    } else if (!useStreams && kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("CollisionAndExchange_noStreams_indexKernel()");
        return CollisionAndExchange_noStreams_indexKernel();

    } else if (!useStreams && !kernelNeedsFluidNodeIndicesToRun) {
        VF_LOG_INFO("CollisionAndExchange_noStreams_oldKernel()");
        return CollisionAndExchange_noStreams_oldKernel();
    }

    throw std::runtime_error("Invalid Configuration for collision and exchange");
}

void CollisionAndExchange_noStreams_indexKernel::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level,
                                                            unsigned int t)
{
    //! \details steps:
    //!
    //! 1. run collision
    //!
    for( CollisionTemplate tag: para->getParH(level)->allocatedBulkFluidNodeTags )
    {
        updateGrid->collisionUsingIndices(  level, t, 
                                            para->getParD(level)->taggedFluidNodeIndices[tag],
                                            para->getParD(level)->numberOfTaggedFluidNodes[tag],
                                            tag,
                                            CudaStreamIndex::Legacy);
    }

    //! 2. exchange information between GPUs
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, false);
}

void CollisionAndExchange_noStreams_oldKernel::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level,
                                                          unsigned int t)
{
    //! \details steps:
    //!
    //! 1. run collision
    //!
    updateGrid->collisionAllNodes(level, t);

    //! 2. exchange information between GPUs
    updateGrid->exchangeMultiGPU_noStreams_withPrepare(level, false);
}

void CollisionAndExchange_streams::operator()(UpdateGrid27 *updateGrid, Parameter *para, int level, unsigned int t)
{
    //! \details steps:
    //!
    //! 1. run collision for nodes which are at the border of the gpus/processes, running with WriteMacroVars in case probes sample on these nodes
    //!    
    updateGrid->collisionUsingIndices(  level, t, 
                                        para->getParD(level)->taggedFluidNodeIndices[CollisionTemplate::SubDomainBorder],
                                        para->getParD(level)->numberOfTaggedFluidNodes[CollisionTemplate::SubDomainBorder], 
                                        CollisionTemplate::WriteMacroVars,  
                                        CudaStreamIndex::SubDomainBorder);

    //! 2. prepare the exchange between gpus (collect the send nodes for communication in a buffer on the gpu) and trigger bulk kernel execution when finished
    //!
    updateGrid->prepareExchangeMultiGPU(level, CudaStreamIndex::SubDomainBorder);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(CudaStreamIndex::SubDomainBorder);

    //! 3. launch the collision kernel for bulk nodes. This includes nodes with \param tag Default, WriteMacroVars, ApplyBodyForce, 
    //!    or AllFeatures. All assigned tags are listed in \param allocatedBulkFluidNodeTags during initialization in Simulation::init

    para->getStreamManager()->waitOnStartBulkKernelEvent(CudaStreamIndex::Bulk);
    
    for( CollisionTemplate tag: para->getParH(level)->allocatedBulkFluidNodeTags )
    {
        updateGrid->collisionUsingIndices(  level, t, 
                                            para->getParD(level)->taggedFluidNodeIndices[tag],
                                            para->getParD(level)->numberOfTaggedFluidNodes[tag], 
                                            tag,
                                            CudaStreamIndex::Bulk);
    }
    //! 4. exchange information between GPUs
    updateGrid->exchangeMultiGPU(level, CudaStreamIndex::SubDomainBorder);
}
