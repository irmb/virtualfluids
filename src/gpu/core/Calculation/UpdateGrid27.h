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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef UPDATEGRID27_H
#define UPDATEGRID27_H

#include <vector>

#include "Cuda/CudaMemoryManager.h"

#include "Calculation/Calculation.h"
#include "Cuda/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include "Kernel/AdvectionDiffusionKernel.h"

namespace vf::parallel
{
class Communicator;
}

class BoundaryConditionKernelManager;
class ADKernelManager;
class GridScalingKernelManager;
class Kernel;
class BoundaryConditionFactory;
class GridScalingFactory;
class TurbulenceModelFactory;
class UpdateGrid27;
using CollisionStrategy = std::function<void (UpdateGrid27* updateGrid, Parameter* para, int level, unsigned int t)>;
using RefinementStrategy = std::function<void (UpdateGrid27* updateGrid, Parameter* para, int level)>;


class UpdateGrid27
{
public:
    UpdateGrid27(SPtr<Parameter> para, vf::parallel::Communicator& comm, SPtr<CudaMemoryManager> cudaMemoryManager, std::vector<SPtr<Kernel>>& kernels,
                 std::vector<SPtr<AdvectionDiffusionKernel>>& adkernels, BoundaryConditionFactory* bcFactory,
                 SPtr<TurbulenceModelFactory> tmFactory, GridScalingFactory* scalingFactory);
    void updateGrid(int level, unsigned int t);
    void exchangeData(int level);

private:
    void collisionAllNodes(int level, unsigned int t);
    void collisionUsingIndices(int level, unsigned int t, uint *taggedFluidNodeIndices = nullptr, uint numberOfTaggedFluidNodes = 0, CollisionTemplate collisionTemplate = CollisionTemplate::Default, CudaStreamIndex streamIndex=CudaStreamIndex::Legacy);
    void collisionAdvectionDiffusion(int level);

    void postCollisionBC(int level, unsigned int t);
    void preCollisionBC(int level, unsigned int t);

    void fineToCoarse(int level, InterpolationCells* fineToCoarse, ICellNeigh &neighborFineToCoarse, CudaStreamIndex streamIndex);
    void coarseToFine(int level, InterpolationCells* coarseToFine, ICellNeigh &neighborCoarseToFine, CudaStreamIndex streamIndex);

    void prepareExchangeMultiGPU(int level, CudaStreamIndex streamIndex);
    void prepareExchangeMultiGPUAfterFtoC(int level, CudaStreamIndex streamIndex);

    void exchangeMultiGPU(int level, CudaStreamIndex streamIndex);
    void exchangeMultiGPUAfterFtoC(int level, CudaStreamIndex streamIndex);
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
    friend class CollisionAndExchange_noStreams_withReadWriteFlags;

    RefinementStrategy refinement;
    friend class RefinementAndExchange_streams_exchangeInterface;
    friend class RefinementAndExchange_streams_exchangeAllNodes;
    friend class RefinementAndExchange_noStreams_exchangeInterface;
    friend class RefinementAndExchange_noStreams_exchangeAllNodes;
    friend class Refinement_noExchange;
    friend class NoRefinement;

private:
    SPtr<Parameter> para;
    vf::parallel::Communicator& comm;
    SPtr<CudaMemoryManager> cudaMemoryManager;
    std::vector<SPtr<Kernel>> kernels;
    //! \property lbKernelManager is a shared pointer to an object of LBKernelManager
    std::shared_ptr<BoundaryConditionKernelManager> bcKernelManager;
    //! \property adKernelManager is a shared pointer to an object of ADKernelManager
    std::shared_ptr<ADKernelManager> adKernelManager;
    //! \property gridScalingKernelManager is a shared pointer to an object of GridScalingKernelManager
    std::shared_ptr<GridScalingKernelManager> gridScalingKernelManager;
    //! \property tmFactory is a shared pointer to an object of TurbulenceModelFactory
    std::shared_ptr<TurbulenceModelFactory> tmFactory;
};

#endif

//! \}
