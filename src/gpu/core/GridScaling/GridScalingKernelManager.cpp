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
#include "GridScaling/GridScalingKernelManager.h"
#include "Cuda/CudaMemoryManager.h"

#include <logger/Logger.h>
#include "Parameter/Parameter.h"
#include "Cuda/CudaStreamManager.h"
#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include "GridScaling/GridScalingFactory.h"
#include <stdexcept>
#include "GridScaling/Scaling.cuh"

GridScalingKernelManager::GridScalingKernelManager(SPtr<Parameter> parameter, GridScalingFactory *gridScalingFactory)
    : para(parameter)
{
    if(para->getMaxLevel() != 0){
        if(!gridScalingFactory){
            throw std::runtime_error("There is more than one level, but no scalingFactory was provided.");
        }
        checkScalingFunction(gridScalingFactory->getGridScalingFC(parameter->getUseTurbulentViscosity()), this->para->getParD(0)->fineToCoarse, "scalingFineToCoarse");
        checkScalingFunction(gridScalingFactory->getGridScalingCF(parameter->getUseTurbulentViscosity()), this->para->getParD(0)->coarseToFine, "scalingCoarseToFine");
        this->scalingFineToCoarse = gridScalingFactory->getGridScalingFC(parameter->getUseTurbulentViscosity());
        this->scalingCoarseToFine = gridScalingFactory->getGridScalingCF(parameter->getUseTurbulentViscosity());
    }
    
    if(this->scalingFineToCoarse == nullptr)
        VF_LOG_TRACE("Function for scalingFineToCoarse is nullptr");
    if(this->scalingCoarseToFine == nullptr)
        VF_LOG_TRACE("Function for scalingCoarseToFine is nullptr");
}

void GridScalingKernelManager::runFineToCoarseKernelLB(const int level, InterpolationCells *fineToCoarse, ICellNeigh &neighborFineToCoarse, CudaStreamIndex streamIndex) const
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    this->scalingFineToCoarse(para->getParD(level).get(), para->getParD(level+1).get(), fineToCoarse, neighborFineToCoarse, stream);
}

void GridScalingKernelManager::runFineToCoarseKernelAD(const int level) const
{
    scaleFineToCoarseAdvectionDiffusion(
        para->getParD(level)->distributions.f[0],
        para->getParD(level+1)->distributions.f[0],
        para->getParD(level)->distributionsAD.f[0],
        para->getParD(level+1)->distributionsAD.f[0],
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level+1)->neighborX,
        para->getParD(level+1)->neighborY,
        para->getParD(level+1)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level+1)->numberOfNodes,
        para->getParD(level)->isEvenTimestep,
        para->getParD(level)->fineToCoarse.coarseCellIndices,
        para->getParD(level)->fineToCoarse.fineCellIndices,
        para->getParD(level)->fineToCoarse.numberOfCells,
        para->getParD(level)->viscosity,
        para->getParD(level)->diffusivity,
        para->getParD(level)->numberofthreads,
        para->getParD(level)->neighborFineToCoarse);
}

void GridScalingKernelManager::runCoarseToFineKernelLB(const int level, InterpolationCells* coarseToFine, ICellNeigh &neighborFineToCoarse, CudaStreamIndex streamIndex) const
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);
    this->scalingCoarseToFine(para->getParD(level).get(), para->getParD(level+1).get(), coarseToFine, neighborFineToCoarse, stream);
}

void GridScalingKernelManager::runCoarseToFineKernelAD(const int level) const
{
    scaleCoarseToFineAdvectionDiffusion(
        para->getParD(level)->distributions.f[0],
        para->getParD(level+1)->distributions.f[0],
        para->getParD(level)->distributionsAD.f[0],
        para->getParD(level+1)->distributionsAD.f[0],
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level+1)->neighborX,
        para->getParD(level+1)->neighborY,
        para->getParD(level+1)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level+1)->numberOfNodes,
        para->getParD(level)->isEvenTimestep,
        para->getParD(level)->coarseToFine.coarseCellIndices,
        para->getParD(level)->coarseToFine.fineCellIndices,
        para->getParD(level)->coarseToFine.numberOfCells,
        para->getParD(level)->viscosity,
        para->getParD(level+1)->diffusivity,
        para->getParD(level)->numberofthreads,
        para->getParD(level)->neighborCoarseToFine);
}
