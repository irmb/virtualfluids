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
//! \file GridScalingKernelManager.h
//! \ingroup KernelManager
//! \author Martin Schoenherr
//=======================================================================================
#include "GridScaling/GridScalingKernelManager.h"
#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include <logger/Logger.h>
#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "PreCollisionInteractor/PreCollisionInteractor.h"
#include "GridScaling/GridScalingFactory.h"
#include <stdexcept>

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

    // ScaleFC_comp_D3Q27F3(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level+1)->distributions.f[0],
    //     para->getParD(level)->g6.g[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellFC->ICellFCC,
    //     icellFC->ICellFCF,
    //     icellFC->kFC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC,
    //     stream);

    // ScaleFC_0817_comp_27(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level+1)->distributions.f[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellFC->ICellFCC,
    //     icellFC->ICellFCF,
    //     icellFC->kFC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC,
    //     stream);

    // ScaleFC_RhoSq_3rdMom_comp_27(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level+1)->distributions.f[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellFC->ICellFCC,
    //     icellFC->ICellFCF,
    //     icellFC->kFC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC,
    //     stream);

    // ScaleFC_AA2016_comp_27(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level+1)->distributions.f[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellFC->ICellFCC,
    //     icellFC->ICellFCF,
    //     icellFC->kFC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC,
    //     stream);


    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////

    //ScaleFC27(  para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //            para->getParD(level)->neighborX,   para->getParD(level)->neighborY,   para->getParD(level)->neighborZ,
    //            para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY, para->getParD(level+1)->neighborZ,
    //            para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,  para->getParD(level)->isEvenTimestep,
    //            para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //            para->getParD(level)->K_FC,           para->getParD(level)->omega,          para->getParD(level+1)->omega,
    //            para->getParD(level)->vis,            para->getParD(level)->nx,             para->getParD(level)->ny,
    //            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,           para->getParD(level)->gridNX);
    //getLastCudaError("ScaleFC27 execution failed");

    //ScaleFCEff27(  para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC27 execution failed");

    //ScaleFCLast27( para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC27 execution failed");

    //ScaleFCpress27(para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC27 execution failed");

    // ScaleFC_Fix_comp_27(    para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //                      para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //                      para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //                      para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //                      para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                      para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                      para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                      para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                      para->getParD(level)->offFC);
    // getLastCudaError("ScaleFC27 execution failed");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // incompressible
    //ScaleFC_Fix_27(para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC27 execution failed");

    //ScaleFC_NSPress_27(para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //                   para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //                   para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //                   para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //                   para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                   para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                   para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                   para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                   para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC27 execution failed");
}

void GridScalingKernelManager::runFineToCoarseKernelAD(const int level) const
{
    //A D V E C T I O N    D I F F U S I O N

        ScaleFCThS27(
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

    // ScaleCF_comp_D3Q27F3(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level+1)->distributions.f[0],
    //     para->getParD(level+1)->g6.g[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellCF->ICellCFC,
    //     icellCF->ICellCFF,
    //     icellCF->kCF,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     offCF,
    //     stream);

    // ScaleCF_0817_comp_27(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level + 1)->distributions.f[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellCF->ICellCFC,
    //     icellCF->ICellCFF,
    //     icellCF->kCF,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     offCF,
    //     stream);

    // ScaleCF_RhoSq_3rdMom_comp_27(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level+1)->distributions.f[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellCF->ICellCFC,
    //     icellCF->ICellCFF,
    //     icellCF->kCF,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     offCF,
    //     stream);

    // ScaleCF_AA2016_comp_27(
    //     para->getParD(level)->distributions.f[0],
    //     para->getParD(level+1)->distributions.f[0],
    //     para->getParD(level)->neighborX,
    //     para->getParD(level)->neighborY,
    //     para->getParD(level)->neighborZ,
    //     para->getParD(level+1)->neighborX,
    //     para->getParD(level+1)->neighborY,
    //     para->getParD(level+1)->neighborZ,
    //     para->getParD(level)->numberOfNodes,
    //     para->getParD(level+1)->numberOfNodes,
    //     para->getParD(level)->isEvenTimestep,
    //     icellCF->ICellCFC,
    //     icellCF->ICellCFF,
    //     icellCF->kCF,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     offCF,
    //     stream);


    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////

    //ScaleCF27(  para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //            para->getParD(level)->neighborX,   para->getParD(level)->neighborY,   para->getParD(level)->neighborZ,
    //            para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY, para->getParD(level+1)->neighborZ,
    //            para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,  para->getParD(level)->isEvenTimestep,
    //            para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //            para->getParD(level)->K_CF,           para->getParD(level)->omega,          para->getParD(level+1)->omega,
    //            para->getParD(level)->vis,            para->getParD(level)->nx,             para->getParD(level)->ny,
    //            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,           para->getParD(level)->gridNX);
    //getLastCudaError("ScaleCF27 execution failed");

    //ScaleCFEff27(  para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF27 execution failed");

    //ScaleCFLast27( para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF27 execution failed");

    //ScaleCFpress27(para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF27 execution failed");

    // ScaleCF_Fix_comp_27(    para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //                      para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //                      para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //                      para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //                      para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                      para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                      para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                      para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                      para->getParD(level)->offCF);
    // getLastCudaError("ScaleCF27 execution failed");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // incompressible
    //ScaleCF_Fix_27(para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //               para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //               para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //               para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //               para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF27 execution failed");

    //ScaleCF_NSPress_27(para->getParD(level)->distributions.f[0],      para->getParD(level+1)->distributions.f[0],
    //                   para->getParD(level)->neighborX,   para->getParD(level)->neighborY,    para->getParD(level)->neighborZ,
    //                   para->getParD(level+1)->neighborX, para->getParD(level+1)->neighborY,  para->getParD(level+1)->neighborZ,
    //                   para->getParD(level)->numberOfNodes,    para->getParD(level+1)->numberOfNodes,   para->getParD(level)->isEvenTimestep,
    //                   para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                   para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                   para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                   para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                   para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF27 execution failed");
}

void GridScalingKernelManager::runCoarseToFineKernelAD(const int level) const
{
    // A D V E C T I O N    D I F F U S I O N

        ScaleCFThS27(
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
