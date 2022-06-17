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
#include "KernelManager/GridScalingKernelManager.h"
#include "GPU/CudaMemoryManager.h"
#include "GPU/GPU_Interface.h"
#include "Parameter/Parameter.h"

SPtr<GridScalingKernelManager> GridScalingKernelManager::make(SPtr<Parameter> parameter){
    return SPtr<GridScalingKernelManager>(new GridScalingKernelManager(parameter));
}

GridScalingKernelManager::GridScalingKernelManager(SPtr<Parameter> parameter): para(parameter){}

GridScalingKernelManager::GridScalingKernelManager(const GridScalingKernelManager&){}

void GridScalingKernelManager::runFineToCoarseKernelLB(int level){
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
    //     para->getParD(level)->intFC.ICellFCC,
    //     para->getParD(level)->intFC.ICellFCF,
    //     para->getParD(level)->K_FC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC);

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
    //     para->getParD(level)->intFC.ICellFCC,
    //     para->getParD(level)->intFC.ICellFCF,
    //     para->getParD(level)->K_FC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC);

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
    //     para->getParD(level)->intFC.ICellFCC,
    //     para->getParD(level)->intFC.ICellFCF,
    //     para->getParD(level)->K_FC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC);

    ScaleFC_RhoSq_comp_27(
        para->getParD(level)->distributions.f[0],
        para->getParD(level+1)->distributions.f[0],
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level+1)->neighborX,
        para->getParD(level+1)->neighborY,
        para->getParD(level+1)->neighborZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level+1)->numberOfNodes,
        para->getParD(level)->isEvenTimestep,
        para->getParD(level)->intFC.ICellFCC,
        para->getParD(level)->intFC.ICellFCF,
        para->getParD(level)->K_FC,
        para->getParD(level)->omega,
        para->getParD(level+1)->omega,
        para->getParD(level)->vis,
        para->getParD(level)->nx,
        para->getParD(level)->ny,
        para->getParD(level+1)->nx,
        para->getParD(level+1)->ny,
        para->getParD(level)->numberofthreads,
        para->getParD(level)->offFC,
        CU_STREAM_LEGACY);

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
    //     para->getParD(level)->intFC.ICellFCC,
    //     para->getParD(level)->intFC.ICellFCF,
    //     para->getParD(level)->K_FC,
    //     para->getParD(level)->omega,
    //     para->getParD(level+1)->omega,
    //     para->getParD(level)->vis,
    //     para->getParD(level)->nx,
    //     para->getParD(level)->ny,
    //     para->getParD(level+1)->nx,
    //     para->getParD(level+1)->ny,
    //     para->getParD(level)->numberofthreads,
    //     para->getParD(level)->offFC);



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

void GridScalingKernelManager::runFineToCoarseKernelAD(int level){
    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    if (para->getDiffOn())
    {
        if (para->getDiffMod() == 7)
        {
            // ScaleFCThS7(
            //     para->getParD(level)->distributions.f[0],
            //     para->getParD(level+1)->distributions.f[0],
            //     para->getParD(level)->distributionsAD7.f[0],
            //     para->getParD(level+1)->distributionsAD7.f[0],
            //     para->getParD(level)->neighborX,
            //     para->getParD(level)->neighborY,
            //     para->getParD(level)->neighborZ,
            //     para->getParD(level+1)->neighborX,
            //     para->getParD(level+1)->neighborY,
            //     para->getParD(level+1)->neighborZ,
            //     para->getParD(level)->numberOfNodes,
            //     para->getParD(level+1)->numberOfNodes,
            //     para->getParD(level)->isEvenTimestep,
            //     para->getParD(level)->intFC.ICellFCC,
            //     para->getParD(level)->intFC.ICellFCF,
            //     para->getParD(level)->K_FC,
            //     para->getParD(level)->vis,
            //     para->getParD(level)->diffusivity,
            //     para->getParD(level)->numberofthreads);

            ScaleFCThSMG7(
                para->getParD(level)->distributions.f[0],
                para->getParD(level+1)->distributions.f[0],
                para->getParD(level)->distributionsAD7.f[0],
                para->getParD(level+1)->distributionsAD7.f[0],
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level+1)->neighborX,
                para->getParD(level+1)->neighborY,
                para->getParD(level+1)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level+1)->numberOfNodes,
                para->getParD(level)->isEvenTimestep,
                para->getParD(level)->intFC.ICellFCC,
                para->getParD(level)->intFC.ICellFCF,
                para->getParD(level)->K_FC,
                para->getParD(level)->vis,
                para->getParD(level)->diffusivity,
                para->getParD(level)->numberofthreads,
                para->getParD(level)->offFC);

        }
        else if (para->getDiffMod() == 27)
        {
            ScaleFCThS27(
                para->getParD(level)->distributions.f[0],
                para->getParD(level+1)->distributions.f[0],
                para->getParD(level)->distributionsAD27.f[0],
                para->getParD(level+1)->distributionsAD27.f[0],
                para->getParD(level)->neighborX,
                para->getParD(level)->neighborY,
                para->getParD(level)->neighborZ,
                para->getParD(level+1)->neighborX,
                para->getParD(level+1)->neighborY,
                para->getParD(level+1)->neighborZ,
                para->getParD(level)->numberOfNodes,
                para->getParD(level+1)->numberOfNodes,
                para->getParD(level)->isEvenTimestep,
                para->getParD(level)->intFC.ICellFCC,
                para->getParD(level)->intFC.ICellFCF,
                para->getParD(level)->K_FC,
                para->getParD(level)->vis,
                para->getParD(level)->diffusivity,
                para->getParD(level)->numberofthreads,
                para->getParD(level)->offFC);
        }
    }
}