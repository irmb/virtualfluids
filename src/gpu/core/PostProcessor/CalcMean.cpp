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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \
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
#include "CalcMean.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Cuda/CudaMemoryManager.h"
#include "Parameter/Parameter.h"
#include "PostProcessor/MacroscopicQuantities.cuh"

void allocMean(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaMemoryManager->cudaAllocMeanOut(lev);
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++) {
            para->getParH(lev)->vx_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->vy_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->vz_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->rho_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->press_SP_Med_Out[pos] = (real)0.0;
        }
    }
}

void calcMean(Parameter* para, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++) {
            para->getParH(lev)->vx_SP_Med_Out[pos] = para->getParH(lev)->vx_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->vy_SP_Med_Out[pos] = para->getParH(lev)->vy_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->vz_SP_Med_Out[pos] = para->getParH(lev)->vz_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->rho_SP_Med_Out[pos] = para->getParH(lev)->rho_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->press_SP_Med_Out[pos] = para->getParH(lev)->press_SP_Med[pos] / (real)tdiff;
        }
    }
}

void resetMean(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        ResetMeanValuesSP27(para->getParD(lev)->vx_SP_Med, para->getParD(lev)->vy_SP_Med, para->getParD(lev)->vz_SP_Med,
                            para->getParD(lev)->rho_SP_Med, para->getParD(lev)->press_SP_Med,
                            para->getParD(lev)->numberOfNodes, para->getParD(lev)->numberofthreads,
                            para->getParD(lev)->isEvenTimestep);
        getLastCudaError("ResetMeanValuesSP27 execution failed");
    }
}

// Advection-Diffusion
void allocMeanAD(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        cudaMemoryManager->cudaAllocMeanOutAD(lev);
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++) {
            para->getParH(lev)->vx_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->vy_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->vz_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->rho_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->press_SP_Med_Out[pos] = (real)0.0;
            para->getParH(lev)->Conc_Med_Out[pos] = (real)0.0;
        }
    }
}

void calcMeanAD(Parameter* para, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++) {
            para->getParH(lev)->vx_SP_Med_Out[pos] = para->getParH(lev)->vx_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->vy_SP_Med_Out[pos] = para->getParH(lev)->vy_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->vz_SP_Med_Out[pos] = para->getParH(lev)->vz_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->rho_SP_Med_Out[pos] = para->getParH(lev)->rho_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->press_SP_Med_Out[pos] = para->getParH(lev)->press_SP_Med[pos] / (real)tdiff;
            para->getParH(lev)->Conc_Med_Out[pos] = para->getParH(lev)->Conc_Med[pos] / (real)tdiff;
        }
    }
}

void resetMeanAD(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        ResetMeanValuesAD27(para->getParD(lev)->vx_SP_Med, para->getParD(lev)->vy_SP_Med, para->getParD(lev)->vz_SP_Med,
                            para->getParD(lev)->rho_SP_Med, para->getParD(lev)->press_SP_Med, para->getParD(lev)->Conc_Med,
                            para->getParD(lev)->numberOfNodes, para->getParD(lev)->numberofthreads,
                            para->getParD(lev)->isEvenTimestep);
        getLastCudaError("ResetMeanValuesAD27 execution failed");
    }
}
