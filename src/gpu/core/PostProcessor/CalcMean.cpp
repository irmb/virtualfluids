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
//! \addtogroup gpu_PostProcessor PostProcessor
//! \ingroup gpu_core core
//! \{
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
            para->getParH(lev)->meanVelocityInXdirectionOut[pos] = (real)0.0;
            para->getParH(lev)->meanVelocityInYdirectionOut[pos] = (real)0.0;
            para->getParH(lev)->meanVelocityInZdirectionOut[pos] = (real)0.0;
            para->getParH(lev)->meanDensityOut[pos] = (real)0.0;
            para->getParH(lev)->meanPressureOut[pos] = (real)0.0;
        }
    }
}

void calcMean(Parameter* para, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++) {
            para->getParH(lev)->meanVelocityInXdirectionOut[pos] = para->getParH(lev)->meanVelocityInXdirection[pos] / (real)tdiff;
            para->getParH(lev)->meanVelocityInYdirectionOut[pos] = para->getParH(lev)->meanVelocityInYdirection[pos] / (real)tdiff;
            para->getParH(lev)->meanVelocityInZdirectionOut[pos] = para->getParH(lev)->meanVelocityInZdirection[pos] / (real)tdiff;
            para->getParH(lev)->meanDensityOut[pos] = para->getParH(lev)->meanDensity[pos] / (real)tdiff;
            para->getParH(lev)->meanPressureOut[pos] = para->getParH(lev)->meanPressure[pos] / (real)tdiff;
        }
    }
}

void resetMean(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        ResetMeanValuesSP27(para->getParD(lev)->meanVelocityInXdirection, para->getParD(lev)->meanVelocityInYdirection, para->getParD(lev)->meanVelocityInZdirection,
                            para->getParD(lev)->meanDensity, para->getParD(lev)->meanPressure,
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
            para->getParH(lev)->meanVelocityInXdirectionOut[pos] = (real)0.0;
            para->getParH(lev)->meanVelocityInYdirectionOut[pos] = (real)0.0;
            para->getParH(lev)->meanVelocityInZdirectionOut[pos] = (real)0.0;
            para->getParH(lev)->meanDensityOut[pos] = (real)0.0;
            para->getParH(lev)->meanPressureOut[pos] = (real)0.0;
            para->getParH(lev)->meanConcentrationOut[pos] = (real)0.0;
        }
    }
}

void calcMeanAD(Parameter* para, uint tdiff)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        for (size_t pos = 0; pos < para->getParH(lev)->numberOfNodes; pos++) {
            para->getParH(lev)->meanVelocityInXdirectionOut[pos] = para->getParH(lev)->meanVelocityInXdirection[pos] / (real)tdiff;
            para->getParH(lev)->meanVelocityInYdirectionOut[pos] = para->getParH(lev)->meanVelocityInYdirection[pos] / (real)tdiff;
            para->getParH(lev)->meanVelocityInZdirectionOut[pos] = para->getParH(lev)->meanVelocityInZdirection[pos] / (real)tdiff;
            para->getParH(lev)->meanDensityOut[pos] = para->getParH(lev)->meanDensity[pos] / (real)tdiff;
            para->getParH(lev)->meanPressureOut[pos] = para->getParH(lev)->meanPressure[pos] / (real)tdiff;
            para->getParH(lev)->meanConcentrationOut[pos] = para->getParH(lev)->meanConcentration[pos] / (real)tdiff;
        }
    }
}

void resetMeanAD(Parameter* para)
{
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++) {
        ResetMeanValuesAD27(para->getParD(lev)->meanVelocityInXdirection, para->getParD(lev)->meanVelocityInYdirection, para->getParD(lev)->meanVelocityInZdirection,
                            para->getParD(lev)->meanDensity, para->getParD(lev)->meanPressure, para->getParD(lev)->meanConcentration,
                            para->getParD(lev)->numberOfNodes, para->getParD(lev)->numberofthreads,
                            para->getParD(lev)->isEvenTimestep);
        getLastCudaError("ResetMeanValuesAD27 execution failed");
    }
}

//! \}
