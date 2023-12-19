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
//! \addtogroup gpu_PreProcessor PreProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "PreProcessor/InitLattice.h"

#include "Cuda/CudaMemoryManager.h"

#include "Parameter/Parameter.h"
#include "PostProcessor/MacroscopicQuantities.cuh"
#include "PreProcessor/PreProcessor.h"

void initLattice(SPtr<Parameter> para, SPtr<PreProcessor> preProcessor, SPtr<PreProcessor> preProcessorAD, SPtr<CudaMemoryManager> cudaMemoryManager)
{
    for (int lev = para->getFine(); lev >= para->getCoarse(); lev--) {
        preProcessor->init(para, lev);

        calculateMacroscopicQuantitiesCompressible(
            para->getParD(lev)->velocityX, 
            para->getParD(lev)->velocityY, 
            para->getParD(lev)->velocityZ, 
            para->getParD(lev)->rho,
            para->getParD(lev)->pressure, 
            para->getParD(lev)->typeOfGridNode, 
            para->getParD(lev)->neighborX,
            para->getParD(lev)->neighborY, 
            para->getParD(lev)->neighborZ, 
            para->getParD(lev)->numberOfNodes,
            para->getParD(lev)->numberofthreads, 
            para->getParD(lev)->distributions.f[0], 
            para->getParD(lev)->isEvenTimestep);

        if (para->getCalcMean()) {
            constexpr uint tdiff = 1;
            calculateMacrosopicMean(
                para->getParD(lev)->meanVelocityInXdirection, 
                para->getParD(lev)->meanVelocityInYdirection, 
                para->getParD(lev)->meanVelocityInZdirection,
                para->getParD(lev)->meanDensity, 
                para->getParD(lev)->meanPressure, 
                para->getParD(lev)->typeOfGridNode,
                para->getParD(lev)->neighborX, 
                para->getParD(lev)->neighborY,
                para->getParD(lev)->neighborZ, 
                tdiff, 
                para->getParD(lev)->numberOfNodes,
                para->getParD(lev)->numberofthreads, 
                para->getParD(lev)->isEvenTimestep);
        }
        // advection - diffusion
        if (para->getDiffOn()) {

            cudaMemoryManager->cudaAllocConcentration(lev);

            for (size_t index = 0; index < para->getParH(lev)->numberOfNodes; index++) {
                para->getParH(lev)->concentration[index] = para->getConcentrationInit();
            }

            preProcessorAD->init(para, lev);
        }
    }
}

//! \}
