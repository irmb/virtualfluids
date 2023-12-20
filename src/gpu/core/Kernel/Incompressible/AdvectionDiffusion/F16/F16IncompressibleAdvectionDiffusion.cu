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
//! \addtogroup gpu_Kernel Kernel
//! \ingroup gpu_core core
//! \{
#include "F16IncompressibleAdvectionDiffusion.h"

#include "F16IncompressibleAdvectionDiffusion_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<F16IncompressibleAdvectionDiffusion> F16IncompressibleAdvectionDiffusion::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<F16IncompressibleAdvectionDiffusion>(new F16IncompressibleAdvectionDiffusion(para, level));
}

void F16IncompressibleAdvectionDiffusion::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    F16IncompressibleAdvectionDiffusion_Device<<< grid.grid, grid.threads >>>(
        para->getParD(level)->diffusivity, 
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY, 
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0], 
        para->getParD(level)->distributionsAD.f[0], 
        para->getParD(level)->numberOfNodes, 
        para->getParD(level)->forcing,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("F16IncompressibleAdvectionDiffusion_Device execution failed");
}

F16IncompressibleAdvectionDiffusion::F16IncompressibleAdvectionDiffusion(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitAdvectionDiffusionIncompressible);

}

F16IncompressibleAdvectionDiffusion::F16IncompressibleAdvectionDiffusion()
{
}

//! \}
