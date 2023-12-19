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
#include "K15CompressibleNavierStokes.h"

#include "K15CompressibleNavierStokes_Device.cuh"

#include "Parameter/Parameter.h"

std::shared_ptr<K15CompressibleNavierStokes> K15CompressibleNavierStokes::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<K15CompressibleNavierStokes>(new K15CompressibleNavierStokes(para, level));
}

void K15CompressibleNavierStokes::run()
{
    int numberOfThreads = para->getParD(level)->numberofthreads;
    int size_Mat = (int)para->getParD(level)->numberOfNodes;

    int Grid = (size_Mat / numberOfThreads) + 1;
    int Grid1, Grid2;
    if (Grid > 512)
    {
        Grid1 = 512;
        Grid2 = (Grid / Grid1) + 1;
    }
    else
    {
        Grid1 = 1;
        Grid2 = Grid;
    }
    dim3 grid(Grid1, Grid2, 1);
    dim3 threads(numberOfThreads, 1, 1);

    K15CompressibleNavierStokes_Device <<< grid, threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        level,
        para->getForcesDev(),
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_CumulantK15Comp execution failed");
}

K15CompressibleNavierStokes::K15CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitNavierStokesCompressible);

    
}
//! \}
