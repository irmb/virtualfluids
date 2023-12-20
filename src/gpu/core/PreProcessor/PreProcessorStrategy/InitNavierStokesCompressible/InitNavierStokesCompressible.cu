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
#include "InitNavierStokesCompressible.h"

#include "InitNavierStokesCompressible_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<PreProcessorStrategy> InitNavierStokesCompressible::getNewInstance(std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<PreProcessorStrategy>(new InitNavierStokesCompressible(para));
}

void InitNavierStokesCompressible::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    if( ! para->getUseInitNeq() )
    {
        InitNavierStokesCompressible_Device <<< grid.grid, grid.threads >>> (
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->typeOfGridNode,
            para->getParD(level)->rho,
            para->getParD(level)->velocityX,
            para->getParD(level)->velocityY,
            para->getParD(level)->velocityZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->isEvenTimestep);
        getLastCudaError("InitNavierStokesCompressible_Device execution failed");
    }
    else
    {
        InitNavierStokesCompressibleNonEquilibrium_Device <<< grid.grid, grid.threads >>> (
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->neighborInverse,
            para->getParD(level)->typeOfGridNode,
            para->getParD(level)->rho,
            para->getParD(level)->velocityX,
            para->getParD(level)->velocityY,
            para->getParD(level)->velocityZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->omega,
            para->getParD(level)->isEvenTimestep);
        cudaDeviceSynchronize();
        getLastCudaError("InitNavierStokesCompressibleNonEquilibrium_Device execution failed");
    }



}

bool InitNavierStokesCompressible::checkParameter()
{
    return false;
}

InitNavierStokesCompressible::InitNavierStokesCompressible(std::shared_ptr<Parameter> para)
{
    this->para = para;
}

InitNavierStokesCompressible::InitNavierStokesCompressible()
{
}

//! \}
