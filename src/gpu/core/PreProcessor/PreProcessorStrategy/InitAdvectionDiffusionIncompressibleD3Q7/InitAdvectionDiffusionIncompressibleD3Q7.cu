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
#include "InitAdvectionDiffusionIncompressibleD3Q7.h"

#include "InitAdvectionDiffusionIncompressibleD3Q7_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<PreProcessorStrategy> InitAdvectionDiffusionIncompressibleD3Q7::getNewInstance(std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<PreProcessorStrategy>(new InitAdvectionDiffusionIncompressibleD3Q7(para));
}

void InitAdvectionDiffusionIncompressibleD3Q7::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    InitAdvectionDiffusionIncompressibleD3Q7_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->concentration,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->distributionsAD.f[0],
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("InitAdvectionDiffusionIncompressibleD3Q7_Device execution failed");
}

bool InitAdvectionDiffusionIncompressibleD3Q7::checkParameter()
{
    return false;
}

InitAdvectionDiffusionIncompressibleD3Q7::InitAdvectionDiffusionIncompressibleD3Q7(std::shared_ptr<Parameter> para)
{
    this->para = para;
}

InitAdvectionDiffusionIncompressibleD3Q7::InitAdvectionDiffusionIncompressibleD3Q7()
{
}
