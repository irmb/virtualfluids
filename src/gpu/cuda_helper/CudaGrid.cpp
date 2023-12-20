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
//! \addtogroup gpu_cuda_helper cuda_helper
//! \{
//! \author Soeren Peters
//=======================================================================================
#include "CudaGrid.h"

#include <logger/Logger.h>

namespace vf::cuda
{

CudaGrid::CudaGrid(unsigned int numberOfThreads, unsigned int numberOfEntities)
    : threads { dim3(numberOfThreads, 1, 1) }, grid { getCudaGrid(numberOfThreads, numberOfEntities) }
{
}
void CudaGrid::print() const
{
    VF_LOG_INFO("blocks: ({},{},{}), threads: ({},{},{})", grid.x, grid.y, grid.z, threads.x, threads.y, threads.z);
}

dim3 getCudaGrid(unsigned int numberOfThreads, unsigned int numberOfEntities)
{
    unsigned int Grid = (numberOfEntities / numberOfThreads) + 1;
    unsigned int Grid1, Grid2;
    if (Grid > 512) {
        Grid1 = 512;
        Grid2 = (Grid / Grid1) + 1;
    } else {
        Grid1 = 1;
        Grid2 = Grid;
    }
    dim3 cudaGrid(Grid1, Grid2);
    return cudaGrid;
}

} // namespace vf::cuda

//! \}
