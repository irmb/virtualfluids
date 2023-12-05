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
//! \author Anna Wellmann
//=======================================================================================
#include "TurbulenceIntensity.h"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/MacroscopicQuantities.h>

#include "LBM/GPUHelperFunctions/KernelUtilities.h"

__global__ void CalcTurbulenceIntensity(real* vxx, real* vyy, real* vzz, real* vxy, real* vxz, real* vyz, real* vx_mean,
                                        real* vy_mean, real* vz_mean, real* distributions, uint* typeOfGridNode,
                                        unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                                        unsigned long long numberOfLBnodes, bool isEvenTimestep)
{
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if (nodeIndex >= numberOfLBnodes)
        return;

    if (!vf::gpu::isValidFluidNode(typeOfGridNode[nodeIndex]))
        return;

    Distributions27 dist;
    vf::gpu::getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
    vf::gpu::ListIndices listIndices(nodeIndex, neighborX, neighborY, neighborZ);

    real distribution[27];
    vf::gpu::getPreCollisionDistribution(distribution, dist, listIndices);

    // analogue to LBCalcMacCompSP27
    real rho = vf::lbm::getDensity(distribution);
    real vx = vf::lbm::getCompressibleVelocityX1(distribution, rho);
    real vy = vf::lbm::getCompressibleVelocityX2(distribution, rho);
    real vz = vf::lbm::getCompressibleVelocityX3(distribution, rho);

    // compute subtotals:
    // fluctuations
    vxx[nodeIndex] = vxx[nodeIndex] + vx * vx;
    vyy[nodeIndex] = vyy[nodeIndex] + vy * vy;
    vzz[nodeIndex] = vzz[nodeIndex] + vz * vz;
    vxy[nodeIndex] = vxy[nodeIndex] + vx * vy;
    vxz[nodeIndex] = vxz[nodeIndex] + vx * vz;
    vyz[nodeIndex] = vyz[nodeIndex] + vy * vz;

    // velocity (for mean velocity)
    vx_mean[nodeIndex] = vx_mean[nodeIndex] + vx;
    vy_mean[nodeIndex] = vy_mean[nodeIndex] + vy;
    vz_mean[nodeIndex] = vz_mean[nodeIndex] + vz;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CalcTurbulenceIntensityDevice(real* vxx, real* vyy, real* vzz, real* vxy, real* vxz, real* vyz, real* vx_mean,
                                   real* vy_mean, real* vz_mean, real* DD, uint* typeOfGridNode, unsigned int* neighborX,
                                   unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                                   bool isEvenTimestep, uint numberOfThreads)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);
    CalcTurbulenceIntensity<<<grid.grid, grid.threads>>>(vxx, vyy, vzz, vxy, vxz, vyz, vx_mean, vy_mean, vz_mean, DD,
                                                         typeOfGridNode, neighborX, neighborY, neighborZ, numberOfLBnodes,
                                                         isEvenTimestep);
    getLastCudaError("CalcTurbulenceIntensity execution failed");
}
