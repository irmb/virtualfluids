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
//! \author Henry Korb, Henrik Asmuth
//======================================================================================

#include "TurbulentViscosityKernels.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

#include "Utilities/KernelUtilities.h"
#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"

using namespace vf::basics::constant;

__host__ __device__ __forceinline__ void calcDerivatives(const uint& k, uint& kM, uint& kP, uint* typeOfGridNode, real* vx, real* vy, real* vz, real& dvx, real& dvy, real& dvz)
{
    bool fluidP = (typeOfGridNode[kP] == GEO_FLUID);
    bool fluidM = (typeOfGridNode[kM] == GEO_FLUID);
    real div = (fluidM & fluidP) ? c1o2 : c1o1;

    dvx = ((fluidP ? vx[kP] : vx[k])-(fluidM ? vx[kM] : vx[k]))*div;
    dvy = ((fluidP ? vy[kP] : vy[k])-(fluidM ? vy[kM] : vy[k]))*div;
    dvz = ((fluidP ? vz[kP] : vz[k])-(fluidM ? vz[kM] : vz[k]))*div;
}

__global__ void calcAMD(
    real* vx,
    real* vy,
    real* vz,
    real* turbulentViscosity,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    uint* neighborWSB,
    uint* typeOfGridNode,
    unsigned long long numberOfLBnodes,
    real SGSConstant)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if(nodeIndex >= numberOfLBnodes) return;
    if(typeOfGridNode[nodeIndex] != GEO_FLUID) return;

    uint kPx = neighborX[nodeIndex];
    uint kPy = neighborY[nodeIndex];
    uint kPz = neighborZ[nodeIndex];
    uint kMxyz = neighborWSB[nodeIndex];
    uint kMx = neighborZ[neighborY[kMxyz]];
    uint kMy = neighborZ[neighborX[kMxyz]];
    uint kMz = neighborY[neighborX[kMxyz]];

    real dvxdx, dvxdy, dvxdz,
         dvydx, dvydy, dvydz,
         dvzdx, dvzdy, dvzdz;

    calcDerivatives(nodeIndex, kMx, kPx, typeOfGridNode, vx, vy, vz, dvxdx, dvydx, dvzdx);
    calcDerivatives(nodeIndex, kMy, kPy, typeOfGridNode, vx, vy, vz, dvxdy, dvydy, dvzdy);
    calcDerivatives(nodeIndex, kMz, kPz, typeOfGridNode, vx, vy, vz, dvxdz, dvydz, dvzdz);

    real denominator =  dvxdx*dvxdx + dvydx*dvydx + dvzdx*dvzdx + 
                        dvxdy*dvxdy + dvydy*dvydy + dvzdy*dvzdy +
                        dvxdz*dvxdz + dvydz*dvydz + dvzdz*dvzdz;
    real enumerator =   (dvxdx*dvxdx + dvxdy*dvxdy + dvxdz*dvxdz) * dvxdx + 
                        (dvydx*dvydx + dvydy*dvydy + dvydz*dvydz) * dvydy + 
                        (dvzdx*dvzdx + dvzdy*dvzdy + dvzdz*dvzdz) * dvzdz +
                        (dvxdx*dvydx + dvxdy*dvydy + dvxdz*dvydz) * (dvxdy+dvydx) +
                        (dvxdx*dvzdx + dvxdy*dvzdy + dvxdz*dvzdz) * (dvxdz+dvzdx) + 
                        (dvydx*dvzdx + dvydy*dvzdy + dvydz*dvzdz) * (dvydz+dvzdy);

    turbulentViscosity[nodeIndex] = denominator != c0o1 ? max(c0o1,-SGSConstant*enumerator)/denominator : c0o1;
}

void calcTurbulentViscosityAMD(Parameter* para, int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, para->getParH(level)->numberOfNodes);
    calcAMD<<<grid.grid, grid.threads>>>(
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->turbViscosity,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->neighborInverse,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->numberOfNodes,
        para->getSGSConstant()
    );
    getLastCudaError("calcAMD execution failed");
}