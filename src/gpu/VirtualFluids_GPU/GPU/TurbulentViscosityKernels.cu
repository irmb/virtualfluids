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
//! \file TurbulentViscosityKernels.cu
//! \ingroup GPU
//! \author Henry Korb, Henrik Asmuth
//======================================================================================

#include "TurbulentViscosityKernels.h"
#include "lbm/constants/NumericConstants.h"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "LBM/LB.h"

using namespace vf::lbm::constant;

__host__ __device__ __forceinline__ void calcDerivatives(const uint& k, uint& kM, uint& kP, uint* typeOfGridNode, real* vx, real* vy, real* vz, real& dvx, real& dvy, real& dvz)
{
    bool fluidP = (typeOfGridNode[kP] == GEO_FLUID);
    bool fluidM = (typeOfGridNode[kM] == GEO_FLUID);
    real div = (fluidM & fluidP) ? c1o2 : c1o1;

    dvx = ((fluidP ? vx[kP] : vx[k])-(fluidM ? vx[kM] : vx[k]))*div;
    dvy = ((fluidP ? vy[kP] : vy[k])-(fluidM ? vy[kM] : vy[k]))*div;
    dvz = ((fluidP ? vz[kP] : vz[k])-(fluidM ? vz[kM] : vz[k]))*div;
}

__global__ void calcAMD(real* vx,
                        real* vy,
                        real* vz,
                        real* turbulentViscosity,
                        uint* neighborX,
                        uint* neighborY,
                        uint* neighborZ,
                        uint* neighborWSB,
                        uint* typeOfGridNode,
                        uint size_Mat,
                        real SGSConstant)
{

    const uint x = threadIdx.x; 
    const uint y = blockIdx.x; 
    const uint z = blockIdx.y; 

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint k = nx*(ny*z + y) + x;
    if(k >= size_Mat) return;
    if(typeOfGridNode[k] != GEO_FLUID) return;

    uint kPx = neighborX[k];
    uint kPy = neighborY[k];
    uint kPz = neighborZ[k];
    uint kMxyz = neighborWSB[k];
    uint kMx = neighborZ[neighborY[kMxyz]];
    uint kMy = neighborZ[neighborX[kMxyz]];
    uint kMz = neighborY[neighborX[kMxyz]];

    real dvxdx, dvxdy, dvxdz,
         dvydx, dvydy, dvydz,
         dvzdx, dvzdy, dvzdz;

    calcDerivatives(k, kMx, kPx, typeOfGridNode, vx, vy, vz, dvxdx, dvydx, dvzdx);
    calcDerivatives(k, kMy, kPy, typeOfGridNode, vx, vy, vz, dvxdy, dvydy, dvzdy);
    calcDerivatives(k, kMz, kPz, typeOfGridNode, vx, vy, vz, dvxdz, dvydz, dvzdz);

    real denominator =  dvxdx*dvxdx + dvydx*dvydx + dvzdx*dvzdx + 
                        dvxdy*dvxdy + dvydy*dvydy + dvzdy*dvzdy +
                        dvxdz*dvxdz + dvydz*dvydz + dvzdz*dvzdz;
    real enumerator =   (dvxdx*dvxdx + dvxdy*dvxdy + dvxdz*dvxdz) * dvxdx + 
                        (dvydx*dvydx + dvydy*dvydy + dvydz*dvydz) * dvydy + 
                        (dvzdx*dvzdx + dvzdy*dvzdy + dvzdz*dvzdz) * dvzdz +
                        (dvxdx*dvydx + dvxdy*dvydy + dvxdz*dvydz) * (dvxdy+dvydx) +
                        (dvxdx*dvzdx + dvxdy*dvzdy + dvxdz*dvzdz) * (dvxdz+dvzdx) + 
                        (dvydx*dvzdx + dvydy*dvzdy + dvydz*dvzdz) * (dvydz+dvzdy);

    turbulentViscosity[k] = max(c0o1,-SGSConstant*enumerator)/denominator;
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