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

#include "KineticEnergyAnalyzer.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cmath>
#include <iomanip>
#include <sstream>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>

#include <cuda_helper/CudaGrid.h>

#include <basics/constants/NumericConstants.h>

#include "GPU/GPU_Kernels.cuh"
#include "Parameter/Parameter.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

__global__                 void kineticEnergyKernel  (real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, unsigned long long numberOfLBnodes);

__host__ __device__ inline void kineticEnergyFunction(real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, uint index);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool KineticEnergyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

    int lev = 0;
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(lev)->numberofthreads, para->getParD(lev)->numberOfNodes);

    thrust::device_vector<real> kineticEnergy( this->para->getParD(lev)->numberOfNodes, c0o1);
    thrust::device_vector<uint> isFluid      ( this->para->getParD(lev)->numberOfNodes, 0);

    LBCalcMacCompSP27<<< grid.grid, grid.threads >>>(
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
        para->getParD(lev)->distributions.f[0],
        para->getParD(lev)->isEvenTimestep); 
    getLastCudaError("LBCalcMacCompSP27 execution failed"); 

    kineticEnergyKernel<<< grid.grid, grid.threads >>>(
        para->getParD(lev)->velocityX, 
        para->getParD(lev)->velocityY, 
        para->getParD(lev)->velocityZ, 
        para->getParD(lev)->rho, 
        para->getParD(lev)->neighborX,
        para->getParD(lev)->neighborY,
        para->getParD(lev)->neighborZ,
        para->getParD(lev)->neighborInverse,
        para->getParD(lev)->typeOfGridNode,
        kineticEnergy.data().get(), 
        isFluid.data().get(),
        para->getParD(lev)->numberOfNodes);
    cudaDeviceSynchronize();

    getLastCudaError("kineticEnergyKernel execution failed");

     real EKin               = thrust::reduce(kineticEnergy.begin(), kineticEnergy.end(), c0o1, thrust::plus<real>());
     uint numberOfFluidNodes = thrust::reduce(isFluid.begin(),       isFluid.end(),       0,    thrust::plus<uint>());

    //std::cout << "EKin " << EKin << "   " << numberOfFluidNodes << std::endl;

    this->kineticEnergyTimeSeries.push_back( EKin / real(numberOfFluidNodes) );

    //TODO: Should this function probably return nothing?
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void kineticEnergyKernel(real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, unsigned long long numberOfLBnodes)
{
    //////////////////////////////////////////////////////////////////////////
    const uint x = threadIdx.x;  // Globaler x-Index 
    const uint y = blockIdx.x;   // Globaler y-Index 
    const uint z = blockIdx.y;   // Globaler z-Index 

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint index = nx*(ny*z + y) + x;
    ////////////////////////////////////////////////////////////////////////////////
    //printf("%d\n", index);

    //if( index % 34 == 0 || index % 34 == 33 ) return;

    if( index >= (uint)numberOfLBnodes) return;

    unsigned int BC;
    BC = geo[index];
    if (BC != GEO_FLUID) return;

    kineticEnergyFunction( vx, vy, vz, rho, neighborX, neighborY, neighborZ, neighborWSB, geo, kineticEnergy, isFluid, index );
}

__host__ __device__ void kineticEnergyFunction(real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, uint index)
{
    //////////////////////////////////////////////////////////////////////////////
    //neighbor index                                
    uint k     = index;                             
    uint kPx   = neighborX[k];                      if( geo[ kPx   ] != GEO_FLUID ) return;
    uint kPy   = neighborY[k];                      if( geo[ kPy   ] != GEO_FLUID ) return;
    uint kPz   = neighborZ[k];                      if( geo[ kPz   ] != GEO_FLUID ) return;
    uint kMxyz = neighborWSB[k];                    if( geo[ kMxyz ] != GEO_FLUID ) return;
    uint kMx   = neighborZ[neighborY[kMxyz]];       if( geo[ kMx   ] != GEO_FLUID ) return;
    uint kMy   = neighborZ[neighborX[kMxyz]];       if( geo[ kMy   ] != GEO_FLUID ) return;
    uint kMz   = neighborY[neighborX[kMxyz]];       if( geo[ kMz   ] != GEO_FLUID ) return;
    //////////////////////////////////////////////////////////////////////////

    isFluid[ index ] = 1;

    kineticEnergy[ index ] = c1o2 * ( vx[index] * vx[index] + vy[index] * vy[index] + vz[index] * vz[index] ) * (rho[index] + c1o1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KineticEnergyAnalyzer::KineticEnergyAnalyzer(SPtr<Parameter> para, uint analyzeIter)
{
    this->para = para;
    this->analyzeIter = analyzeIter;
}

void KineticEnergyAnalyzer::writeToFile(std::string filename)
{
    std::cout << "KineticEnergyAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + "_KineticEnergyData.dat" );

    for( auto& EKin : this->kineticEnergyTimeSeries )
        file << std::setprecision(15) << EKin << std::endl;

    file.close();

    std::cout << "done!\n";
}


