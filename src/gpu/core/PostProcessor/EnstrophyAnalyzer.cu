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
//! \addtogroup gpu_PostProcessor PostProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================

#include "EnstrophyAnalyzer.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cmath>
#include <iomanip>
#include <sstream>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>

#include <cuda_helper/CudaGrid.h>

#include <basics/constants/NumericConstants.h>

#include "Parameter/Parameter.h"
#include "PostProcessor/MacroscopicQuantities.cuh"

using namespace vf::basics::constant;

__global__                 void enstrophyKernel  ( real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* enstrophy, uint* isFluid, unsigned long long numberOfLBnodes );

__host__ __device__ inline void enstrophyFunction( real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* enstrophy, uint* isFluid, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool EnstrophyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

    int lev = 0;
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(lev)->numberofthreads, para->getParD(lev)->numberOfNodes);

    thrust::device_vector<real> enstrophy( this->para->getParD(lev)->numberOfNodes, c0o1);
    thrust::device_vector<uint> isFluid  ( this->para->getParD(lev)->numberOfNodes, 0);

    calculateMacroscopicQuantitiesCompressible(para->getParD(lev)->velocityX, para->getParD(lev)->velocityY, para->getParD(lev)->velocityZ,
                    para->getParD(lev)->rho, para->getParD(lev)->pressure, para->getParD(lev)->typeOfGridNode,
                    para->getParD(lev)->neighborX, para->getParD(lev)->neighborY, para->getParD(lev)->neighborZ,
                    para->getParD(lev)->numberOfNodes, para->getParD(lev)->numberofthreads, para->getParD(lev)->distributions.f[0],
                    para->getParD(lev)->isEvenTimestep);

    enstrophyKernel<<< grid.grid, grid.threads >>>(
        para->getParD(lev)->velocityX,
        para->getParD(lev)->velocityY, 
        para->getParD(lev)->velocityZ, 
        para->getParD(lev)->rho, 
        para->getParD(lev)->neighborX,
        para->getParD(lev)->neighborY,
        para->getParD(lev)->neighborZ,
        para->getParD(lev)->neighborInverse,
        para->getParD(lev)->typeOfGridNode,
        enstrophy.data().get(), 
        isFluid.data().get(),
        para->getParD(lev)->numberOfNodes);
    cudaDeviceSynchronize(); 
    getLastCudaError("enstrophyKernel execution failed");

    real EnstrophyTmp       = thrust::reduce(enstrophy.begin(), enstrophy.end(), c0o1, thrust::plus<real>());
    uint numberOfFluidNodes = thrust::reduce(isFluid.begin(),   isFluid.end(),   0,    thrust::plus<uint>());

    //std::cout << "Enstrophy " << EnstrophyTmp << "   " << numberOfFluidNodes << std::endl;

    this->enstrophyTimeSeries.push_back( EnstrophyTmp / real(numberOfFluidNodes) );

    //TODO: Should this function probably return nothing?
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void enstrophyKernel(real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* enstrophy, uint* isFluid, unsigned long long numberOfLBnodes)
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

    enstrophyFunction( veloX, veloY, veloZ, rho, neighborX, neighborY, neighborZ, neighborWSB, geo, enstrophy, isFluid, index );
}

__host__ __device__ void enstrophyFunction(real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* enstrophy, uint* isFluid, uint index)
{
    uint maxOrderX = 8;
    uint maxOrderY = 8;
    uint maxOrderZ = 8;

    //////////////////////////////////////////////////////////////////////////
    //neighbor index                                
    uint k       = index;                             
    uint kPx     = neighborX[k];                    if( geo[ kPx     ] != GEO_FLUID ) return;
    uint kPy     = neighborY[k];                    if( geo[ kPy     ] != GEO_FLUID ) return;
    uint kPz     = neighborZ[k];                    if( geo[ kPz     ] != GEO_FLUID ) return;
    uint kMxyz   = neighborWSB[k];                    
    uint kMx     = neighborZ[neighborY[kMxyz]];     if( geo[ kMx     ] != GEO_FLUID ) return;
    uint kMy     = neighborZ[neighborX[kMxyz]];     if( geo[ kMy     ] != GEO_FLUID ) return;
    uint kMz     = neighborY[neighborX[kMxyz]];     if( geo[ kMz     ] != GEO_FLUID ) return;
    //////////////////////////////////////////////////////////////////////////
    uint kPx2    = neighborX[kPx];                  if( geo[ kPx2    ] != GEO_FLUID ) maxOrderX = 2;
    uint kPy2    = neighborY[kPy];                  if( geo[ kPy2    ] != GEO_FLUID ) maxOrderY = 2;
    uint kPz2    = neighborZ[kPz];                  if( geo[ kPz2    ] != GEO_FLUID ) maxOrderZ = 2;
    uint kMxWSB  = neighborWSB[kMx];                 
    uint kMyWSB  = neighborWSB[kMy];                 
    uint kMzWSB  = neighborWSB[kMz];                 
    uint kMx2    = neighborZ[neighborY[kMxWSB]];    if( geo[ kMx2    ] != GEO_FLUID ) maxOrderX = 2;
    uint kMy2    = neighborZ[neighborX[kMyWSB]];    if( geo[ kMy2    ] != GEO_FLUID ) maxOrderY = 2;
    uint kMz2    = neighborY[neighborX[kMzWSB]];    if( geo[ kMz2    ] != GEO_FLUID ) maxOrderZ = 2;
    //////////////////////////////////////////////////////////////////////////
    uint kPx3    = neighborX[kPx2];                 if( geo[ kPx3    ] != GEO_FLUID && maxOrderX > 4 ) maxOrderX = 4;
    uint kPy3    = neighborY[kPy2];                 if( geo[ kPy3    ] != GEO_FLUID && maxOrderY > 4 ) maxOrderY = 4;
    uint kPz3    = neighborZ[kPz2];                 if( geo[ kPz3    ] != GEO_FLUID && maxOrderZ > 4 ) maxOrderZ = 4;
    uint kMx2WSB = neighborWSB[kMx2];               
    uint kMy2WSB = neighborWSB[kMy2];               
    uint kMz2WSB = neighborWSB[kMz2];               
    uint kMx3    = neighborZ[neighborY[kMx2WSB]];   if( geo[ kMx3    ] != GEO_FLUID && maxOrderX > 4 ) maxOrderX = 4;
    uint kMy3    = neighborZ[neighborX[kMy2WSB]];   if( geo[ kMy3    ] != GEO_FLUID && maxOrderY > 4 ) maxOrderY = 4;
    uint kMz3    = neighborY[neighborX[kMz2WSB]];   if( geo[ kMz3    ] != GEO_FLUID && maxOrderZ > 4 ) maxOrderZ = 4;
    //////////////////////////////////////////////////////////////////////////
    uint kPx4    = neighborX[kPx3];                 if( geo[ kPx4    ] != GEO_FLUID && maxOrderX > 6 ) maxOrderX = 6;
    uint kPy4    = neighborY[kPy3];                 if( geo[ kPy4    ] != GEO_FLUID && maxOrderY > 6 ) maxOrderY = 6;
    uint kPz4    = neighborZ[kPz3];                 if( geo[ kPz4    ] != GEO_FLUID && maxOrderZ > 6 ) maxOrderZ = 6;
    uint kMx3WSB = neighborWSB[kMx3];               
    uint kMy3WSB = neighborWSB[kMy3];               
    uint kMz3WSB = neighborWSB[kMz3];               
    uint kMx4    = neighborZ[neighborY[kMx3WSB]];   if( geo[ kMx4    ] != GEO_FLUID && maxOrderX > 6 ) maxOrderX = 6;
    uint kMy4    = neighborZ[neighborX[kMy3WSB]];   if( geo[ kMy4    ] != GEO_FLUID && maxOrderY > 6 ) maxOrderY = 6;
    uint kMz4    = neighborY[neighborX[kMz3WSB]];   if( geo[ kMz4    ] != GEO_FLUID && maxOrderZ > 6 ) maxOrderZ = 6;
    //////////////////////////////////////////////////////////////////////////
    //getVeloX//
    //real veloXNeighborPx = veloX[kPx];
    //real veloXNeighborMx = veloX[kMx];
    real veloXNeighborPy = veloX[kPy];
    real veloXNeighborMy = veloX[kMy];
    real veloXNeighborPz = veloX[kPz];
    real veloXNeighborMz = veloX[kMz];
    //getVeloY//
    real veloYNeighborPx = veloY[kPx];
    real veloYNeighborMx = veloY[kMx];
    //real veloYNeighborPy = veloY[kPy];
    //real veloYNeighborMy = veloY[kMy];
    real veloYNeighborPz = veloY[kPz];
    real veloYNeighborMz = veloY[kMz];
    //getVeloZ//
    real veloZNeighborPx = veloZ[kPx];
    real veloZNeighborMx = veloZ[kMx];
    real veloZNeighborPy = veloZ[kPy];
    real veloZNeighborMy = veloZ[kMy];
    //real veloZNeighborPz = veloZ[kPz];
    //real veloZNeighborMz = veloZ[kMz];
    //////////////////////////////////////////////////////////////////////////////
    //getVeloX//
    //real veloXNeighborPx2 = veloX[kPx2];
    //real veloXNeighborMx2 = veloX[kMx2];
    real veloXNeighborPy2 = veloX[kPy2];
    real veloXNeighborMy2 = veloX[kMy2];
    real veloXNeighborPz2 = veloX[kPz2];
    real veloXNeighborMz2 = veloX[kMz2];
    //getVeloY//
    real veloYNeighborPx2 = veloY[kPx2];
    real veloYNeighborMx2 = veloY[kMx2];
    //real veloYNeighborPy2 = veloY[kPy2];
    //real veloYNeighborMy2 = veloY[kMy2];
    real veloYNeighborPz2 = veloY[kPz2];
    real veloYNeighborMz2 = veloY[kMz2];
    //getVeloZ//
    real veloZNeighborPx2 = veloZ[kPx2];
    real veloZNeighborMx2 = veloZ[kMx2];
    real veloZNeighborPy2 = veloZ[kPy2];
    real veloZNeighborMy2 = veloZ[kMy2];
    //real veloZNeighborPz2 = veloZ[kPz2];
    //real veloZNeighborMz2 = veloZ[kMz2];
    //////////////////////////////////////////////////////////////////////////////
    //getVeloX//
    //real veloXNeighborPx3 = veloX[kPx3];
    //real veloXNeighborMx3 = veloX[kMx3];
    real veloXNeighborPy3 = veloX[kPy3];
    real veloXNeighborMy3 = veloX[kMy3];
    real veloXNeighborPz3 = veloX[kPz3];
    real veloXNeighborMz3 = veloX[kMz3];
    //getVeloY//
    real veloYNeighborPx3 = veloY[kPx3];
    real veloYNeighborMx3 = veloY[kMx3];
    //real veloYNeighborPy3 = veloY[kPy3];
    //real veloYNeighborMy3 = veloY[kMy3];
    real veloYNeighborPz3 = veloY[kPz3];
    real veloYNeighborMz3 = veloY[kMz3];
    //getVeloZ//
    real veloZNeighborPx3 = veloZ[kPx3];
    real veloZNeighborMx3 = veloZ[kMx3];
    real veloZNeighborPy3 = veloZ[kPy3];
    real veloZNeighborMy3 = veloZ[kMy3];
    //real veloZNeighborPz3 = veloZ[kPz3];
    //real veloZNeighborMz3 = veloZ[kMz3];
    //////////////////////////////////////////////////////////////////////////////
    //getVeloX//
    //real veloXNeighborPx4 = veloX[kPx4];
    //real veloXNeighborMx4 = veloX[kMx4];
    real veloXNeighborPy4 = veloX[kPy4];
    real veloXNeighborMy4 = veloX[kMy4];
    real veloXNeighborPz4 = veloX[kPz4];
    real veloXNeighborMz4 = veloX[kMz4];
    //getVeloY//
    real veloYNeighborPx4 = veloY[kPx4];
    real veloYNeighborMx4 = veloY[kMx4];
    //real veloYNeighborPy4 = veloY[kPy4];
    //real veloYNeighborMy4 = veloY[kMy4];
    real veloYNeighborPz4 = veloY[kPz4];
    real veloYNeighborMz4 = veloY[kMz4];
    //getVeloZ//
    real veloZNeighborPx4 = veloZ[kPx4];
    real veloZNeighborMx4 = veloZ[kMx4];
    real veloZNeighborPy4 = veloZ[kPy4];
    real veloZNeighborMy4 = veloZ[kMy4];
    //real veloZNeighborPz4 = veloZ[kPz4];
    //real veloZNeighborMz4 = veloZ[kMz4];
    //////////////////////////////////////////////////////////////////////////////
    //real dxvx = c0o1;
    real dyvx = c0o1;
    real dzvx = c0o1;
    real dxvy = c0o1;
    //real dyvy = c0o1;
    real dzvy = c0o1;
    real dxvz = c0o1;
    real dyvz = c0o1;
    //real dzvz = c0o1;
    //////////////////////////////////////////////////////////////////////////

    //maxOrderX = 2;
    //maxOrderY = 2;
    //maxOrderZ = 2;

    //if( maxOrder == 2 )
    {
        if( maxOrderX == 2 ) dxvy = (veloYNeighborPx - veloYNeighborMx) / c2o1;
        if( maxOrderX == 2 ) dxvz = (veloZNeighborPx - veloZNeighborMx) / c2o1;
        if( maxOrderY == 2 ) dyvx = (veloXNeighborPy - veloXNeighborMy) / c2o1;
        if( maxOrderY == 2 ) dyvz = (veloZNeighborPy - veloZNeighborMy) / c2o1;
        if( maxOrderZ == 2 ) dzvx = (veloXNeighborPz - veloXNeighborMz) / c2o1;
        if( maxOrderZ == 2 ) dzvy = (veloYNeighborPz - veloYNeighborMz) / c2o1;
    }

    //////////////////////////////////////////////////////////////////////////

    //if( maxOrder == 4 )
    {
        if( maxOrderX == 4 ) dxvy = ((c8o1 * veloYNeighborPx - c8o1 * veloYNeighborMx) - (veloYNeighborPx2 - veloYNeighborMx2)) / c12o1;
        if( maxOrderX == 4 ) dxvz = ((c8o1 * veloZNeighborPx - c8o1 * veloZNeighborMx) - (veloZNeighborPx2 - veloZNeighborMx2)) / c12o1;
        if( maxOrderY == 4 ) dyvx = ((c8o1 * veloXNeighborPy - c8o1 * veloXNeighborMy) - (veloXNeighborPy2 - veloXNeighborMy2)) / c12o1;
        if( maxOrderY == 4 ) dyvz = ((c8o1 * veloZNeighborPy - c8o1 * veloZNeighborMy) - (veloZNeighborPy2 - veloZNeighborMy2)) / c12o1;
        if( maxOrderZ == 4 ) dzvx = ((c8o1 * veloXNeighborPz - c8o1 * veloXNeighborMz) - (veloXNeighborPz2 - veloXNeighborMz2)) / c12o1;
        if( maxOrderZ == 4 ) dzvy = ((c8o1 * veloYNeighborPz - c8o1 * veloYNeighborMz) - (veloYNeighborPz2 - veloYNeighborMz2)) / c12o1;

    }
    //////////////////////////////////////////////////////////////////////////

    //if( maxOrder == 6 )
    {
        if( maxOrderX == 6 ) dxvy = ((c5o1 * c9o1) * (veloYNeighborPx - veloYNeighborMx) - c9o1 * (veloYNeighborPx2 - veloYNeighborMx2) + (veloYNeighborPx3 - veloYNeighborMx3)) / (c6o1 * c10o1);
        if( maxOrderX == 6 ) dxvz = ((c5o1 * c9o1) * (veloZNeighborPx - veloZNeighborMx) - c9o1 * (veloZNeighborPx2 - veloZNeighborMx2) + (veloZNeighborPx3 - veloZNeighborMx3)) / (c6o1 * c10o1);
        if( maxOrderY == 6 ) dyvx = ((c5o1 * c9o1) * (veloXNeighborPy - veloXNeighborMy) - c9o1 * (veloXNeighborPy2 - veloXNeighborMy2) + (veloXNeighborPy3 - veloXNeighborMy3)) / (c6o1 * c10o1);
        if( maxOrderY == 6 ) dyvz = ((c5o1 * c9o1) * (veloZNeighborPy - veloZNeighborMy) - c9o1 * (veloZNeighborPy2 - veloZNeighborMy2) + (veloZNeighborPy3 - veloZNeighborMy3)) / (c6o1 * c10o1);
        if( maxOrderZ == 6 ) dzvx = ((c5o1 * c9o1) * (veloXNeighborPz - veloXNeighborMz) - c9o1 * (veloXNeighborPz2 - veloXNeighborMz2) + (veloXNeighborPz3 - veloXNeighborMz3)) / (c6o1 * c10o1);
        if( maxOrderZ == 6 ) dzvy = ((c5o1 * c9o1) * (veloYNeighborPz - veloYNeighborMz) - c9o1 * (veloYNeighborPz2 - veloYNeighborMz2) + (veloYNeighborPz3 - veloYNeighborMz3)) / (c6o1 * c10o1);
    }

    //////////////////////////////////////////////////////////////////////////

    //if( maxOrder == 8 )
    {
        if( maxOrderX == 8 ) dxvy = ((c28o1 * c8o1) * (veloYNeighborPx - veloYNeighborMx) - (c7o1 * c8o1) * (veloYNeighborPx2 - veloYNeighborMx2) + (c8o1 * c4o1 * c1o3) * (veloYNeighborPx3 - veloYNeighborMx3) - (veloYNeighborPx4 - veloYNeighborMx4)) / (c7o1 * c10o1 * c4o1);
        if( maxOrderX == 8 ) dxvz = ((c28o1 * c8o1) * (veloZNeighborPx - veloZNeighborMx) - (c7o1 * c8o1) * (veloZNeighborPx2 - veloZNeighborMx2) + (c8o1 * c4o1 * c1o3) * (veloZNeighborPx3 - veloZNeighborMx3) - (veloZNeighborPx4 - veloZNeighborMx4)) / (c7o1 * c10o1 * c4o1);
        if( maxOrderY == 8 ) dyvx = ((c28o1 * c8o1) * (veloXNeighborPy - veloXNeighborMy) - (c7o1 * c8o1) * (veloXNeighborPy2 - veloXNeighborMy2) + (c8o1 * c4o1 * c1o3) * (veloXNeighborPy3 - veloXNeighborMy3) - (veloXNeighborPy4 - veloXNeighborMy4)) / (c7o1 * c10o1 * c4o1);
        if( maxOrderY == 8 ) dyvz = ((c28o1 * c8o1) * (veloZNeighborPy - veloZNeighborMy) - (c7o1 * c8o1) * (veloZNeighborPy2 - veloZNeighborMy2) + (c8o1 * c4o1 * c1o3) * (veloZNeighborPy3 - veloZNeighborMy3) - (veloZNeighborPy4 - veloZNeighborMy4)) / (c7o1 * c10o1 * c4o1);
        if( maxOrderZ == 8 ) dzvx = ((c28o1 * c8o1) * (veloXNeighborPz - veloXNeighborMz) - (c7o1 * c8o1) * (veloXNeighborPz2 - veloXNeighborMz2) + (c8o1 * c4o1 * c1o3) * (veloXNeighborPz3 - veloXNeighborMz3) - (veloXNeighborPz4 - veloXNeighborMz4)) / (c7o1 * c10o1 * c4o1);
        if( maxOrderZ == 8 ) dzvy = ((c28o1 * c8o1) * (veloYNeighborPz - veloYNeighborMz) - (c7o1 * c8o1) * (veloYNeighborPz2 - veloYNeighborMz2) + (c8o1 * c4o1 * c1o3) * (veloYNeighborPz3 - veloYNeighborMz3) - (veloYNeighborPz4 - veloYNeighborMz4)) / (c7o1 * c10o1 * c4o1);
    }

    //////////////////////////////////////////////////////////////////////////

    real tmpX = dyvz - dzvy;
    real tmpY = dzvx - dxvz;
    real tmpZ = dxvy - dyvx;
    //////////////////////////////////////////////////////////////////////////

    isFluid[ index ] = 1;

    enstrophy[ index ] = c1o2 * (rho[index] + c1o1) * ( tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EnstrophyAnalyzer::EnstrophyAnalyzer(SPtr<Parameter> para, uint analyzeIter)
{
    this->para = para;
    this->analyzeIter = analyzeIter;
}

void EnstrophyAnalyzer::writeToFile( std::string filename )
{
    std::cout << "EnstrophyAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + "_EnstrophyData.dat" );

    for( auto& EKin : this->enstrophyTimeSeries )
        file << std::setprecision(15) << EKin << std::endl;

    file.close();

    std::cout << "done!\n";
}



//! \}
