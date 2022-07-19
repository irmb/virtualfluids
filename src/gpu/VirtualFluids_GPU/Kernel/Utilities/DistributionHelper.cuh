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
//! \file Cumulant27chim.cu
//! \ingroup GPU
//! \author Martin Schoenherr, Soeren Peters
//=======================================================================================
#ifndef DISTRIBUTUION_HELPER_CUH
#define DISTRIBUTUION_HELPER_CUH

#include "LBM/LB.h" 

#include "lbm/KernelParameter.h"
#include "lbm/constants/D3Q27.h"

using namespace vf::lbm::dir;

namespace vf::gpu
{

__inline__ __device__ __host__ void getPointersToDistributions(Distributions27 &dist, real *distributionArray, const uint numberOfLBnodes, const bool isEvenTimestep)
{
    if (isEvenTimestep)
    {
        dist.f[E   ] = &distributionArray[E   *numberOfLBnodes];
        dist.f[W   ] = &distributionArray[W   *numberOfLBnodes];
        dist.f[N   ] = &distributionArray[N   *numberOfLBnodes];
        dist.f[S   ] = &distributionArray[S   *numberOfLBnodes];
        dist.f[T   ] = &distributionArray[T   *numberOfLBnodes];
        dist.f[B   ] = &distributionArray[B   *numberOfLBnodes];
        dist.f[NE  ] = &distributionArray[NE  *numberOfLBnodes];
        dist.f[SW  ] = &distributionArray[SW  *numberOfLBnodes];
        dist.f[SE  ] = &distributionArray[SE  *numberOfLBnodes];
        dist.f[NW  ] = &distributionArray[NW  *numberOfLBnodes];
        dist.f[TE  ] = &distributionArray[TE  *numberOfLBnodes];
        dist.f[BW  ] = &distributionArray[BW  *numberOfLBnodes];
        dist.f[BE  ] = &distributionArray[BE  *numberOfLBnodes];
        dist.f[TW  ] = &distributionArray[TW  *numberOfLBnodes];
        dist.f[TN  ] = &distributionArray[TN  *numberOfLBnodes];
        dist.f[BS  ] = &distributionArray[BS  *numberOfLBnodes];
        dist.f[BN  ] = &distributionArray[BN  *numberOfLBnodes];
        dist.f[TS  ] = &distributionArray[TS  *numberOfLBnodes];
        dist.f[REST] = &distributionArray[REST*numberOfLBnodes];
        dist.f[TNE ] = &distributionArray[TNE *numberOfLBnodes];
        dist.f[TSW ] = &distributionArray[TSW *numberOfLBnodes];
        dist.f[TSE ] = &distributionArray[TSE *numberOfLBnodes];
        dist.f[TNW ] = &distributionArray[TNW *numberOfLBnodes];
        dist.f[BNE ] = &distributionArray[BNE *numberOfLBnodes];
        dist.f[BSW ] = &distributionArray[BSW *numberOfLBnodes];
        dist.f[BSE ] = &distributionArray[BSE *numberOfLBnodes];
        dist.f[BNW ] = &distributionArray[BNW *numberOfLBnodes];
    }
    else
    {
         dist.f[W   ] = &distributionArray[E   *numberOfLBnodes];
         dist.f[E   ] = &distributionArray[W   *numberOfLBnodes];
         dist.f[S   ] = &distributionArray[N   *numberOfLBnodes];
         dist.f[N   ] = &distributionArray[S   *numberOfLBnodes];
         dist.f[B   ] = &distributionArray[T   *numberOfLBnodes];
         dist.f[T   ] = &distributionArray[B   *numberOfLBnodes];
         dist.f[SW  ] = &distributionArray[NE  *numberOfLBnodes];
         dist.f[NE  ] = &distributionArray[SW  *numberOfLBnodes];
         dist.f[NW  ] = &distributionArray[SE  *numberOfLBnodes];
         dist.f[SE  ] = &distributionArray[NW  *numberOfLBnodes];
         dist.f[BW  ] = &distributionArray[TE  *numberOfLBnodes];
         dist.f[TE  ] = &distributionArray[BW  *numberOfLBnodes];
         dist.f[TW  ] = &distributionArray[BE  *numberOfLBnodes];
         dist.f[BE  ] = &distributionArray[TW  *numberOfLBnodes];
         dist.f[BS  ] = &distributionArray[TN  *numberOfLBnodes];
         dist.f[TN  ] = &distributionArray[BS  *numberOfLBnodes];
         dist.f[TS  ] = &distributionArray[BN  *numberOfLBnodes];
         dist.f[BN  ] = &distributionArray[TS  *numberOfLBnodes];
         dist.f[REST] = &distributionArray[REST*numberOfLBnodes];
         dist.f[TNE ] = &distributionArray[BSW *numberOfLBnodes];
         dist.f[TSW ] = &distributionArray[BNE *numberOfLBnodes];
         dist.f[TSE ] = &distributionArray[BNW *numberOfLBnodes];
         dist.f[TNW ] = &distributionArray[BSE *numberOfLBnodes];
         dist.f[BNE ] = &distributionArray[TSW *numberOfLBnodes];
         dist.f[BSW ] = &distributionArray[TNE *numberOfLBnodes];
         dist.f[BSE ] = &distributionArray[TNW *numberOfLBnodes];
         dist.f[BNW ] = &distributionArray[TSE *numberOfLBnodes];
    }
}

/**
*  Getting references to the 27 directions.
*  @params distributions 1D real* array containing all data (number of elements = 27 * matrix_size)
*  @params matrix_size number of discretizations nodes
*  @params isEvenTimestep: stored data dependent on timestep is based on the esoteric twist algorithm
*  @return a data struct containing the addresses to the 27 directions within the 1D distribution array
*/
__inline__ __device__ __host__ DistributionReferences27 getDistributionReferences27(real* distributions, unsigned int numberOfLBnodes, bool isEvenTimestep){
    DistributionReferences27 distribution_references;
    getPointersToDistributions(distribution_references, distributions, numberOfLBnodes, isEvenTimestep);
    return distribution_references;
}


/**
*  Holds the references to all directions and the concrete distributions for a single node.
*  After instantiation the distributions are read to the member "distribution" from "distribution_references".
*  After computation the data can be written back to "distribution_references".
*/
struct DistributionWrapper
{
    __device__ DistributionWrapper(
        real* distributions,
        unsigned int size_Mat,
        bool isEvenTimestep,
        uint k,
        uint* neighborX,
        uint* neighborY,
        uint* neighborZ);

    __device__ void read();

    __device__ void write();

    // origin distributions to read from and write to after computation
    DistributionReferences27 distribution_references;

    // distribution pass to kernel computation
    vf::lbm::Distribution27 distribution;

    const uint k;
    const uint kw;
    const uint ks;
    const uint kb;
    const uint ksw;
    const uint kbw;
    const uint kbs;
    const uint kbsw;
};

__inline__ __device__ unsigned int getNodeIndex()
{
    const unsigned x = threadIdx.x;
    const unsigned y = blockIdx.x;
    const unsigned z = blockIdx.y;

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    return nx * (ny * z + y) + x;
}

__device__ bool isValidFluidNode(uint nodeType);

}

#endif
