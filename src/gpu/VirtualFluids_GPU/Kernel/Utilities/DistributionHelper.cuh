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
        dist.f[DIR_000] = &distributionArray[DIR_000 * numberOfLBnodes];
        dist.f[DIR_P00] = &distributionArray[DIR_P00 * numberOfLBnodes];
        dist.f[DIR_M00] = &distributionArray[DIR_M00 * numberOfLBnodes];
        dist.f[DIR_0P0] = &distributionArray[DIR_0P0 * numberOfLBnodes];
        dist.f[DIR_0M0] = &distributionArray[DIR_0M0 * numberOfLBnodes];
        dist.f[DIR_00P] = &distributionArray[DIR_00P * numberOfLBnodes];
        dist.f[DIR_00M] = &distributionArray[DIR_00M * numberOfLBnodes];
        dist.f[DIR_PP0] = &distributionArray[DIR_PP0 * numberOfLBnodes];
        dist.f[DIR_MM0] = &distributionArray[DIR_MM0 * numberOfLBnodes];
        dist.f[DIR_PM0] = &distributionArray[DIR_PM0 * numberOfLBnodes];
        dist.f[DIR_MP0] = &distributionArray[DIR_MP0 * numberOfLBnodes];
        dist.f[DIR_P0P] = &distributionArray[DIR_P0P * numberOfLBnodes];
        dist.f[DIR_M0M] = &distributionArray[DIR_M0M * numberOfLBnodes];
        dist.f[DIR_P0M] = &distributionArray[DIR_P0M * numberOfLBnodes];
        dist.f[DIR_M0P] = &distributionArray[DIR_M0P * numberOfLBnodes];
        dist.f[DIR_0PP] = &distributionArray[DIR_0PP * numberOfLBnodes];
        dist.f[DIR_0MM] = &distributionArray[DIR_0MM * numberOfLBnodes];
        dist.f[DIR_0PM] = &distributionArray[DIR_0PM * numberOfLBnodes];
        dist.f[DIR_0MP] = &distributionArray[DIR_0MP * numberOfLBnodes];
        dist.f[DIR_PPP] = &distributionArray[DIR_PPP * numberOfLBnodes];
        dist.f[DIR_MMP] = &distributionArray[DIR_MMP * numberOfLBnodes];
        dist.f[DIR_PMP] = &distributionArray[DIR_PMP * numberOfLBnodes];
        dist.f[DIR_MPP] = &distributionArray[DIR_MPP * numberOfLBnodes];
        dist.f[DIR_PPM] = &distributionArray[DIR_PPM * numberOfLBnodes];
        dist.f[DIR_MMM] = &distributionArray[DIR_MMM * numberOfLBnodes];
        dist.f[DIR_PMM] = &distributionArray[DIR_PMM * numberOfLBnodes];
        dist.f[DIR_MPM] = &distributionArray[DIR_MPM * numberOfLBnodes];
    }
    else
    {
         dist.f[DIR_M00] = &distributionArray[DIR_P00 * numberOfLBnodes];
         dist.f[DIR_P00] = &distributionArray[DIR_M00 * numberOfLBnodes];
         dist.f[DIR_0M0] = &distributionArray[DIR_0P0 * numberOfLBnodes];
         dist.f[DIR_0P0] = &distributionArray[DIR_0M0 * numberOfLBnodes];
         dist.f[DIR_00M] = &distributionArray[DIR_00P * numberOfLBnodes];
         dist.f[DIR_00P] = &distributionArray[DIR_00M * numberOfLBnodes];
         dist.f[DIR_MM0] = &distributionArray[DIR_PP0 * numberOfLBnodes];
         dist.f[DIR_PP0] = &distributionArray[DIR_MM0 * numberOfLBnodes];
         dist.f[DIR_MP0] = &distributionArray[DIR_PM0 * numberOfLBnodes];
         dist.f[DIR_PM0] = &distributionArray[DIR_MP0 * numberOfLBnodes];
         dist.f[DIR_M0M] = &distributionArray[DIR_P0P * numberOfLBnodes];
         dist.f[DIR_P0P] = &distributionArray[DIR_M0M * numberOfLBnodes];
         dist.f[DIR_M0P] = &distributionArray[DIR_P0M * numberOfLBnodes];
         dist.f[DIR_P0M] = &distributionArray[DIR_M0P * numberOfLBnodes];
         dist.f[DIR_0MM] = &distributionArray[DIR_0PP * numberOfLBnodes];
         dist.f[DIR_0PP] = &distributionArray[DIR_0MM * numberOfLBnodes];
         dist.f[DIR_0MP] = &distributionArray[DIR_0PM * numberOfLBnodes];
         dist.f[DIR_0PM] = &distributionArray[DIR_0MP * numberOfLBnodes];
         dist.f[DIR_000] = &distributionArray[DIR_000 * numberOfLBnodes];
         dist.f[DIR_PPP] = &distributionArray[DIR_MMM * numberOfLBnodes];
         dist.f[DIR_MMP] = &distributionArray[DIR_PPM * numberOfLBnodes];
         dist.f[DIR_PMP] = &distributionArray[DIR_MPM * numberOfLBnodes];
         dist.f[DIR_MPP] = &distributionArray[DIR_PMM * numberOfLBnodes];
         dist.f[DIR_PPM] = &distributionArray[DIR_MMP * numberOfLBnodes];
         dist.f[DIR_MMM] = &distributionArray[DIR_PPP * numberOfLBnodes];
         dist.f[DIR_PMM] = &distributionArray[DIR_MPP * numberOfLBnodes];
         dist.f[DIR_MPM] = &distributionArray[DIR_PMP * numberOfLBnodes];
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
