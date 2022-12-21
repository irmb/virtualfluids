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
        dist.f[DIR_000] = &distributionArray[(unsigned long long)(DIR_000) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_P00] = &distributionArray[(unsigned long long)(DIR_P00) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_M00] = &distributionArray[(unsigned long long)(DIR_M00) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_0P0] = &distributionArray[(unsigned long long)(DIR_0P0) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_0M0] = &distributionArray[(unsigned long long)(DIR_0M0) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_00P] = &distributionArray[(unsigned long long)(DIR_00P) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_00M] = &distributionArray[(unsigned long long)(DIR_00M) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_PP0] = &distributionArray[(unsigned long long)(DIR_PP0) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_MM0] = &distributionArray[(unsigned long long)(DIR_MM0) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_PM0] = &distributionArray[(unsigned long long)(DIR_PM0) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_MP0] = &distributionArray[(unsigned long long)(DIR_MP0) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_P0P] = &distributionArray[(unsigned long long)(DIR_P0P) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_M0M] = &distributionArray[(unsigned long long)(DIR_M0M) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_P0M] = &distributionArray[(unsigned long long)(DIR_P0M) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_M0P] = &distributionArray[(unsigned long long)(DIR_M0P) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_0PP] = &distributionArray[(unsigned long long)(DIR_0PP) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_0MM] = &distributionArray[(unsigned long long)(DIR_0MM) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_0PM] = &distributionArray[(unsigned long long)(DIR_0PM) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_0MP] = &distributionArray[(unsigned long long)(DIR_0MP) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_PPP] = &distributionArray[(unsigned long long)(DIR_PPP) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_MMP] = &distributionArray[(unsigned long long)(DIR_MMP) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_PMP] = &distributionArray[(unsigned long long)(DIR_PMP) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_MPP] = &distributionArray[(unsigned long long)(DIR_MPP) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_PPM] = &distributionArray[(unsigned long long)(DIR_PPM) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_MMM] = &distributionArray[(unsigned long long)(DIR_MMM) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_PMM] = &distributionArray[(unsigned long long)(DIR_PMM) * (unsigned long long)(numberOfLBnodes)];
        dist.f[DIR_MPM] = &distributionArray[(unsigned long long)(DIR_MPM) * (unsigned long long)(numberOfLBnodes)];
    }
    else
    {
         dist.f[DIR_M00] = &distributionArray[(unsigned long long)(DIR_P00) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_P00] = &distributionArray[(unsigned long long)(DIR_M00) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_0M0] = &distributionArray[(unsigned long long)(DIR_0P0) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_0P0] = &distributionArray[(unsigned long long)(DIR_0M0) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_00M] = &distributionArray[(unsigned long long)(DIR_00P) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_00P] = &distributionArray[(unsigned long long)(DIR_00M) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_MM0] = &distributionArray[(unsigned long long)(DIR_PP0) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_PP0] = &distributionArray[(unsigned long long)(DIR_MM0) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_MP0] = &distributionArray[(unsigned long long)(DIR_PM0) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_PM0] = &distributionArray[(unsigned long long)(DIR_MP0) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_M0M] = &distributionArray[(unsigned long long)(DIR_P0P) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_P0P] = &distributionArray[(unsigned long long)(DIR_M0M) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_M0P] = &distributionArray[(unsigned long long)(DIR_P0M) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_P0M] = &distributionArray[(unsigned long long)(DIR_M0P) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_0MM] = &distributionArray[(unsigned long long)(DIR_0PP) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_0PP] = &distributionArray[(unsigned long long)(DIR_0MM) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_0MP] = &distributionArray[(unsigned long long)(DIR_0PM) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_0PM] = &distributionArray[(unsigned long long)(DIR_0MP) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_000] = &distributionArray[(unsigned long long)(DIR_000) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_PPP] = &distributionArray[(unsigned long long)(DIR_MMM) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_MMP] = &distributionArray[(unsigned long long)(DIR_PPM) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_PMP] = &distributionArray[(unsigned long long)(DIR_MPM) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_MPP] = &distributionArray[(unsigned long long)(DIR_PMM) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_PPM] = &distributionArray[(unsigned long long)(DIR_MMP) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_MMM] = &distributionArray[(unsigned long long)(DIR_PPP) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_PMM] = &distributionArray[(unsigned long long)(DIR_MPP) * (unsigned long long)(numberOfLBnodes)];
         dist.f[DIR_MPM] = &distributionArray[(unsigned long long)(DIR_PMP) * (unsigned long long)(numberOfLBnodes)];
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
