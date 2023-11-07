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
//! \file KernelUtilities.h
//! \ingroup LBM/GPUHelperFunctions
//! \author Martin Schoenherr, Anna Wellmann, Soeren Peters
//=======================================================================================
#ifndef GPU_DISTRIBUTION_HELPER_H
#define GPU_DISTRIBUTION_HELPER_H

#include "LBM/LB.h"

#include <array>
#include <cassert>
#include <optional>

#include <cuda_runtime.h>

#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

namespace vf::gpu
{

inline real getForceFactor(int level)
{
    real factor = c1o1;
    for (size_t i = 1; i <= level; i++) {
        factor *= c2o1;
    }
    factor = 1 / factor;
    return factor;
}

__inline__ __device__ __host__ void getPointersToDistributions(Distributions27 &dist, real *distributionArray, const unsigned long long numberOfLBnodes, const bool isEvenTimestep)
{
    if (isEvenTimestep) {
        dist.f[d000] = &distributionArray[d000 * numberOfLBnodes];
        dist.f[dP00] = &distributionArray[dP00 * numberOfLBnodes];
        dist.f[dM00] = &distributionArray[dM00 * numberOfLBnodes];
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
    } else {
        dist.f[dM00] = &distributionArray[dP00 * numberOfLBnodes];
        dist.f[dP00] = &distributionArray[dM00 * numberOfLBnodes];
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
        dist.f[d000] = &distributionArray[d000 * numberOfLBnodes];
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
__inline__ __device__ __host__ DistributionReferences27 getDistributionReferences27(real* distributions, const unsigned long long numberOfLBnodes, const bool isEvenTimestep){
    DistributionReferences27 distribution_references;
    getPointersToDistributions(distribution_references, distributions, numberOfLBnodes, isEvenTimestep);
    return distribution_references;
}

__inline__ __device__ void getPointersToSubgridDistances(SubgridDistances27& subgridD, real* subgridDistances, const unsigned int numberOfSubgridIndices)
{
    subgridD.q[dP00] = &subgridDistances[dP00 * numberOfSubgridIndices];
    subgridD.q[dM00] = &subgridDistances[dM00 * numberOfSubgridIndices];
    subgridD.q[DIR_0P0] = &subgridDistances[DIR_0P0 * numberOfSubgridIndices];
    subgridD.q[DIR_0M0] = &subgridDistances[DIR_0M0 * numberOfSubgridIndices];
    subgridD.q[DIR_00P] = &subgridDistances[DIR_00P * numberOfSubgridIndices];
    subgridD.q[DIR_00M] = &subgridDistances[DIR_00M * numberOfSubgridIndices];
    subgridD.q[DIR_PP0] = &subgridDistances[DIR_PP0 * numberOfSubgridIndices];
    subgridD.q[DIR_MM0] = &subgridDistances[DIR_MM0 * numberOfSubgridIndices];
    subgridD.q[DIR_PM0] = &subgridDistances[DIR_PM0 * numberOfSubgridIndices];
    subgridD.q[DIR_MP0] = &subgridDistances[DIR_MP0 * numberOfSubgridIndices];
    subgridD.q[DIR_P0P] = &subgridDistances[DIR_P0P * numberOfSubgridIndices];
    subgridD.q[DIR_M0M] = &subgridDistances[DIR_M0M * numberOfSubgridIndices];
    subgridD.q[DIR_P0M] = &subgridDistances[DIR_P0M * numberOfSubgridIndices];
    subgridD.q[DIR_M0P] = &subgridDistances[DIR_M0P * numberOfSubgridIndices];
    subgridD.q[DIR_0PP] = &subgridDistances[DIR_0PP * numberOfSubgridIndices];
    subgridD.q[DIR_0MM] = &subgridDistances[DIR_0MM * numberOfSubgridIndices];
    subgridD.q[DIR_0PM] = &subgridDistances[DIR_0PM * numberOfSubgridIndices];
    subgridD.q[DIR_0MP] = &subgridDistances[DIR_0MP * numberOfSubgridIndices];
    subgridD.q[d000] = &subgridDistances[d000 * numberOfSubgridIndices];
    subgridD.q[DIR_PPP] = &subgridDistances[DIR_PPP * numberOfSubgridIndices];
    subgridD.q[DIR_MMP] = &subgridDistances[DIR_MMP * numberOfSubgridIndices];
    subgridD.q[DIR_PMP] = &subgridDistances[DIR_PMP * numberOfSubgridIndices];
    subgridD.q[DIR_MPP] = &subgridDistances[DIR_MPP * numberOfSubgridIndices];
    subgridD.q[DIR_PPM] = &subgridDistances[DIR_PPM * numberOfSubgridIndices];
    subgridD.q[DIR_MMM] = &subgridDistances[DIR_MMM * numberOfSubgridIndices];
    subgridD.q[DIR_PMM] = &subgridDistances[DIR_PMM * numberOfSubgridIndices];
    subgridD.q[DIR_MPM] = &subgridDistances[DIR_MPM * numberOfSubgridIndices];
}

__inline__ __device__ real getEquilibriumForBC(const real& drho, const real& velocity, const real& cu_sq, const real weight)
{
    return weight * (drho + c9o2 * velocity * velocity * (c1o1 + drho) - cu_sq);
}

__inline__ __device__ real getInterpolatedDistributionForVeloBC(const real& q, const real& f, const real& fInverse, const real& feq,
                                                                const real& omega, const real& velocity, const real weight)
{

    return (c1o1-q) / (c1o1+q) * (f - fInverse + (f + fInverse - c2o1 * feq * omega) / (c1o1 - omega)) * c1o2
           + (q * (f + fInverse) - c6o1 * weight * velocity) / (c1o1 + q);
}

__inline__ __device__ real getBounceBackDistributionForVeloBC(  const real& f,
                                                                const real& velocity, const real weight)
{

    return f - (c6o1 * weight * velocity);
}

__inline__ __device__ real getInterpolatedDistributionForNoSlipBC(const real& q, const real& f, const real& fInverse, const real& feq,
                                                                  const real& omega)
{

    return (c1o1-q) / (c1o1+q) * (f - fInverse + (f + fInverse - c2o1 * feq * omega) / (c1o1 - omega)) * c1o2
           + (q * (f + fInverse)) / (c1o1 + q);
}

__inline__ __device__ real getInterpolatedDistributionForNoSlipWithPressureBC(const real& q, const real& f, const real& fInverse, const real& feq, 
                                                                  const real& omega, const real& drho, const real weight)
{

    return (c1o1-q) / (c1o1+q) * (f - fInverse + (f + fInverse - c2o1 * feq * omega) / (c1o1 - omega)) * c1o2 
           + (q * (f + fInverse)) / (c1o1 + q) - weight * drho;
}


__inline__ __device__ real getInterpolatedDistributionForVeloWithPressureBC(const real& q, const real& f, const real& fInverse, const real& feq,
                                                                            const real& omega, const real& drho, const real& velocity, const real weight)
{

    return (c1o1-q) / (c1o1+q) * (f - fInverse + (f + fInverse - c2o1 * feq * omega) / (c1o1 - omega)) * c1o2
           + (q * (f + fInverse) - c6o1 * weight * velocity) / (c1o1 + q) - weight * drho;
}

__inline__ __device__ unsigned int getNodeIndex()
{
    const unsigned x = threadIdx.x;
    const unsigned y = blockIdx.x;
    const unsigned z = blockIdx.y;

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    return nx * (ny * z + y) + x;
}

__inline__ __device__ bool isValidFluidNode(uint nodeType)
{
    return (nodeType == GEO_FLUID || nodeType == GEO_PM_0 || nodeType == GEO_PM_1 || nodeType == GEO_PM_2);
}

struct ListIndices
{
    __device__ ListIndices() {} ;
    __device__ ListIndices(unsigned int index, const unsigned int* neighborX, const unsigned int* neighborY,
                           const unsigned int* neighborZ)
    {
        k_000 = index;
        k_M00 = neighborX[k_000];
        k_0M0 = neighborY[k_000];
        k_00M = neighborZ[k_000];
        k_MM0 = neighborY[k_M00];
        k_M0M = neighborZ[k_M00];
        k_0MM = neighborZ[k_0M0];
        k_MMM = neighborZ[k_MM0];
    }

    unsigned int k_000 { 0 };
    unsigned int k_M00 { 0 };
    unsigned int k_0M0 { 0 };
    unsigned int k_00M { 0 };
    unsigned int k_MM0 { 0 };
    unsigned int k_M0M { 0 };
    unsigned int k_0MM { 0 };
    unsigned int k_MMM { 0 };
};

////////////////////////////////////////////////////////////////////////////////////
//! - Read distributions: style of reading and writing the distributions from/to
//! stored arrays dependent on timestep is based on the esoteric twist algorithm
//! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
//! DOI:10.3390/computation5020019 ]</b></a>
__device__ __inline__ void getPreCollisionDistribution(real* destination, const Distributions27& source, const ListIndices& indices)
{
    destination[d000] = (source.f[d000])[indices.k_000];
    destination[dP00] = (source.f[dP00])[indices.k_000];
    destination[dM00] = (source.f[dM00])[indices.k_M00];
    destination[DIR_0P0] = (source.f[DIR_0P0])[indices.k_000];
    destination[DIR_0M0] = (source.f[DIR_0M0])[indices.k_0M0];
    destination[DIR_00P] = (source.f[DIR_00P])[indices.k_000];
    destination[DIR_00M] = (source.f[DIR_00M])[indices.k_00M];
    destination[DIR_PP0] = (source.f[DIR_PP0])[indices.k_000];
    destination[DIR_MM0] = (source.f[DIR_MM0])[indices.k_MM0];
    destination[DIR_PM0] = (source.f[DIR_PM0])[indices.k_0M0];
    destination[DIR_MP0] = (source.f[DIR_MP0])[indices.k_M00];
    destination[DIR_P0P] = (source.f[DIR_P0P])[indices.k_000];
    destination[DIR_M0M] = (source.f[DIR_M0M])[indices.k_M0M];
    destination[DIR_P0M] = (source.f[DIR_P0M])[indices.k_00M];
    destination[DIR_M0P] = (source.f[DIR_M0P])[indices.k_M00];
    destination[DIR_0PP] = (source.f[DIR_0PP])[indices.k_000];
    destination[DIR_0MM] = (source.f[DIR_0MM])[indices.k_0MM];
    destination[DIR_0PM] = (source.f[DIR_0PM])[indices.k_00M];
    destination[DIR_0MP] = (source.f[DIR_0MP])[indices.k_0M0];
    destination[DIR_PPP] = (source.f[DIR_PPP])[indices.k_000];
    destination[DIR_MPP] = (source.f[DIR_MPP])[indices.k_M00];
    destination[DIR_PMP] = (source.f[DIR_PMP])[indices.k_0M0];
    destination[DIR_MMP] = (source.f[DIR_MMP])[indices.k_MM0];
    destination[DIR_PPM] = (source.f[DIR_PPM])[indices.k_00M];
    destination[DIR_MPM] = (source.f[DIR_MPM])[indices.k_M0M];
    destination[DIR_PMM] = (source.f[DIR_PMM])[indices.k_0MM];
    destination[DIR_MMM] = (source.f[DIR_MMM])[indices.k_MMM];
}

__device__ __inline__ void getPostCollisionDistribution(real* destination, const Distributions27& source, const ListIndices& indices)
{
    destination[d000] = (source.f[d000])[indices.k_000];
    destination[dP00] = (source.f[dP00])[indices.k_000];
    destination[dM00] = (source.f[dM00])[indices.k_M00];
    destination[DIR_0P0] = (source.f[DIR_0P0])[indices.k_000];
    destination[DIR_0M0] = (source.f[DIR_0M0])[indices.k_0M0];
    destination[DIR_00P] = (source.f[DIR_00P])[indices.k_000];
    destination[DIR_00M] = (source.f[DIR_00M])[indices.k_00M];
    destination[DIR_PP0] = (source.f[DIR_PP0])[indices.k_000];
    destination[DIR_MM0] = (source.f[DIR_MM0])[indices.k_MM0];
    destination[DIR_PM0] = (source.f[DIR_PM0])[indices.k_0M0];
    destination[DIR_MP0] = (source.f[DIR_MP0])[indices.k_M00];
    destination[DIR_P0P] = (source.f[DIR_P0P])[indices.k_000];
    destination[DIR_M0M] = (source.f[DIR_M0M])[indices.k_M0M];
    destination[DIR_P0M] = (source.f[DIR_P0M])[indices.k_00M];
    destination[DIR_M0P] = (source.f[DIR_M0P])[indices.k_M00];
    destination[DIR_0PP] = (source.f[DIR_0PP])[indices.k_000];
    destination[DIR_0MM] = (source.f[DIR_0MM])[indices.k_0MM];
    destination[DIR_0PM] = (source.f[DIR_0PM])[indices.k_00M];
    destination[DIR_0MP] = (source.f[DIR_0MP])[indices.k_0M0];
    destination[DIR_PPP] = (source.f[DIR_PPP])[indices.k_000];
    destination[DIR_MPP] = (source.f[DIR_MPP])[indices.k_M00];
    destination[DIR_PMP] = (source.f[DIR_PMP])[indices.k_0M0];
    destination[DIR_MMP] = (source.f[DIR_MMP])[indices.k_MM0];
    destination[DIR_PPM] = (source.f[DIR_PPM])[indices.k_00M];
    destination[DIR_MPM] = (source.f[DIR_MPM])[indices.k_M0M];
    destination[DIR_PMM] = (source.f[DIR_PMM])[indices.k_0MM];
    destination[DIR_MMM] = (source.f[DIR_MMM])[indices.k_MMM];
}

////////////////////////////////////////////////////////////////////////////////////
//! - Write distributions: style of reading and writing the distributions from/to
//! stored arrays dependent on timestep is based on the esoteric twist algorithm
//! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
//! DOI:10.3390/computation5020019 ]</b></a>
__inline__ __device__ void setPreCollisionDistribution(Distributions27& destination, const ListIndices& indices, const real* source)
{
    (destination.f[d000])[indices.k_000] = source[d000];
    (destination.f[dP00])[indices.k_000] = source[dP00];
    (destination.f[dM00])[indices.k_M00] = source[dM00];
    (destination.f[DIR_0P0])[indices.k_000] = source[DIR_0P0];
    (destination.f[DIR_0M0])[indices.k_0M0] = source[DIR_0M0];
    (destination.f[DIR_00P])[indices.k_000] = source[DIR_00P];
    (destination.f[DIR_00M])[indices.k_00M] = source[DIR_00M];
    (destination.f[DIR_PP0])[indices.k_000] = source[DIR_PP0];
    (destination.f[DIR_MM0])[indices.k_MM0] = source[DIR_MM0];
    (destination.f[DIR_PM0])[indices.k_0M0] = source[DIR_PM0];
    (destination.f[DIR_MP0])[indices.k_M00] = source[DIR_MP0];
    (destination.f[DIR_P0P])[indices.k_000] = source[DIR_P0P];
    (destination.f[DIR_M0M])[indices.k_M0M] = source[DIR_M0M];
    (destination.f[DIR_P0M])[indices.k_00M] = source[DIR_P0M];
    (destination.f[DIR_M0P])[indices.k_M00] = source[DIR_M0P];
    (destination.f[DIR_0PP])[indices.k_000] = source[DIR_0PP];
    (destination.f[DIR_0MM])[indices.k_0MM] = source[DIR_0MM];
    (destination.f[DIR_0PM])[indices.k_00M] = source[DIR_0PM];
    (destination.f[DIR_0MP])[indices.k_0M0] = source[DIR_0MP];
    (destination.f[DIR_PPP])[indices.k_000] = source[DIR_PPP];
    (destination.f[DIR_MPP])[indices.k_M00] = source[DIR_MPP];
    (destination.f[DIR_PMP])[indices.k_0M0] = source[DIR_PMP];
    (destination.f[DIR_MMP])[indices.k_MM0] = source[DIR_MMP];
    (destination.f[DIR_PPM])[indices.k_00M] = source[DIR_PPM];
    (destination.f[DIR_MPM])[indices.k_M0M] = source[DIR_MPM];
    (destination.f[DIR_PMM])[indices.k_0MM] = source[DIR_PMM];
    (destination.f[DIR_MMM])[indices.k_MMM] = source[DIR_MMM];
}

__inline__ __device__ void setPostCollisionDistribution(Distributions27& destination, const ListIndices& indices, const real* source)
{
    (destination.f[d000])[indices.k_000] = source[d000];
    (destination.f[dP00])[indices.k_000] = source[dM00];
    (destination.f[dM00])[indices.k_M00] = source[dP00];
    (destination.f[DIR_0P0])[indices.k_000] = source[DIR_0M0];
    (destination.f[DIR_0M0])[indices.k_0M0] = source[DIR_0P0];
    (destination.f[DIR_00P])[indices.k_000] = source[DIR_00M];
    (destination.f[DIR_00M])[indices.k_00M] = source[DIR_00P];
    (destination.f[DIR_PP0])[indices.k_000] = source[DIR_MM0];
    (destination.f[DIR_MM0])[indices.k_MM0] = source[DIR_PP0];
    (destination.f[DIR_PM0])[indices.k_0M0] = source[DIR_MP0];
    (destination.f[DIR_MP0])[indices.k_M00] = source[DIR_PM0];
    (destination.f[DIR_P0P])[indices.k_000] = source[DIR_M0M];
    (destination.f[DIR_M0M])[indices.k_M0M] = source[DIR_P0P];
    (destination.f[DIR_P0M])[indices.k_00M] = source[DIR_M0P];
    (destination.f[DIR_M0P])[indices.k_M00] = source[DIR_P0M];
    (destination.f[DIR_0PP])[indices.k_000] = source[DIR_0MM];
    (destination.f[DIR_0MM])[indices.k_0MM] = source[DIR_0PP];
    (destination.f[DIR_0PM])[indices.k_00M] = source[DIR_0MP];
    (destination.f[DIR_0MP])[indices.k_0M0] = source[DIR_0PM];
    (destination.f[DIR_PPP])[indices.k_000] = source[DIR_MMM];
    (destination.f[DIR_MPP])[indices.k_M00] = source[DIR_PMM];
    (destination.f[DIR_PMP])[indices.k_0M0] = source[DIR_MPM];
    (destination.f[DIR_MMP])[indices.k_MM0] = source[DIR_PPM];
    (destination.f[DIR_PPM])[indices.k_00M] = source[DIR_MMP];
    (destination.f[DIR_MPM])[indices.k_M0M] = source[DIR_PMP];
    (destination.f[DIR_PMM])[indices.k_0MM] = source[DIR_MPP];
    (destination.f[DIR_MMM])[indices.k_MMM] = source[DIR_PPP];
}

}

#endif
