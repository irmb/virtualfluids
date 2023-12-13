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

#include "Calculation/Calculation.h"

#include <array>
#include <cassert>
#include <optional>

#include <cuda.h>
#include <cuda_runtime.h>

#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

namespace vf::gpu
{

inline real getForceFactor(size_t level)
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
        dist.f[d0P0] = &distributionArray[d0P0 * numberOfLBnodes];
        dist.f[d0M0] = &distributionArray[d0M0 * numberOfLBnodes];
        dist.f[d00P] = &distributionArray[d00P * numberOfLBnodes];
        dist.f[d00M] = &distributionArray[d00M * numberOfLBnodes];
        dist.f[dPP0] = &distributionArray[dPP0 * numberOfLBnodes];
        dist.f[dMM0] = &distributionArray[dMM0 * numberOfLBnodes];
        dist.f[dPM0] = &distributionArray[dPM0 * numberOfLBnodes];
        dist.f[dMP0] = &distributionArray[dMP0 * numberOfLBnodes];
        dist.f[dP0P] = &distributionArray[dP0P * numberOfLBnodes];
        dist.f[dM0M] = &distributionArray[dM0M * numberOfLBnodes];
        dist.f[dP0M] = &distributionArray[dP0M * numberOfLBnodes];
        dist.f[dM0P] = &distributionArray[dM0P * numberOfLBnodes];
        dist.f[d0PP] = &distributionArray[d0PP * numberOfLBnodes];
        dist.f[d0MM] = &distributionArray[d0MM * numberOfLBnodes];
        dist.f[d0PM] = &distributionArray[d0PM * numberOfLBnodes];
        dist.f[d0MP] = &distributionArray[d0MP * numberOfLBnodes];
        dist.f[dPPP] = &distributionArray[dPPP * numberOfLBnodes];
        dist.f[dMMP] = &distributionArray[dMMP * numberOfLBnodes];
        dist.f[dPMP] = &distributionArray[dPMP * numberOfLBnodes];
        dist.f[dMPP] = &distributionArray[dMPP * numberOfLBnodes];
        dist.f[dPPM] = &distributionArray[dPPM * numberOfLBnodes];
        dist.f[dMMM] = &distributionArray[dMMM * numberOfLBnodes];
        dist.f[dPMM] = &distributionArray[dPMM * numberOfLBnodes];
        dist.f[dMPM] = &distributionArray[dMPM * numberOfLBnodes];
    } else {
        dist.f[dM00] = &distributionArray[dP00 * numberOfLBnodes];
        dist.f[dP00] = &distributionArray[dM00 * numberOfLBnodes];
        dist.f[d0M0] = &distributionArray[d0P0 * numberOfLBnodes];
        dist.f[d0P0] = &distributionArray[d0M0 * numberOfLBnodes];
        dist.f[d00M] = &distributionArray[d00P * numberOfLBnodes];
        dist.f[d00P] = &distributionArray[d00M * numberOfLBnodes];
        dist.f[dMM0] = &distributionArray[dPP0 * numberOfLBnodes];
        dist.f[dPP0] = &distributionArray[dMM0 * numberOfLBnodes];
        dist.f[dMP0] = &distributionArray[dPM0 * numberOfLBnodes];
        dist.f[dPM0] = &distributionArray[dMP0 * numberOfLBnodes];
        dist.f[dM0M] = &distributionArray[dP0P * numberOfLBnodes];
        dist.f[dP0P] = &distributionArray[dM0M * numberOfLBnodes];
        dist.f[dM0P] = &distributionArray[dP0M * numberOfLBnodes];
        dist.f[dP0M] = &distributionArray[dM0P * numberOfLBnodes];
        dist.f[d0MM] = &distributionArray[d0PP * numberOfLBnodes];
        dist.f[d0PP] = &distributionArray[d0MM * numberOfLBnodes];
        dist.f[d0MP] = &distributionArray[d0PM * numberOfLBnodes];
        dist.f[d0PM] = &distributionArray[d0MP * numberOfLBnodes];
        dist.f[d000] = &distributionArray[d000 * numberOfLBnodes];
        dist.f[dPPP] = &distributionArray[dMMM * numberOfLBnodes];
        dist.f[dMMP] = &distributionArray[dPPM * numberOfLBnodes];
        dist.f[dPMP] = &distributionArray[dMPM * numberOfLBnodes];
        dist.f[dMPP] = &distributionArray[dPMM * numberOfLBnodes];
        dist.f[dPPM] = &distributionArray[dMMP * numberOfLBnodes];
        dist.f[dMMM] = &distributionArray[dPPP * numberOfLBnodes];
        dist.f[dPMM] = &distributionArray[dMPP * numberOfLBnodes];
        dist.f[dMPM] = &distributionArray[dPMP * numberOfLBnodes];
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
    subgridD.q[d0P0] = &subgridDistances[d0P0 * numberOfSubgridIndices];
    subgridD.q[d0M0] = &subgridDistances[d0M0 * numberOfSubgridIndices];
    subgridD.q[d00P] = &subgridDistances[d00P * numberOfSubgridIndices];
    subgridD.q[d00M] = &subgridDistances[d00M * numberOfSubgridIndices];
    subgridD.q[dPP0] = &subgridDistances[dPP0 * numberOfSubgridIndices];
    subgridD.q[dMM0] = &subgridDistances[dMM0 * numberOfSubgridIndices];
    subgridD.q[dPM0] = &subgridDistances[dPM0 * numberOfSubgridIndices];
    subgridD.q[dMP0] = &subgridDistances[dMP0 * numberOfSubgridIndices];
    subgridD.q[dP0P] = &subgridDistances[dP0P * numberOfSubgridIndices];
    subgridD.q[dM0M] = &subgridDistances[dM0M * numberOfSubgridIndices];
    subgridD.q[dP0M] = &subgridDistances[dP0M * numberOfSubgridIndices];
    subgridD.q[dM0P] = &subgridDistances[dM0P * numberOfSubgridIndices];
    subgridD.q[d0PP] = &subgridDistances[d0PP * numberOfSubgridIndices];
    subgridD.q[d0MM] = &subgridDistances[d0MM * numberOfSubgridIndices];
    subgridD.q[d0PM] = &subgridDistances[d0PM * numberOfSubgridIndices];
    subgridD.q[d0MP] = &subgridDistances[d0MP * numberOfSubgridIndices];
    subgridD.q[d000] = &subgridDistances[d000 * numberOfSubgridIndices];
    subgridD.q[dPPP] = &subgridDistances[dPPP * numberOfSubgridIndices];
    subgridD.q[dMMP] = &subgridDistances[dMMP * numberOfSubgridIndices];
    subgridD.q[dPMP] = &subgridDistances[dPMP * numberOfSubgridIndices];
    subgridD.q[dMPP] = &subgridDistances[dMPP * numberOfSubgridIndices];
    subgridD.q[dPPM] = &subgridDistances[dPPM * numberOfSubgridIndices];
    subgridD.q[dMMM] = &subgridDistances[dMMM * numberOfSubgridIndices];
    subgridD.q[dPMM] = &subgridDistances[dPMM * numberOfSubgridIndices];
    subgridD.q[dMPM] = &subgridDistances[dMPM * numberOfSubgridIndices];
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
    __device__ ListIndices() {};
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
__device__ __inline__ void getPreCollisionDistribution(real* local, const Distributions27& global, const ListIndices& indices)
{
    local[d000] = (global.f[d000])[indices.k_000];
    local[dP00] = (global.f[dP00])[indices.k_000];
    local[dM00] = (global.f[dM00])[indices.k_M00];
    local[d0P0] = (global.f[d0P0])[indices.k_000];
    local[d0M0] = (global.f[d0M0])[indices.k_0M0];
    local[d00P] = (global.f[d00P])[indices.k_000];
    local[d00M] = (global.f[d00M])[indices.k_00M];
    local[dPP0] = (global.f[dPP0])[indices.k_000];
    local[dMM0] = (global.f[dMM0])[indices.k_MM0];
    local[dPM0] = (global.f[dPM0])[indices.k_0M0];
    local[dMP0] = (global.f[dMP0])[indices.k_M00];
    local[dP0P] = (global.f[dP0P])[indices.k_000];
    local[dM0M] = (global.f[dM0M])[indices.k_M0M];
    local[dP0M] = (global.f[dP0M])[indices.k_00M];
    local[dM0P] = (global.f[dM0P])[indices.k_M00];
    local[d0PP] = (global.f[d0PP])[indices.k_000];
    local[d0MM] = (global.f[d0MM])[indices.k_0MM];
    local[d0PM] = (global.f[d0PM])[indices.k_00M];
    local[d0MP] = (global.f[d0MP])[indices.k_0M0];
    local[dPPP] = (global.f[dPPP])[indices.k_000];
    local[dMPP] = (global.f[dMPP])[indices.k_M00];
    local[dPMP] = (global.f[dPMP])[indices.k_0M0];
    local[dMMP] = (global.f[dMMP])[indices.k_MM0];
    local[dPPM] = (global.f[dPPM])[indices.k_00M];
    local[dMPM] = (global.f[dMPM])[indices.k_M0M];
    local[dPMM] = (global.f[dPMM])[indices.k_0MM];
    local[dMMM] = (global.f[dMMM])[indices.k_MMM];
}

__device__ __inline__ void getPostCollisionDistribution(real* local, const Distributions27& global, const ListIndices& indices)
{
    local[d000] = (global.f[d000])[indices.k_000];
    local[dM00] = (global.f[dP00])[indices.k_000];
    local[dP00] = (global.f[dM00])[indices.k_M00];
    local[d0M0] = (global.f[d0P0])[indices.k_000];
    local[d0P0] = (global.f[d0M0])[indices.k_0M0];
    local[d00M] = (global.f[d00P])[indices.k_000];
    local[d00P] = (global.f[d00M])[indices.k_00M];
    local[dMM0] = (global.f[dPP0])[indices.k_000];
    local[dPP0] = (global.f[dMM0])[indices.k_MM0];
    local[dMP0] = (global.f[dPM0])[indices.k_0M0];
    local[dPM0] = (global.f[dMP0])[indices.k_M00];
    local[dM0M] = (global.f[dP0P])[indices.k_000];
    local[dP0P] = (global.f[dM0M])[indices.k_M0M];
    local[dM0P] = (global.f[dP0M])[indices.k_00M];
    local[dP0M] = (global.f[dM0P])[indices.k_M00];
    local[d0MM] = (global.f[d0PP])[indices.k_000];
    local[d0PP] = (global.f[d0MM])[indices.k_0MM];
    local[d0MP] = (global.f[d0PM])[indices.k_00M];
    local[d0PM] = (global.f[d0MP])[indices.k_0M0];
    local[dMMM] = (global.f[dPPP])[indices.k_000];
    local[dPMM] = (global.f[dMPP])[indices.k_M00];
    local[dMPM] = (global.f[dPMP])[indices.k_0M0];
    local[dPPM] = (global.f[dMMP])[indices.k_MM0];
    local[dMMP] = (global.f[dPPM])[indices.k_00M];
    local[dPMP] = (global.f[dMPM])[indices.k_M0M];
    local[dMPP] = (global.f[dPMM])[indices.k_0MM];
    local[dPPP] = (global.f[dMMM])[indices.k_MMM];
}

////////////////////////////////////////////////////////////////////////////////////
//! - Write distributions: style of reading and writing the distributions from/to
//! stored arrays dependent on timestep is based on the esoteric twist algorithm
//! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
//! DOI:10.3390/computation5020019 ]</b></a>
__inline__ __device__ void setPreCollisionDistribution(Distributions27& global, const ListIndices& indices, const real* local)
{
    (global.f[d000])[indices.k_000] = local[d000];
    (global.f[dP00])[indices.k_000] = local[dP00];
    (global.f[dM00])[indices.k_M00] = local[dM00];
    (global.f[d0P0])[indices.k_000] = local[d0P0];
    (global.f[d0M0])[indices.k_0M0] = local[d0M0];
    (global.f[d00P])[indices.k_000] = local[d00P];
    (global.f[d00M])[indices.k_00M] = local[d00M];
    (global.f[dPP0])[indices.k_000] = local[dPP0];
    (global.f[dMM0])[indices.k_MM0] = local[dMM0];
    (global.f[dPM0])[indices.k_0M0] = local[dPM0];
    (global.f[dMP0])[indices.k_M00] = local[dMP0];
    (global.f[dP0P])[indices.k_000] = local[dP0P];
    (global.f[dM0M])[indices.k_M0M] = local[dM0M];
    (global.f[dP0M])[indices.k_00M] = local[dP0M];
    (global.f[dM0P])[indices.k_M00] = local[dM0P];
    (global.f[d0PP])[indices.k_000] = local[d0PP];
    (global.f[d0MM])[indices.k_0MM] = local[d0MM];
    (global.f[d0PM])[indices.k_00M] = local[d0PM];
    (global.f[d0MP])[indices.k_0M0] = local[d0MP];
    (global.f[dPPP])[indices.k_000] = local[dPPP];
    (global.f[dMPP])[indices.k_M00] = local[dMPP];
    (global.f[dPMP])[indices.k_0M0] = local[dPMP];
    (global.f[dMMP])[indices.k_MM0] = local[dMMP];
    (global.f[dPPM])[indices.k_00M] = local[dPPM];
    (global.f[dMPM])[indices.k_M0M] = local[dMPM];
    (global.f[dPMM])[indices.k_0MM] = local[dPMM];
    (global.f[dMMM])[indices.k_MMM] = local[dMMM];
}

__inline__ __device__ void setPostCollisionDistribution(Distributions27& global, const ListIndices& indices, const real* local)
{
    (global.f[d000])[indices.k_000] = local[d000];
    (global.f[dP00])[indices.k_000] = local[dM00];
    (global.f[dM00])[indices.k_M00] = local[dP00];
    (global.f[d0P0])[indices.k_000] = local[d0M0];
    (global.f[d0M0])[indices.k_0M0] = local[d0P0];
    (global.f[d00P])[indices.k_000] = local[d00M];
    (global.f[d00M])[indices.k_00M] = local[d00P];
    (global.f[dPP0])[indices.k_000] = local[dMM0];
    (global.f[dMM0])[indices.k_MM0] = local[dPP0];
    (global.f[dPM0])[indices.k_0M0] = local[dMP0];
    (global.f[dMP0])[indices.k_M00] = local[dPM0];
    (global.f[dP0P])[indices.k_000] = local[dM0M];
    (global.f[dM0M])[indices.k_M0M] = local[dP0P];
    (global.f[dP0M])[indices.k_00M] = local[dM0P];
    (global.f[dM0P])[indices.k_M00] = local[dP0M];
    (global.f[d0PP])[indices.k_000] = local[d0MM];
    (global.f[d0MM])[indices.k_0MM] = local[d0PP];
    (global.f[d0PM])[indices.k_00M] = local[d0MP];
    (global.f[d0MP])[indices.k_0M0] = local[d0PM];
    (global.f[dPPP])[indices.k_000] = local[dMMM];
    (global.f[dMPP])[indices.k_M00] = local[dPMM];
    (global.f[dPMP])[indices.k_0M0] = local[dMPM];
    (global.f[dMMP])[indices.k_MM0] = local[dPPM];
    (global.f[dPPM])[indices.k_00M] = local[dMMP];
    (global.f[dMPM])[indices.k_M0M] = local[dPMP];
    (global.f[dPMM])[indices.k_0MM] = local[dMPP];
    (global.f[dMMM])[indices.k_MMM] = local[dPPP];
}

}

#endif
