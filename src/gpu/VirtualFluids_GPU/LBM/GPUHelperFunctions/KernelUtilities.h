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
#ifndef KERNEL_UTILITIES_H
#define KERNEL_UTILITIES_H

#include "LBM/LB.h"

#include <lbm/constants/D3Q27.h>
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

namespace vf::gpu
{

__inline__ __device__ __host__ void getPointersToDistributions(Distributions27 &dist, real *distributionArray, const unsigned long long numberOfLBnodes, const bool isEvenTimestep)
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

__inline__ __device__ void getPointersToSubgridDistances(SubgridDistances27& subgridD, real* subgridDistances, const unsigned int numberOfSubgridIndices)
{
    subgridD.q[DIR_P00] = &subgridDistances[DIR_P00 * numberOfSubgridIndices];
    subgridD.q[DIR_M00] = &subgridDistances[DIR_M00 * numberOfSubgridIndices];
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
    subgridD.q[DIR_000] = &subgridDistances[DIR_000 * numberOfSubgridIndices];
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
    __device__ ListIndices(unsigned int k, unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ)
    {
        k_000 = k;
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
__device__ __inline__ void read(real *f, const Distributions27 &dist, const ListIndices &indices)
{
    f[DIR_000] = (dist.f[DIR_000])[indices.k_000];
    f[DIR_P00] = (dist.f[DIR_P00])[indices.k_000];
    f[DIR_M00] = (dist.f[DIR_M00])[indices.k_M00];
    f[DIR_0P0] = (dist.f[DIR_0P0])[indices.k_000];
    f[DIR_0M0] = (dist.f[DIR_0M0])[indices.k_0M0];
    f[DIR_00P] = (dist.f[DIR_00P])[indices.k_000];
    f[DIR_00M] = (dist.f[DIR_00M])[indices.k_00M];
    f[DIR_PP0] = (dist.f[DIR_PP0])[indices.k_000];
    f[DIR_MM0] = (dist.f[DIR_MM0])[indices.k_MM0];
    f[DIR_PM0] = (dist.f[DIR_PM0])[indices.k_0M0];
    f[DIR_MP0] = (dist.f[DIR_MP0])[indices.k_M00];
    f[DIR_P0P] = (dist.f[DIR_P0P])[indices.k_000];
    f[DIR_M0M] = (dist.f[DIR_M0M])[indices.k_M0M];
    f[DIR_P0M] = (dist.f[DIR_P0M])[indices.k_00M];
    f[DIR_M0P] = (dist.f[DIR_M0P])[indices.k_M00];
    f[DIR_0PP] = (dist.f[DIR_0PP])[indices.k_000];
    f[DIR_0MM] = (dist.f[DIR_0MM])[indices.k_0MM];
    f[DIR_0PM] = (dist.f[DIR_0PM])[indices.k_00M];
    f[DIR_0MP] = (dist.f[DIR_0MP])[indices.k_0M0];
    f[DIR_PPP] = (dist.f[DIR_PPP])[indices.k_000];
    f[DIR_MPP] = (dist.f[DIR_MPP])[indices.k_M00];
    f[DIR_PMP] = (dist.f[DIR_PMP])[indices.k_0M0];
    f[DIR_MMP] = (dist.f[DIR_MMP])[indices.k_MM0];
    f[DIR_PPM] = (dist.f[DIR_PPM])[indices.k_00M];
    f[DIR_MPM] = (dist.f[DIR_MPM])[indices.k_M0M];
    f[DIR_PMM] = (dist.f[DIR_PMM])[indices.k_0MM];
    f[DIR_MMM] = (dist.f[DIR_MMM])[indices.k_MMM];
}

__device__ __inline__ void readInverse(real *f, const Distributions27 &dist, const ListIndices &indices)
{
    //TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/101
    assert(false) || !fprintf(stderr, "Not implemented yet.\n")); 
}


////////////////////////////////////////////////////////////////////////////////////
//! - Write distributions: style of reading and writing the distributions from/to
//! stored arrays dependent on timestep is based on the esoteric twist algorithm
//! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
//! DOI:10.3390/computation5020019 ]</b></a>
__inline__ __device__ void write(Distributions27 &destination, const ListIndices &indices, const real* f)
{
    (destination.f[DIR_000])[indices.k_000] = f[DIR_000];
    (destination.f[DIR_P00])[indices.k_000] = f[DIR_P00];
    (destination.f[DIR_M00])[indices.k_M00] = f[DIR_M00];
    (destination.f[DIR_0P0])[indices.k_000] = f[DIR_0P0];
    (destination.f[DIR_0M0])[indices.k_0M0] = f[DIR_0M0];
    (destination.f[DIR_00P])[indices.k_000] = f[DIR_00P];
    (destination.f[DIR_00M])[indices.k_00M] = f[DIR_00M];
    (destination.f[DIR_PP0])[indices.k_000] = f[DIR_PP0];
    (destination.f[DIR_MM0])[indices.k_MM0] = f[DIR_MM0];
    (destination.f[DIR_PM0])[indices.k_0M0] = f[DIR_PM0];
    (destination.f[DIR_MP0])[indices.k_M00] = f[DIR_MP0];
    (destination.f[DIR_P0P])[indices.k_000] = f[DIR_P0P];
    (destination.f[DIR_M0M])[indices.k_M0M] = f[DIR_M0M];
    (destination.f[DIR_P0M])[indices.k_00M] = f[DIR_P0M];
    (destination.f[DIR_M0P])[indices.k_M00] = f[DIR_M0P];
    (destination.f[DIR_0PP])[indices.k_000] = f[DIR_0PP];
    (destination.f[DIR_0MM])[indices.k_0MM] = f[DIR_0MM];
    (destination.f[DIR_0PM])[indices.k_00M] = f[DIR_0PM];
    (destination.f[DIR_0MP])[indices.k_0M0] = f[DIR_0MP];
    (destination.f[DIR_PPP])[indices.k_000] = f[DIR_PPP];
    (destination.f[DIR_MPP])[indices.k_M00] = f[DIR_MPP];
    (destination.f[DIR_PMP])[indices.k_0M0] = f[DIR_PMP];
    (destination.f[DIR_MMP])[indices.k_MM0] = f[DIR_MMP];
    (destination.f[DIR_PPM])[indices.k_00M] = f[DIR_PPM];
    (destination.f[DIR_MPM])[indices.k_M0M] = f[DIR_MPM];
    (destination.f[DIR_PMM])[indices.k_0MM] = f[DIR_PMM];
    (destination.f[DIR_MMM])[indices.k_MMM] = f[DIR_MMM];
}

__inline__ __device__ void writeInverse(Distributions27 &destination, const ListIndices &indices, const real* f)
{
    //TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/101
    assert(false) || !fprintf(stderr, "Not implemented yet.\n")); 
}




__inline__ __device__ real calculateOmega(const real omega_old, real turbulenceViscosity)
{
    return omega_old / (c1o1 + c3o1 * omega_old * turbulenceViscosity);
}

}

#endif
