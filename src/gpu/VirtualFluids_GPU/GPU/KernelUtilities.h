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
//! \ingroup GPU
//! \author Martin Schoenherr, Anna Wellmann
//======================================================================================
#ifndef KERNELUTILS_H
#define KERNELUTILS_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

__inline__ __device__ void getPointersToDistributions(Distributions27 &dist, real *distributionArray, const uint numberOfLBnodes, const bool isEvenTimestep)
{
    if (isEvenTimestep)
    {
        dist.f[DIR_P00   ] = &distributionArray[DIR_P00   *numberOfLBnodes];
        dist.f[DIR_M00   ] = &distributionArray[DIR_M00   *numberOfLBnodes];
        dist.f[DIR_0P0   ] = &distributionArray[DIR_0P0   *numberOfLBnodes];
        dist.f[DIR_0M0   ] = &distributionArray[DIR_0M0   *numberOfLBnodes];
        dist.f[DIR_00P   ] = &distributionArray[DIR_00P   *numberOfLBnodes];
        dist.f[DIR_00M   ] = &distributionArray[DIR_00M   *numberOfLBnodes];
        dist.f[DIR_PP0  ] = &distributionArray[DIR_PP0  *numberOfLBnodes];
        dist.f[DIR_MM0  ] = &distributionArray[DIR_MM0  *numberOfLBnodes];
        dist.f[DIR_PM0  ] = &distributionArray[DIR_PM0  *numberOfLBnodes];
        dist.f[DIR_MP0  ] = &distributionArray[DIR_MP0  *numberOfLBnodes];
        dist.f[DIR_P0P  ] = &distributionArray[DIR_P0P  *numberOfLBnodes];
        dist.f[DIR_M0M  ] = &distributionArray[DIR_M0M  *numberOfLBnodes];
        dist.f[DIR_P0M  ] = &distributionArray[DIR_P0M  *numberOfLBnodes];
        dist.f[DIR_M0P  ] = &distributionArray[DIR_M0P  *numberOfLBnodes];
        dist.f[DIR_0PP  ] = &distributionArray[DIR_0PP  *numberOfLBnodes];
        dist.f[DIR_0MM  ] = &distributionArray[DIR_0MM  *numberOfLBnodes];
        dist.f[DIR_0PM  ] = &distributionArray[DIR_0PM  *numberOfLBnodes];
        dist.f[DIR_0MP  ] = &distributionArray[DIR_0MP  *numberOfLBnodes];
        dist.f[DIR_000] = &distributionArray[DIR_000*numberOfLBnodes];
        dist.f[DIR_PPP ] = &distributionArray[DIR_PPP *numberOfLBnodes];
        dist.f[DIR_MMP ] = &distributionArray[DIR_MMP *numberOfLBnodes];
        dist.f[DIR_PMP ] = &distributionArray[DIR_PMP *numberOfLBnodes];
        dist.f[DIR_MPP ] = &distributionArray[DIR_MPP *numberOfLBnodes];
        dist.f[DIR_PPM ] = &distributionArray[DIR_PPM *numberOfLBnodes];
        dist.f[DIR_MMM ] = &distributionArray[DIR_MMM *numberOfLBnodes];
        dist.f[DIR_PMM ] = &distributionArray[DIR_PMM *numberOfLBnodes];
        dist.f[DIR_MPM ] = &distributionArray[DIR_MPM *numberOfLBnodes];
    }
    else
    {
         dist.f[DIR_M00   ] = &distributionArray[DIR_P00   *numberOfLBnodes];
         dist.f[DIR_P00   ] = &distributionArray[DIR_M00   *numberOfLBnodes];
         dist.f[DIR_0M0   ] = &distributionArray[DIR_0P0   *numberOfLBnodes];
         dist.f[DIR_0P0   ] = &distributionArray[DIR_0M0   *numberOfLBnodes];
         dist.f[DIR_00M   ] = &distributionArray[DIR_00P   *numberOfLBnodes];
         dist.f[DIR_00P   ] = &distributionArray[DIR_00M   *numberOfLBnodes];
         dist.f[DIR_MM0  ] = &distributionArray[DIR_PP0  *numberOfLBnodes];
         dist.f[DIR_PP0  ] = &distributionArray[DIR_MM0  *numberOfLBnodes];
         dist.f[DIR_MP0  ] = &distributionArray[DIR_PM0  *numberOfLBnodes];
         dist.f[DIR_PM0  ] = &distributionArray[DIR_MP0  *numberOfLBnodes];
         dist.f[DIR_M0M  ] = &distributionArray[DIR_P0P  *numberOfLBnodes];
         dist.f[DIR_P0P  ] = &distributionArray[DIR_M0M  *numberOfLBnodes];
         dist.f[DIR_M0P  ] = &distributionArray[DIR_P0M  *numberOfLBnodes];
         dist.f[DIR_P0M  ] = &distributionArray[DIR_M0P  *numberOfLBnodes];
         dist.f[DIR_0MM  ] = &distributionArray[DIR_0PP  *numberOfLBnodes];
         dist.f[DIR_0PP  ] = &distributionArray[DIR_0MM  *numberOfLBnodes];
         dist.f[DIR_0MP  ] = &distributionArray[DIR_0PM  *numberOfLBnodes];
         dist.f[DIR_0PM  ] = &distributionArray[DIR_0MP  *numberOfLBnodes];
         dist.f[DIR_000] = &distributionArray[DIR_000*numberOfLBnodes];
         dist.f[DIR_PPP ] = &distributionArray[DIR_MMM *numberOfLBnodes];
         dist.f[DIR_MMP ] = &distributionArray[DIR_PPM *numberOfLBnodes];
         dist.f[DIR_PMP ] = &distributionArray[DIR_MPM *numberOfLBnodes];
         dist.f[DIR_MPP ] = &distributionArray[DIR_PMM *numberOfLBnodes];
         dist.f[DIR_PPM ] = &distributionArray[DIR_MMP *numberOfLBnodes];
         dist.f[DIR_MMM ] = &distributionArray[DIR_PPP *numberOfLBnodes];
         dist.f[DIR_PMM ] = &distributionArray[DIR_MPP *numberOfLBnodes];
         dist.f[DIR_MPM ] = &distributionArray[DIR_PMP *numberOfLBnodes];
    }
}

__inline__ __device__ void getPointersToSubgridDistances(SubgridDistances27& subgridD, real* subgridDistances, const unsigned int numberOfSubgridIndices)
{
    subgridD.q[DIR_P00   ] = &subgridDistances[DIR_P00    *numberOfSubgridIndices];
    subgridD.q[DIR_M00   ] = &subgridDistances[DIR_M00    *numberOfSubgridIndices];
    subgridD.q[DIR_0P0   ] = &subgridDistances[DIR_0P0    *numberOfSubgridIndices];
    subgridD.q[DIR_0M0   ] = &subgridDistances[DIR_0M0    *numberOfSubgridIndices];
    subgridD.q[DIR_00P   ] = &subgridDistances[DIR_00P    *numberOfSubgridIndices];
    subgridD.q[DIR_00M   ] = &subgridDistances[DIR_00M    *numberOfSubgridIndices];
    subgridD.q[DIR_PP0  ] = &subgridDistances[DIR_PP0   *numberOfSubgridIndices];
    subgridD.q[DIR_MM0  ] = &subgridDistances[DIR_MM0   *numberOfSubgridIndices];
    subgridD.q[DIR_PM0  ] = &subgridDistances[DIR_PM0   *numberOfSubgridIndices];
    subgridD.q[DIR_MP0  ] = &subgridDistances[DIR_MP0   *numberOfSubgridIndices];
    subgridD.q[DIR_P0P  ] = &subgridDistances[DIR_P0P   *numberOfSubgridIndices];
    subgridD.q[DIR_M0M  ] = &subgridDistances[DIR_M0M   *numberOfSubgridIndices];
    subgridD.q[DIR_P0M  ] = &subgridDistances[DIR_P0M   *numberOfSubgridIndices];
    subgridD.q[DIR_M0P  ] = &subgridDistances[DIR_M0P   *numberOfSubgridIndices];
    subgridD.q[DIR_0PP  ] = &subgridDistances[DIR_0PP   *numberOfSubgridIndices];
    subgridD.q[DIR_0MM  ] = &subgridDistances[DIR_0MM   *numberOfSubgridIndices];
    subgridD.q[DIR_0PM  ] = &subgridDistances[DIR_0PM   *numberOfSubgridIndices];
    subgridD.q[DIR_0MP  ] = &subgridDistances[DIR_0MP   *numberOfSubgridIndices];
    subgridD.q[DIR_000] = &subgridDistances[DIR_000 *numberOfSubgridIndices];
    subgridD.q[DIR_PPP ] = &subgridDistances[DIR_PPP  *numberOfSubgridIndices];
    subgridD.q[DIR_MMP ] = &subgridDistances[DIR_MMP  *numberOfSubgridIndices];
    subgridD.q[DIR_PMP ] = &subgridDistances[DIR_PMP  *numberOfSubgridIndices];
    subgridD.q[DIR_MPP ] = &subgridDistances[DIR_MPP  *numberOfSubgridIndices];
    subgridD.q[DIR_PPM ] = &subgridDistances[DIR_PPM  *numberOfSubgridIndices];
    subgridD.q[DIR_MMM ] = &subgridDistances[DIR_MMM  *numberOfSubgridIndices];
    subgridD.q[DIR_PMM ] = &subgridDistances[DIR_PMM  *numberOfSubgridIndices];
    subgridD.q[DIR_MPM ] = &subgridDistances[DIR_MPM  *numberOfSubgridIndices];
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

__inline__ __device__ real getInterpolatedDistributionForNoSlipBC(const real& q, const real& f, const real& fInverse, const real& feq, 
                                                                  const real& omega)
{

    return (c1o1-q) / (c1o1+q) * (f - fInverse + (f + fInverse - c2o1 * feq * omega) / (c1o1 - omega)) * c1o2 
           + (q * (f + fInverse)) / (c1o1 + q);
}


__inline__ __device__ real getInterpolatedDistributionForVeloWithPressureBC(const real& q, const real& f, const real& fInverse, const real& feq, 
                                                                            const real& omega, const real& drho, const real& velocity, const real weight)
{

    return (c1o1-q) / (c1o1+q) * (f - fInverse + (f + fInverse - c2o1 * feq * omega) / (c1o1 - omega)) * c1o2 
           + (q * (f + fInverse) - c6o1 * weight * velocity) / (c1o1 + q) - weight * drho;
}



#endif
