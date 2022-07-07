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
#include "LBM/D3Q27.h"
#include "lbm/constants/NumericConstants.h"

using namespace vf::lbm::constant;

__inline__ __device__ void getPointersToDistributions(Distributions27 &dist, real *distributionArray, const uint numberOfLBnodes, const bool isEvenTimestep)
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
        dist.f[dirREST] = &distributionArray[dirREST*numberOfLBnodes];
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
         dist.f[dirREST] = &distributionArray[dirREST*numberOfLBnodes];
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

__inline__ __device__ void getPointersToSubgridDistances(SubgridDistances27& subgridD, real* subgridDistances, const unsigned int numberOfSubgridIndices)
{
    subgridD.q[E   ] = &subgridDistances[E    *numberOfSubgridIndices];
    subgridD.q[W   ] = &subgridDistances[W    *numberOfSubgridIndices];
    subgridD.q[N   ] = &subgridDistances[N    *numberOfSubgridIndices];
    subgridD.q[S   ] = &subgridDistances[S    *numberOfSubgridIndices];
    subgridD.q[T   ] = &subgridDistances[T    *numberOfSubgridIndices];
    subgridD.q[B   ] = &subgridDistances[B    *numberOfSubgridIndices];
    subgridD.q[NE  ] = &subgridDistances[NE   *numberOfSubgridIndices];
    subgridD.q[SW  ] = &subgridDistances[SW   *numberOfSubgridIndices];
    subgridD.q[SE  ] = &subgridDistances[SE   *numberOfSubgridIndices];
    subgridD.q[NW  ] = &subgridDistances[NW   *numberOfSubgridIndices];
    subgridD.q[TE  ] = &subgridDistances[TE   *numberOfSubgridIndices];
    subgridD.q[BW  ] = &subgridDistances[BW   *numberOfSubgridIndices];
    subgridD.q[BE  ] = &subgridDistances[BE   *numberOfSubgridIndices];
    subgridD.q[TW  ] = &subgridDistances[TW   *numberOfSubgridIndices];
    subgridD.q[TN  ] = &subgridDistances[TN   *numberOfSubgridIndices];
    subgridD.q[BS  ] = &subgridDistances[BS   *numberOfSubgridIndices];
    subgridD.q[BN  ] = &subgridDistances[BN   *numberOfSubgridIndices];
    subgridD.q[TS  ] = &subgridDistances[TS   *numberOfSubgridIndices];
    subgridD.q[dirREST] = &subgridDistances[dirREST *numberOfSubgridIndices];
    subgridD.q[TNE ] = &subgridDistances[TNE  *numberOfSubgridIndices];
    subgridD.q[TSW ] = &subgridDistances[TSW  *numberOfSubgridIndices];
    subgridD.q[TSE ] = &subgridDistances[TSE  *numberOfSubgridIndices];
    subgridD.q[TNW ] = &subgridDistances[TNW  *numberOfSubgridIndices];
    subgridD.q[BNE ] = &subgridDistances[BNE  *numberOfSubgridIndices];
    subgridD.q[BSW ] = &subgridDistances[BSW  *numberOfSubgridIndices];
    subgridD.q[BSE ] = &subgridDistances[BSE  *numberOfSubgridIndices];
    subgridD.q[BNW ] = &subgridDistances[BNW  *numberOfSubgridIndices];
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
