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
        dist.f[dirE   ] = &distributionArray[dirE   *numberOfLBnodes];
        dist.f[dirW   ] = &distributionArray[dirW   *numberOfLBnodes];
        dist.f[dirN   ] = &distributionArray[dirN   *numberOfLBnodes];
        dist.f[dirS   ] = &distributionArray[dirS   *numberOfLBnodes];
        dist.f[dirT   ] = &distributionArray[dirT   *numberOfLBnodes];
        dist.f[dirB   ] = &distributionArray[dirB   *numberOfLBnodes];
        dist.f[dirNE  ] = &distributionArray[dirNE  *numberOfLBnodes];
        dist.f[dirSW  ] = &distributionArray[dirSW  *numberOfLBnodes];
        dist.f[dirSE  ] = &distributionArray[dirSE  *numberOfLBnodes];
        dist.f[dirNW  ] = &distributionArray[dirNW  *numberOfLBnodes];
        dist.f[dirTE  ] = &distributionArray[dirTE  *numberOfLBnodes];
        dist.f[dirBW  ] = &distributionArray[dirBW  *numberOfLBnodes];
        dist.f[dirBE  ] = &distributionArray[dirBE  *numberOfLBnodes];
        dist.f[dirTW  ] = &distributionArray[dirTW  *numberOfLBnodes];
        dist.f[dirTN  ] = &distributionArray[dirTN  *numberOfLBnodes];
        dist.f[dirBS  ] = &distributionArray[dirBS  *numberOfLBnodes];
        dist.f[dirBN  ] = &distributionArray[dirBN  *numberOfLBnodes];
        dist.f[dirTS  ] = &distributionArray[dirTS  *numberOfLBnodes];
        dist.f[dirZERO] = &distributionArray[dirZERO*numberOfLBnodes];
        dist.f[dirTNE ] = &distributionArray[dirTNE *numberOfLBnodes];
        dist.f[dirTSW ] = &distributionArray[dirTSW *numberOfLBnodes];
        dist.f[dirTSE ] = &distributionArray[dirTSE *numberOfLBnodes];
        dist.f[dirTNW ] = &distributionArray[dirTNW *numberOfLBnodes];
        dist.f[dirBNE ] = &distributionArray[dirBNE *numberOfLBnodes];
        dist.f[dirBSW ] = &distributionArray[dirBSW *numberOfLBnodes];
        dist.f[dirBSE ] = &distributionArray[dirBSE *numberOfLBnodes];
        dist.f[dirBNW ] = &distributionArray[dirBNW *numberOfLBnodes];
    }
    else
    {
         dist.f[dirW   ] = &distributionArray[dirE   *numberOfLBnodes];
         dist.f[dirE   ] = &distributionArray[dirW   *numberOfLBnodes];
         dist.f[dirS   ] = &distributionArray[dirN   *numberOfLBnodes];
         dist.f[dirN   ] = &distributionArray[dirS   *numberOfLBnodes];
         dist.f[dirB   ] = &distributionArray[dirT   *numberOfLBnodes];
         dist.f[dirT   ] = &distributionArray[dirB   *numberOfLBnodes];
         dist.f[dirSW  ] = &distributionArray[dirNE  *numberOfLBnodes];
         dist.f[dirNE  ] = &distributionArray[dirSW  *numberOfLBnodes];
         dist.f[dirNW  ] = &distributionArray[dirSE  *numberOfLBnodes];
         dist.f[dirSE  ] = &distributionArray[dirNW  *numberOfLBnodes];
         dist.f[dirBW  ] = &distributionArray[dirTE  *numberOfLBnodes];
         dist.f[dirTE  ] = &distributionArray[dirBW  *numberOfLBnodes];
         dist.f[dirTW  ] = &distributionArray[dirBE  *numberOfLBnodes];
         dist.f[dirBE  ] = &distributionArray[dirTW  *numberOfLBnodes];
         dist.f[dirBS  ] = &distributionArray[dirTN  *numberOfLBnodes];
         dist.f[dirTN  ] = &distributionArray[dirBS  *numberOfLBnodes];
         dist.f[dirTS  ] = &distributionArray[dirBN  *numberOfLBnodes];
         dist.f[dirBN  ] = &distributionArray[dirTS  *numberOfLBnodes];
         dist.f[dirZERO] = &distributionArray[dirZERO*numberOfLBnodes];
         dist.f[dirTNE ] = &distributionArray[dirBSW *numberOfLBnodes];
         dist.f[dirTSW ] = &distributionArray[dirBNE *numberOfLBnodes];
         dist.f[dirTSE ] = &distributionArray[dirBNW *numberOfLBnodes];
         dist.f[dirTNW ] = &distributionArray[dirBSE *numberOfLBnodes];
         dist.f[dirBNE ] = &distributionArray[dirTSW *numberOfLBnodes];
         dist.f[dirBSW ] = &distributionArray[dirTNE *numberOfLBnodes];
         dist.f[dirBSE ] = &distributionArray[dirTNW *numberOfLBnodes];
         dist.f[dirBNW ] = &distributionArray[dirTSE *numberOfLBnodes];
    }
}

__inline__ __device__ void getPointersToSubgridDistances(SubgridDistances27& subgridD, real* subgridDistances, const unsigned int numberOfSubgridIndices)
{
    subgridD.q[dirE   ] = &subgridDistances[dirE    *numberOfSubgridIndices];
    subgridD.q[dirW   ] = &subgridDistances[dirW    *numberOfSubgridIndices];
    subgridD.q[dirN   ] = &subgridDistances[dirN    *numberOfSubgridIndices];
    subgridD.q[dirS   ] = &subgridDistances[dirS    *numberOfSubgridIndices];
    subgridD.q[dirT   ] = &subgridDistances[dirT    *numberOfSubgridIndices];
    subgridD.q[dirB   ] = &subgridDistances[dirB    *numberOfSubgridIndices];
    subgridD.q[dirNE  ] = &subgridDistances[dirNE   *numberOfSubgridIndices];
    subgridD.q[dirSW  ] = &subgridDistances[dirSW   *numberOfSubgridIndices];
    subgridD.q[dirSE  ] = &subgridDistances[dirSE   *numberOfSubgridIndices];
    subgridD.q[dirNW  ] = &subgridDistances[dirNW   *numberOfSubgridIndices];
    subgridD.q[dirTE  ] = &subgridDistances[dirTE   *numberOfSubgridIndices];
    subgridD.q[dirBW  ] = &subgridDistances[dirBW   *numberOfSubgridIndices];
    subgridD.q[dirBE  ] = &subgridDistances[dirBE   *numberOfSubgridIndices];
    subgridD.q[dirTW  ] = &subgridDistances[dirTW   *numberOfSubgridIndices];
    subgridD.q[dirTN  ] = &subgridDistances[dirTN   *numberOfSubgridIndices];
    subgridD.q[dirBS  ] = &subgridDistances[dirBS   *numberOfSubgridIndices];
    subgridD.q[dirBN  ] = &subgridDistances[dirBN   *numberOfSubgridIndices];
    subgridD.q[dirTS  ] = &subgridDistances[dirTS   *numberOfSubgridIndices];
    subgridD.q[dirZERO] = &subgridDistances[dirZERO *numberOfSubgridIndices];
    subgridD.q[dirTNE ] = &subgridDistances[dirTNE  *numberOfSubgridIndices];
    subgridD.q[dirTSW ] = &subgridDistances[dirTSW  *numberOfSubgridIndices];
    subgridD.q[dirTSE ] = &subgridDistances[dirTSE  *numberOfSubgridIndices];
    subgridD.q[dirTNW ] = &subgridDistances[dirTNW  *numberOfSubgridIndices];
    subgridD.q[dirBNE ] = &subgridDistances[dirBNE  *numberOfSubgridIndices];
    subgridD.q[dirBSW ] = &subgridDistances[dirBSW  *numberOfSubgridIndices];
    subgridD.q[dirBSE ] = &subgridDistances[dirBSE  *numberOfSubgridIndices];
    subgridD.q[dirBNW ] = &subgridDistances[dirBNW  *numberOfSubgridIndices];
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
