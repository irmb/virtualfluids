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


#include <lbm/CumulantChimeraK17.h>

namespace vf
{
namespace gpu
{

__device__ __host__ Distributions27 getDistributions27(real* distributions, unsigned int size_Mat, bool isEvenTimestep);

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

    Distributions27 dist;

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

__device__ unsigned int getNodeIndex();

__device__ bool isValidFluidNode(uint k, int size_Mat, uint nodeType);

__device__ void getLevelForce(real fx, real fy, real fz, int level, real* forces);

}
}

#endif
