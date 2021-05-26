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
/* Device code */
#include <lbm/CumulantChimeraPreCompiled.h>

#include "Kernel/Utilities/DistributionHelper.cuh"


namespace vf
{
namespace gpu 
{

__global__ void LB_Kernel_CumulantK17Unified(
    real omega,
    uint* typeOfGridNode,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    real* distributions,
    int size_Mat,
    int level,
    real* forces,
    bool isEvenTimestep)
{
    const uint k = getNodeIndex();
    const uint nodeType = typeOfGridNode[k];

    if (!isValidFluidNode(k, size_Mat, nodeType))
        return;

    DistributionWrapper distributionWrapper {
        distributions, size_Mat, isEvenTimestep, k, neighborX, neighborY, neighborZ
    };

    real level_forces[3];
    getLevelForce(forces[0], forces[1], forces[2], level, level_forces);

    vf::lbm::cumulantChimeraK17(distributionWrapper.distribution, omega, level_forces);

    distributionWrapper.write();
}


}
}
