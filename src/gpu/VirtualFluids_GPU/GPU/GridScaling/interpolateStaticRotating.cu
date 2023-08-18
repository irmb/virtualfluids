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
//! \file scaleCF_compressible.cu
//! \ingroup GPU/GridScaling
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================

#include "DataTypes.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"
#include "LBM/GPUHelperFunctions/ChimeraTransformation.h"
#include "LBM/GPUHelperFunctions/ScalingUtilities.h"
#include "LBM/GPUHelperFunctions/GridTraversion.h"
#include "LBM/GPUHelperFunctions/CoordinateTransformation.h"

#include <lbm/refinement/Coefficients.h>
#include <lbm/refinement/InterpolationCF.h>

__global__ void interpolateStaticToRotating(
    unsigned int numberOfInterfaceNodes,
    unsigned int * indicesStaticCell,
    unsigned int *indicesRotating,
    const real *coordDestinationX,
    const real *coordDestinationY,
    const real *coordDestinationZ,
    const real *coordSourceX,
    const real *coordSourceY,
    const real *coordSourceZ,
    const uint *neighborXstatic,
    const uint *neighborYstatic,
    const uint *neighborZstatic,
    const uint *neighborMMMstatic,
    real centerCoordX,
    real centerCoordY,
    real centerCoordZ,
    real angleX,
    real angleY,
    real angleZ,
    real angularVelocityX,
    real angularVelocityY,
    real angularVelocityZ,
    real dx)
{
    const unsigned listIndex = vf::gpu::getNodeIndex();

    if (listIndex >= numberOfInterfaceNodes) return;

    uint destinationIndex = indicesRotating[listIndex];
    uint sourceIndex =  indicesStaticCell[listIndex];
    uint indexNeighborMMMsource = neighborMMMstatic[sourceIndex];

    real globalX;
    real globalY;
    real globalZ;

    transformRotatingToGlobal(globalX, globalY, globalZ, coordDestinationX[destinationIndex],
                              coordDestinationY[destinationIndex], coordDestinationZ[destinationIndex], centerCoordX,
                              centerCoordY, centerCoordZ, angleX, angleY, angleZ);

    // printf("coordDestinationX %.4f,  coordSourceX %.4f, d-s %.4f \n", globalX, coordSourceX[sourceIndex],
    //        globalX - coordSourceX[sourceIndex]);

    indicesStaticCell[listIndex] = traverseSourceCell(
        globalX, globalY, globalZ,
        indexNeighborMMMsource, coordSourceX[indexNeighborMMMsource], coordSourceY[indexNeighborMMMsource],
        coordSourceZ[indexNeighborMMMsource], neighborXstatic, neighborYstatic, neighborZstatic, dx);

    // if (sourceIndex !=  indicesStaticCell[listIndex]) {
    //     printf("sourceIndex: old = %d, new = %d\t", indicesStaticCell[listIndex], sourceIndex);
    // }
}

__global__ void updateGlobalCoordinates(
    unsigned int numberOfNodes,
    real *globalX,
    real *globalY,
    real *globalZ,
    const real *localX,
    const real *localY,
    const real *localZ,
    real centerCoordX,
    real centerCoordY,
    real centerCoordZ,
    real angleX,
    real angleY,
    real angleZ)
{
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if (nodeIndex >= numberOfNodes) return;

    real globalXtemp;
    real globalYtemp;
    real globalZtemp;

    transformRotatingToGlobal(globalXtemp, globalYtemp, globalZtemp, localX[nodeIndex], localY[nodeIndex], localZ[nodeIndex], centerCoordX,
                              centerCoordY, centerCoordZ, angleX, angleY, angleZ);

    // printf("%.3f\t", localY[nodeIndex]);

    globalX[nodeIndex] = globalXtemp;
    globalY[nodeIndex] = globalYtemp;
    globalZ[nodeIndex] = globalZtemp;
}