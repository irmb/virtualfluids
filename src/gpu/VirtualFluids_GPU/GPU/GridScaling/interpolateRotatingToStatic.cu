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

#include <lbm/MacroscopicQuantities.h>
#include <lbm/refinement/Coefficients.h>
#include <lbm/refinement/InterpolationCF.h>

using namespace vf::lbm;
using namespace vf::gpu;

__device__ __inline__ void calculateRotationalForces(real &forceX, real &forceY, real &forceZ, real angularVelocityX,
                                                     real angularVelocityY, real angularVelocityZ, real coordSourceXlocal,
                                                     real coordSourceYlocal, real coordSourceZlocal, real vx1, real vx2,
                                                     real vx3)
{
    forceX = -((-(angularVelocityY * angularVelocityY * coordSourceXlocal) -
                angularVelocityZ * angularVelocityZ * coordSourceXlocal +
                angularVelocityX * angularVelocityY * coordSourceYlocal -
                2 * angularVelocityX * angularVelocityX * angularVelocityZ * coordSourceYlocal -
                2 * angularVelocityY * angularVelocityY * angularVelocityZ * coordSourceYlocal -
                2 * angularVelocityZ * angularVelocityZ * angularVelocityZ * coordSourceYlocal +
                2 * angularVelocityX * angularVelocityX * angularVelocityY * coordSourceZlocal +
                2 * angularVelocityY * angularVelocityY * angularVelocityY * coordSourceZlocal +
                angularVelocityX * angularVelocityZ * coordSourceZlocal +
                2 * angularVelocityY * angularVelocityZ * angularVelocityZ * coordSourceZlocal +
                4 * angularVelocityY * angularVelocityY * vx1 + 4 * angularVelocityZ * angularVelocityZ * vx1 -
                4 * angularVelocityX * angularVelocityY * vx2 - 2 * angularVelocityZ * vx2 + 2 * angularVelocityY * vx3 -
                4 * angularVelocityX * angularVelocityZ * vx3) /
               (1 + 4 * angularVelocityX * angularVelocityX + 4 * angularVelocityY * angularVelocityY +
                4 * angularVelocityZ * angularVelocityZ));

    forceY = -((angularVelocityX * angularVelocityY * coordSourceXlocal +
                2 * angularVelocityX * angularVelocityX * angularVelocityZ * coordSourceXlocal +
                2 * angularVelocityY * angularVelocityY * angularVelocityZ * coordSourceXlocal +
                2 * angularVelocityZ * angularVelocityZ * angularVelocityZ * coordSourceXlocal -
                angularVelocityX * angularVelocityX * coordSourceYlocal -
                angularVelocityZ * angularVelocityZ * coordSourceYlocal -
                2 * angularVelocityX * angularVelocityX * angularVelocityX * coordSourceZlocal -
                2 * angularVelocityX * angularVelocityY * angularVelocityY * coordSourceZlocal +
                angularVelocityY * angularVelocityZ * coordSourceZlocal -
                2 * angularVelocityX * angularVelocityZ * angularVelocityZ * coordSourceZlocal -
                4 * angularVelocityX * angularVelocityY * vx1 + 2 * angularVelocityZ * vx1 +
                4 * angularVelocityX * angularVelocityX * vx2 + 4 * angularVelocityZ * angularVelocityZ * vx2 -
                2 * angularVelocityX * vx3 - 4 * angularVelocityY * angularVelocityZ * vx3) /
               (1 + 4 * angularVelocityX * angularVelocityX + 4 * angularVelocityY * angularVelocityY +
                4 * angularVelocityZ * angularVelocityZ));

    forceZ = -((-2 * angularVelocityX * angularVelocityX * angularVelocityY * coordSourceXlocal -
                2 * angularVelocityY * angularVelocityY * angularVelocityY * coordSourceXlocal +
                angularVelocityX * angularVelocityZ * coordSourceXlocal -
                2 * angularVelocityY * angularVelocityZ * angularVelocityZ * coordSourceXlocal +
                2 * angularVelocityX * angularVelocityX * angularVelocityX * coordSourceYlocal +
                2 * angularVelocityX * angularVelocityY * angularVelocityY * coordSourceYlocal +
                angularVelocityY * angularVelocityZ * coordSourceYlocal +
                2 * angularVelocityX * angularVelocityZ * angularVelocityZ * coordSourceYlocal -
                angularVelocityX * angularVelocityX * coordSourceZlocal -
                angularVelocityY * angularVelocityY * coordSourceZlocal - 2 * angularVelocityY * vx1 -
                4 * angularVelocityX * angularVelocityZ * vx1 + 2 * angularVelocityX * vx2 -
                4 * angularVelocityY * angularVelocityZ * vx2 + 4 * angularVelocityX * angularVelocityX * vx3 +
                4 * angularVelocityY * angularVelocityY * vx3) /
               (1 + 4 * angularVelocityX * angularVelocityX + 4 * angularVelocityY * angularVelocityY +
                4 * angularVelocityZ * angularVelocityZ));
}

__device__ __inline__ void calculateHalfRotationalForcesByVertex(real *forcesVertexX, real *forcesVertexY, real *forcesVertexZ,
                                                             uint index, InterpolationVertex direction,
                                                             const Distributions27 &distRotating, const real *coordSourceX,
                                                             const real *coordSourceY, const real *coordSourceZ,
                                                             real angularVelocityX, real angularVelocityY,
                                                             real angularVelocityZ)
{
    real drho, vx1, vx2, vx3;
    real coordSourceXlocal, coordSourceYlocal, coordSourceZlocal;
    real forceX, forceY, forceZ;

    getCompressibleMacroscopicValues(distRotating, index, drho, vx1, vx2, vx3);

    coordSourceXlocal = coordSourceX[index];
    coordSourceYlocal = coordSourceY[index];
    coordSourceZlocal = coordSourceZ[index];

    calculateRotationalForces(forceX, forceY, forceZ, angularVelocityX, angularVelocityY, angularVelocityZ,
                              coordSourceXlocal, coordSourceYlocal, coordSourceZlocal, vx1, vx2, vx3);
    forceX *= c1o2;
    forceY *= c1o2;
    forceZ *= c1o2;

    forcesVertexX[direction] = forceX;
    forcesVertexY[direction] = forceY;
    forcesVertexZ[direction] = forceZ;
}

__device__ __inline__ void calculateTangentialVelocities(real *tangentialVelocitiesX, real *tangentialVelocitiesY,
                                                         real *tangentialVelocitiesZ, uint index,
                                                         InterpolationVertex direction, const real *coordSourceX,
                                                         const real *coordSourceY, const real *coordSourceZ,
                                                         real angularVelocityX, real angularVelocityY, real angularVelocityZ, real dx)
{
    real coordSourceXlocal = coordSourceX[index] / dx;
    real coordSourceYlocal = coordSourceY[index] / dx;
    real coordSourceZlocal = coordSourceZ[index] / dx;
    tangentialVelocitiesX[direction] = -(angularVelocityZ * coordSourceYlocal) + angularVelocityY * coordSourceZlocal;
    tangentialVelocitiesY[direction] = angularVelocityZ * coordSourceXlocal - angularVelocityX * coordSourceZlocal;
    tangentialVelocitiesZ[direction] = -(angularVelocityY * coordSourceXlocal) + angularVelocityX * coordSourceYlocal;
}

__global__ void interpolateRotatingToStatic(
    real *distributionsStatic,
    real *distributionsRotating,
    unsigned int numberOfLBNodesStatic,
    unsigned int numberOfLBNodesRotating,
    unsigned int numberOfInterfaceNodes,
    const unsigned int *indicesStatic,
    unsigned int *indicesRotatingCell,
    const real *coordDestinationX,
    const real *coordDestinationY,
    const real *coordDestinationZ,
    const real *coordSourceX,
    const real *coordSourceY,
    const real *coordSourceZ,
    const uint *neighborXrotating,
    const uint *neighborYrotating,
    const uint *neighborZrotating,
    const uint *neighborMMMrotating,
    real centerCoordX,
    real centerCoordY,
    real centerCoordZ,
    real angleX,
    real angleY,
    real angleZ,
    real angularVelocityX,
    real angularVelocityY,
    real angularVelocityZ,
    real omegaRotating,
    bool isEvenTimestep,
    real dx)
{
    const unsigned listIndex = vf::gpu::getNodeIndex();
    if (listIndex >= numberOfInterfaceNodes) return;
    const uint destinationIndex = indicesStatic[listIndex];
    const uint previousSourceIndex = indicesRotatingCell[listIndex];
    const uint indexNeighborMMMsource = neighborMMMrotating[previousSourceIndex];

    real rotatedCoordDestinationX;
    real rotatedCoordDestinationY;
    real rotatedCoordDestinationZ;
    transformGlobalToRotating(rotatedCoordDestinationX, rotatedCoordDestinationY, rotatedCoordDestinationZ,
                              coordDestinationX[destinationIndex], coordDestinationY[destinationIndex],
                              coordDestinationZ[destinationIndex], centerCoordX, centerCoordY, centerCoordZ, angleX, angleY,
                              angleZ);

    const uint sourceIndex = traverseSourceCell(rotatedCoordDestinationX, rotatedCoordDestinationY, rotatedCoordDestinationZ,
                                                indexNeighborMMMsource, coordSourceX[indexNeighborMMMsource],
                                                coordSourceY[indexNeighborMMMsource], coordSourceZ[indexNeighborMMMsource],
                                                neighborXrotating, neighborYrotating, neighborZrotating, dx);

    indicesRotatingCell[listIndex] = sourceIndex;

    Distributions27 distRotating;
    vf::gpu::getPointersToDistributions(distRotating, distributionsRotating, numberOfLBNodesRotating, isEvenTimestep);

    real forcesVertexX [8];
    real forcesVertexY [8];
    real forcesVertexZ [8];
        real tangentialVelocitiesX [8];
    real tangentialVelocitiesY [8];
    real tangentialVelocitiesZ [8];

    // forces and tangential velocites
    // MMM
    uint indexTemp = sourceIndex;
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::MMM,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::MMM, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);
    // MMP
    indexTemp = neighborZrotating[indexTemp];
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::MMP,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::MMP, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);
    // PMP
    indexTemp = neighborXrotating[indexTemp];
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::PMP,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::PMP, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);
    // PMM
    indexTemp = neighborXrotating[sourceIndex]; // back to the base node of the cell --> sourceIndex
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::PMM,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::PMM, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);
    // PPM
    indexTemp = neighborYrotating[indexTemp];
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::PPM,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::PPM, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);
    // MPM
    indexTemp = neighborYrotating[sourceIndex]; // back to the base node of the cell --> sourceIndex
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::MPM,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::MPM, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);
    // MPP
    indexTemp = neighborZrotating[indexTemp];
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::MPP,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::MPP, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);
    // PPP
    indexTemp = neighborXrotating[indexTemp];
    calculateHalfRotationalForcesByVertex(forcesVertexX, forcesVertexY, forcesVertexZ, indexTemp, InterpolationVertex::PPP,
                                          distRotating, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                          angularVelocityY, angularVelocityZ);
    calculateTangentialVelocities(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ, indexTemp,
                                  InterpolationVertex::PPP, coordSourceX, coordSourceY, coordSourceZ, angularVelocityX,
                                  angularVelocityY, angularVelocityZ, dx);

    // 1. calculate moments for the nodes of the source cell, add half of the rotational force to the velocity during the computation of the moments
    vf::lbm::MomentsOnSourceNodeSet momentsSet;
    vf::gpu::calculateMomentSet<false>(momentsSet, listIndex, distributionsRotating, neighborXrotating, neighborYrotating,
                                       neighborZrotating, indicesRotatingCell, nullptr, numberOfLBNodesStatic, omegaRotating,
                                       isEvenTimestep, forcesVertexX, forcesVertexY, forcesVertexZ);

    // 2. calculate the coefficients for the interpolation
    // For this the relative coordinates of the destination pooint inside the source cell are needed
    // this cell coordinate system is centered at the middle of the source cell. The source cell has a side lenght of one
    // add tangential velocity to velocity for coefficient computaion
    const real cellCoordDestinationX = (coordSourceX[sourceIndex] + rotatedCoordDestinationX) / dx - c1o2;
    const real cellCoordDestinationY = (coordSourceY[sourceIndex] + rotatedCoordDestinationY) / dx - c1o2;
    const real cellCoordDestinationZ = (coordSourceZ[sourceIndex] + rotatedCoordDestinationZ) / dx - c1o2;
    momentsSet.addToVelocity(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ);
    vf::lbm::Coefficients coefficients;
    momentsSet.calculateCoefficients(coefficients, cellCoordDestinationX, cellCoordDestinationY, cellCoordDestinationZ);


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

    globalX[nodeIndex] = globalXtemp;
    globalY[nodeIndex] = globalYtemp;
    globalZ[nodeIndex] = globalZtemp;
}
