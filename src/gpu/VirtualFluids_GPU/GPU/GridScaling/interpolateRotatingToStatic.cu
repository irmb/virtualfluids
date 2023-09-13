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
#include "GPU/GridScaling/interpolateRotatingToStaticInlines.h"

#include <lbm/MacroscopicQuantities.h>
#include <lbm/refinement/Coefficients.h>
#include <lbm/refinement/InterpolationCF.h>


using namespace vf::lbm;

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
    const uint *neighborXstatic,
    const uint *neighborYstatic,
    const uint *neighborZstatic,
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
                                       neighborZrotating, indicesRotatingCell, nullptr, numberOfLBNodesRotating, omegaRotating,
                                       isEvenTimestep, forcesVertexX, forcesVertexY, forcesVertexZ);

    // 2. calculate the coefficients for the interpolation
    // For this the relative coordinates of the source cell in the coordinate system of the destination node ("offsets")
    // this cell coordinate system is centered at the destination node. The source cell has a side lenght of one.
    // Also add the tangential velocity to velocity for coefficient computation
    const real cellCoordDestinationX = (coordSourceX[sourceIndex] - rotatedCoordDestinationX) / dx + c1o2;
    const real cellCoordDestinationY = (coordSourceY[sourceIndex] - rotatedCoordDestinationY) / dx + c1o2;
    const real cellCoordDestinationZ = (coordSourceZ[sourceIndex] - rotatedCoordDestinationZ) / dx + c1o2;
    momentsSet.addToVelocity(tangentialVelocitiesX, tangentialVelocitiesY, tangentialVelocitiesZ);
    vf::lbm::Coefficients coefficients;
    momentsSet.calculateCoefficients(coefficients, cellCoordDestinationX, cellCoordDestinationY, cellCoordDestinationZ);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Set all moments for the destination cell to zero
    //!
    real m111 = c0o1;
    real m211 = c0o1;
    real m011 = c0o1;
    real m121 = c0o1;
    real m101 = c0o1;
    real m112 = c0o1;
    real m110 = c0o1;
    real m221 = c0o1;
    real m001 = c0o1;
    real m201 = c0o1;
    real m021 = c0o1;
    real m212 = c0o1;
    real m010 = c0o1;
    real m210 = c0o1;
    real m012 = c0o1;
    real m122 = c0o1;
    real m100 = c0o1;
    real m120 = c0o1;
    real m102 = c0o1;
    real m222 = c0o1;
    real m022 = c0o1;
    real m202 = c0o1;
    real m002 = c0o1;
    real m220 = c0o1;
    real m020 = c0o1;
    real m200 = c0o1;
    real m000 = c0o1;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Define aliases to use the same variable for the distributions (f's):
    //!
    real& f000 = m111;
    real& fP00 = m211;
    real& fM00 = m011;
    real& f0P0 = m121;
    real& f0M0 = m101;
    real& f00P = m112;
    real& f00M = m110;
    real& fPP0 = m221;
    real& fMM0 = m001;
    real& fPM0 = m201;
    real& fMP0 = m021;
    real& fP0P = m212;
    real& fM0M = m010;
    real& fP0M = m210;
    real& fM0P = m012;
    real& f0PP = m122;
    real& f0MM = m100;
    real& f0PM = m120;
    real& f0MP = m102;
    real& fPPP = m222;
    real& fMPP = m022;
    real& fPMP = m202;
    real& fMMP = m002;
    real& fPPM = m220;
    real& fMPM = m020;
    real& fPMM = m200;
    real& fMMM = m000;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    real vvx, vvy, vvz;

    m000 = coefficients.d000; // m000 is press, if drho is interpolated directly
    vvx = coefficients.a000;
    vvy = coefficients.b000;
    vvz = coefficients.c000;

    ////////////////////////////////////////////////////////////////////////////////
    //! - rotate the velocities
    //!
    rotateDataFromGlobalToRotating(vvx, vvy, vvz, angleX, angleY, angleZ);

    ////////////////////////////////////////////////////////////////////////////////
    // calculate the squares of the velocities
    //
    real vxsq = vvx * vvx;
    real vysq = vvy * vvy;
    real vzsq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    // The second order moments have to be rotated. First they are computed in a rotating frame of reference

    real useNEQ = c1o1; // zero; //one;   //.... one = on ..... zero = off
    real mxxPyyPzz = m000;

    real mxxMyy = -c2o3 * ((coefficients.a100 - coefficients.b010)) / omegaRotating * (c1o1 + m000);
    real mxxMzz = -c2o3 * ((coefficients.a100 - coefficients.c001)) / omegaRotating * (c1o1 + m000);

    m011 = -c1o3 * ((coefficients.b001 + coefficients.c010)) / omegaRotating * (c1o1 + m000) * useNEQ;
    m101 = -c1o3 * ((coefficients.a001 + coefficients.c100)) / omegaRotating * (c1o1 + m000) * useNEQ;
    m110 = -c1o3 * ((coefficients.a010 + coefficients.b100)) / omegaRotating * (c1o1 + m000) * useNEQ;

    // rotate some second order moments
    rotateSecondOrderMomentsGlobalToRotating(m011, m101, m110, mxxMyy, mxxMzz, angleX, angleY, angleZ);

    // calculate remaining second order moments from previously rotated moments
    m200 = c1o3 * (        mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m020 = c1o3 * (-c2o1 * mxxMyy +        mxxMzz + mxxPyyPzz) * useNEQ;
    m002 = c1o3 * (        mxxMyy - c2o1 * mxxMzz + mxxPyyPzz) * useNEQ;

    // linear combinations for third order moments
    real mxxyPyzz, mxxyMyzz, mxxzPyyz, mxxzMyyz, mxyyPxzz, mxyyMxzz;
    m111 = c0o1;

    mxxyPyzz = c0o1;
    mxxyMyzz = c0o1;
    mxxzPyyz = c0o1;
    mxxzMyyz = c0o1;
    mxyyPxzz = c0o1;
    mxyyMxzz = c0o1;

    m210 = ( mxxyMyzz + mxxyPyzz) * c1o2;
    m012 = (-mxxyMyzz + mxxyPyzz) * c1o2;
    m201 = ( mxxzMyyz + mxxzPyyz) * c1o2;
    m021 = (-mxxzMyyz + mxxzPyyz) * c1o2;
    m120 = ( mxyyMxzz + mxyyPxzz) * c1o2;
    m102 = (-mxyyMxzz + mxyyPxzz) * c1o2;

    // fourth order moments
    m022 = m000 * c1o9;
    m202 = m022;
    m220 = m022;

    // fifth order moments

    // sixth order moments
    m222 = m000 * c1o27;

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Chimera transform from central moments to well conditioned distributions as defined in Appendix J in
    //! <a href="https://doi.org/10.1016/j.camwa.2015.05.001"><b>[ M. Geier et al. (2015),
    //! DOI:10.1016/j.camwa.2015.05.001 ]</b></a> see also Eq. (88)-(96) in <a
    //! href="https://doi.org/10.1016/j.jcp.2017.05.040"><b>[ M. Geier et al. (2017), DOI:10.1016/j.jcp.2017.05.040
    //! ]</b></a>
    //!
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    backwardInverseChimeraWithK(m000, m100, m200, vvx, vxsq, c1o1, c1o1);
    backwardChimera(            m010, m110, m210, vvx, vxsq);
    backwardInverseChimeraWithK(m020, m120, m220, vvx, vxsq, c3o1, c1o3);
    backwardChimera(            m001, m101, m201, vvx, vxsq);
    backwardChimera(            m011, m111, m211, vvx, vxsq);
    backwardChimera(            m021, m121, m221, vvx, vxsq);
    backwardInverseChimeraWithK(m002, m102, m202, vvx, vxsq, c3o1, c1o3);
    backwardChimera(            m012, m112, m212, vvx, vxsq);
    backwardInverseChimeraWithK(m022, m122, m222, vvx, vxsq, c9o1, c1o9);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    backwardInverseChimeraWithK(m000, m010, m020, vvy, vysq, c6o1, c1o6);
    backwardChimera(            m001, m011, m021, vvy, vysq);
    backwardInverseChimeraWithK(m002, m012, m022, vvy, vysq, c18o1, c1o18);
    backwardInverseChimeraWithK(m100, m110, m120, vvy, vysq, c3o2, c2o3);
    backwardChimera(            m101, m111, m121, vvy, vysq);
    backwardInverseChimeraWithK(m102, m112, m122, vvy, vysq, c9o2, c2o9);
    backwardInverseChimeraWithK(m200, m210, m220, vvy, vysq, c6o1, c1o6);
    backwardChimera(            m201, m211, m221, vvy, vysq);
    backwardInverseChimeraWithK(m202, m212, m222, vvy, vysq, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    backwardInverseChimeraWithK(m000, m001, m002, vvz, vzsq, c36o1, c1o36);
    backwardInverseChimeraWithK(m010, m011, m012, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m020, m021, m022, vvz, vzsq, c36o1, c1o36);
    backwardInverseChimeraWithK(m100, m101, m102, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m110, m111, m112, vvz, vzsq, c9o4,  c4o9);
    backwardInverseChimeraWithK(m120, m121, m122, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m200, m201, m202, vvz, vzsq, c36o1, c1o36);
    backwardInverseChimeraWithK(m210, m211, m212, vvz, vzsq, c9o1,  c1o9);
    backwardInverseChimeraWithK(m220, m221, m222, vvz, vzsq, c36o1, c1o36);


    // ////////////////////////////////////////////////////////////////////////////////

    real fStatic[27];
    fStatic[DIR_000] = f000;
    fStatic[DIR_P00] = fP00;
    fStatic[DIR_M00] = fM00;
    fStatic[DIR_0P0] = f0P0;
    fStatic[DIR_0M0] = f0M0;
    fStatic[DIR_00P] = f00P;
    fStatic[DIR_00M] = f00M;
    fStatic[DIR_PP0] = fPP0;
    fStatic[DIR_MM0] = fMM0;
    fStatic[DIR_PM0] = fPM0;
    fStatic[DIR_MP0] = fMP0;
    fStatic[DIR_P0P] = fP0P;
    fStatic[DIR_M0M] = fM0M;
    fStatic[DIR_P0M] = fP0M;
    fStatic[DIR_M0P] = fM0P;
    fStatic[DIR_0PP] = f0PP;
    fStatic[DIR_0MM] = f0MM;
    fStatic[DIR_0PM] = f0PM;
    fStatic[DIR_0MP] = f0MP;
    fStatic[DIR_PPP] = fPPP;
    fStatic[DIR_MPP] = fMPP;
    fStatic[DIR_PMP] = fPMP;
    fStatic[DIR_MMP] = fMMP;
    fStatic[DIR_PPM] = fPPM;
    fStatic[DIR_MPM] = fMPM;
    fStatic[DIR_PMM] = fPMM;
    fStatic[DIR_MMM] = fMMM;

    // get distribution pointers for destination node
    Distributions27 distStatic;
    vf::gpu::getPointersToDistributions(distStatic, distributionsStatic, numberOfLBNodesStatic, isEvenTimestep);

    // write
    vf::gpu::ListIndices indicesStaticForWriting(destinationIndex, neighborXstatic, neighborYstatic, neighborZstatic);
    vf::gpu::write(distStatic, indicesStaticForWriting, fStatic);
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

__global__ void traverseRotatingToStatic(
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
}