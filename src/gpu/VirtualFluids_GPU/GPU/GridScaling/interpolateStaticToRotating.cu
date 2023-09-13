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
#include "GPU/GridScaling/interpolateStaticToRotatingInlines.h"

#include <lbm/refinement/Coefficients.h>
#include <lbm/refinement/InterpolationCF.h>

using namespace vf::lbm;

__global__ void interpolateStaticToRotating(
    real *distributionsStatic,
    real *distributionsRotating,
    unsigned int numberOfLBNodesStatic,
    unsigned int numberOfLBNodesRotating,
    unsigned int numberOfInterfaceNodes,
    unsigned int *indicesStaticCell,
    const unsigned int *indicesRotating,
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
    const uint *neighborXrotating,
    const uint *neighborYrotating,
    const uint *neighborZrotating,
    real centerCoordX,
    real centerCoordY,
    real centerCoordZ,
    real angleX,
    real angleY,
    real angleZ,
    real angularVelocityX,
    real angularVelocityY,
    real angularVelocityZ,
    real omegaStatic,
    bool isEvenTimestep,
    real dx)
{
    // Interpolate from a cell on the static grid (source cell) to a node on the rotating grid (destination cell)

    // 1. calculate the indices of involved nodes
    const unsigned listIndex = vf::gpu::getNodeIndex();
    if (listIndex >= numberOfInterfaceNodes) return;
    const uint destinationIndex = indicesRotating[listIndex];
    const uint previousSourceIndex = indicesStaticCell[listIndex];
    const uint indexNeighborMMMsource = neighborMMMstatic[previousSourceIndex];

    // calculate the coordinates of the destination cell in the global coordinate system
    real globalCoordDestinationX;
    real globalCoordDestinationY;
    real globalCoordDestinationZ;
    transformRotatingToGlobal(globalCoordDestinationX, globalCoordDestinationY, globalCoordDestinationZ,
                              coordDestinationX[destinationIndex], coordDestinationY[destinationIndex],
                              coordDestinationZ[destinationIndex], centerCoordX, centerCoordY, centerCoordZ, angleX, angleY,
                              angleZ);

    // find the new index of the source cell (a static cell) after the rotation
    const uint sourceIndex =
        traverseSourceCell(globalCoordDestinationX, globalCoordDestinationY, globalCoordDestinationZ, indexNeighborMMMsource,
                           coordSourceX[indexNeighborMMMsource], coordSourceY[indexNeighborMMMsource],
                           coordSourceZ[indexNeighborMMMsource], neighborXstatic, neighborYstatic, neighborZstatic, dx);

    // write the new source index to the array
    indicesStaticCell[listIndex] = sourceIndex;

    // 1. calculate moments for the nodes of the source cell
    vf::lbm::MomentsOnSourceNodeSet momentsSet;
    vf::gpu::calculateMomentSet<false>(momentsSet, listIndex, distributionsStatic, neighborXstatic, neighborYstatic,
                                       neighborZstatic, indicesStaticCell, nullptr, numberOfLBNodesStatic, omegaStatic,
                                       isEvenTimestep);

    // 2. calculate the coefficients for the interpolation
    // For this the relative coordinates of the source cell in the coordinate system of the destination node ("offsets")
    // this cell coordinate system is centered at the destination node. The source cell has a side lenght of one.
    const real cellCoordDestinationX = (coordSourceX[sourceIndex] - globalCoordDestinationX) / dx + c1o2;
    const real cellCoordDestinationY = (coordSourceY[sourceIndex] - globalCoordDestinationY) / dx + c1o2;
    const real cellCoordDestinationZ = (coordSourceZ[sourceIndex] - globalCoordDestinationZ) / dx + c1o2;
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
    //! - Declare local variables for destination nodes
    //!
    real vvx, vvy, vvz, vxsq, vysq, vzsq;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set macroscopic values on destination node (zeroth and first order moments)
    //!
    m000  = coefficients.d000; // m000 is press, if drho is interpolated directly
    vvx   = coefficients.a000;
    vvy   = coefficients.b000;
    vvz   = coefficients.c000;

    ////////////////////////////////////////////////////////////////////////////////
    //! - The magic for the rotating grid begins here: subtract tangential velocity from the velocity in the inertial frame of reference
    //!
    real coordDestinationXlocal = coordDestinationX[destinationIndex] / dx;
    real coordDestinationYlocal = coordDestinationY[destinationIndex] / dx;
    real coordDestinationZlocal = coordDestinationZ[destinationIndex] / dx;

    vvx -= angularVelocityY * coordDestinationZlocal - angularVelocityZ * coordDestinationYlocal;
    vvy -= angularVelocityZ * coordDestinationXlocal - angularVelocityX * coordDestinationZlocal;
    vvz -= angularVelocityX * coordDestinationYlocal - angularVelocityY * coordDestinationXlocal;

    ////////////////////////////////////////////////////////////////////////////////
    //! - calculate the forces due to the rotating frame of reference
    //!
    // centrifugal force
    real forceX = angularVelocityY * angularVelocityY * coordDestinationXlocal +
                  angularVelocityZ * angularVelocityZ * coordDestinationXlocal -
                  angularVelocityX * angularVelocityY * coordDestinationYlocal -
                  angularVelocityX * angularVelocityZ * coordDestinationZlocal;
    real forceY = -angularVelocityX * angularVelocityY * coordDestinationXlocal +
                   angularVelocityX * angularVelocityX * coordDestinationYlocal +
                   angularVelocityZ * angularVelocityZ * coordDestinationYlocal -
                   angularVelocityY * angularVelocityZ * coordDestinationZlocal;
    real forceZ = -angularVelocityX * angularVelocityZ * coordDestinationXlocal -
                   angularVelocityY * angularVelocityZ * coordDestinationYlocal +
                   angularVelocityX * angularVelocityX * coordDestinationZlocal +
                   angularVelocityY * angularVelocityY * coordDestinationZlocal;

    // Coriolis force
    forceX +=  c2o1 * (angularVelocityZ * vvy - angularVelocityY * vvz);
    forceY += -c2o1 * (angularVelocityZ * vvx - angularVelocityX * vvz);
    forceZ +=  c2o1 * (angularVelocityY * vvx - angularVelocityX * vvy);

    ////////////////////////////////////////////////////////////////////////////////
    //! - subtract half the force from the velocities
    //!
    m100 -= c1o2 * forceX;
    m010 -= c1o2 * forceY;
    m001 -= c1o2 * forceZ;

    ////////////////////////////////////////////////////////////////////////////////
    //! - rotate the velocities
    //!
    rotateDataFromRotatingToGlobal(vvx, vvy, vvz, angleX, angleY, angleZ);

    ////////////////////////////////////////////////////////////////////////////////
    // calculate the squares of the velocities
    //
    vxsq = vvx * vvx;
    vysq = vvy * vvy;
    vzsq = vvz * vvz;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (second to sixth order) on destination node
    //!
    // linear combinations for second order moments
    // The second order moments have to be rotated. First they are computed in a static frame of reference

    real useNEQ = c1o1; // c0o1; // c1o1;    //.... one = on ..... zero = off
    real mxxPyyPzz = m000;

    real mxxMyy = -c2o3 * ((coefficients.a100 - coefficients.b010)) / omegaStatic * (c1o1 + m000);
    real mxxMzz = -c2o3 * ((coefficients.a100 - coefficients.c001)) / omegaStatic * (c1o1 + m000);

    m011 = -c1o3 * ((coefficients.b001 + coefficients.c010)) / omegaStatic * (c1o1 + m000) * useNEQ;
    m101 = -c1o3 * ((coefficients.a001 + coefficients.c100)) / omegaStatic * (c1o1 + m000) * useNEQ;
    m110 = -c1o3 * ((coefficients.a010 + coefficients.b100)) / omegaStatic * (c1o1 + m000) * useNEQ;

    // rotate some second order moments
    rotateSecondOrderMomentsRotatingToGlobal(m011, m101, m110, mxxMyy, mxxMzz, angleX, angleY, angleZ);

    // calculate the remaining second order moments from previously rotated moments
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

    real fRotating[27];
    fRotating[DIR_000] = f000;
    fRotating[DIR_P00] = fP00;
    fRotating[DIR_M00] = fM00;
    fRotating[DIR_0P0] = f0P0;
    fRotating[DIR_0M0] = f0M0;
    fRotating[DIR_00P] = f00P;
    fRotating[DIR_00M] = f00M;
    fRotating[DIR_PP0] = fPP0;
    fRotating[DIR_MM0] = fMM0;
    fRotating[DIR_PM0] = fPM0;
    fRotating[DIR_MP0] = fMP0;
    fRotating[DIR_P0P] = fP0P;
    fRotating[DIR_M0M] = fM0M;
    fRotating[DIR_P0M] = fP0M;
    fRotating[DIR_M0P] = fM0P;
    fRotating[DIR_0PP] = f0PP;
    fRotating[DIR_0MM] = f0MM;
    fRotating[DIR_0PM] = f0PM;
    fRotating[DIR_0MP] = f0MP;
    fRotating[DIR_PPP] = fPPP;
    fRotating[DIR_MPP] = fMPP;
    fRotating[DIR_PMP] = fPMP;
    fRotating[DIR_MMP] = fMMP;
    fRotating[DIR_PPM] = fPPM;
    fRotating[DIR_MPM] = fMPM;
    fRotating[DIR_PMM] = fPMM;
    fRotating[DIR_MMM] = fMMM;

    // get distribution pointers for destination node
    Distributions27 distRoating;
    vf::gpu::getPointersToDistributions(distRoating, distributionsRotating, numberOfLBNodesRotating, isEvenTimestep);

    // write
    vf::gpu::ListIndices indicesRotatingForWriting(destinationIndex, neighborXrotating, neighborYrotating,
                                                   neighborZrotating);
    vf::gpu::write(distRoating, indicesRotatingForWriting, fRotating);
}

__global__ void traverseStaticToRotating(
    unsigned int numberOfInterfaceNodes,
    unsigned int *indicesStaticCell,
    const unsigned int *indicesRotating,
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
    real dx)
{
    // Interpolate from a cell on the static grid (source cell) to a node on the rotating grid (destination cell)

    // 1. calculate the indices of involved nodes
    const unsigned listIndex = vf::gpu::getNodeIndex();
    if (listIndex >= numberOfInterfaceNodes) return;
    const uint destinationIndex = indicesRotating[listIndex];
    const uint previousSourceIndex = indicesStaticCell[listIndex];
    const uint indexNeighborMMMsource = neighborMMMstatic[previousSourceIndex];

    // calculate the coordinates of the destination cell in the global coordinate system
    real globalCoordDestinationX;
    real globalCoordDestinationY;
    real globalCoordDestinationZ;
    transformRotatingToGlobal(globalCoordDestinationX, globalCoordDestinationY, globalCoordDestinationZ,
                              coordDestinationX[destinationIndex], coordDestinationY[destinationIndex],
                              coordDestinationZ[destinationIndex], centerCoordX, centerCoordY, centerCoordZ, angleX, angleY,
                              angleZ);

    // find the new index of the source cell (a static cell) after the rotation
    const uint sourceIndex =
        traverseSourceCell(globalCoordDestinationX, globalCoordDestinationY, globalCoordDestinationZ, indexNeighborMMMsource,
                           coordSourceX[indexNeighborMMMsource], coordSourceY[indexNeighborMMMsource],
                           coordSourceZ[indexNeighborMMMsource], neighborXstatic, neighborYstatic, neighborZstatic, dx);

    // write the new source index to the array
    indicesStaticCell[listIndex] = sourceIndex;
}
