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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_PostProcessor PostProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Soeren Peters
//======================================================================================
#include "MacroscopicQuantities.cuh"

#include <helper_cuda.h>

#include <cuda_helper/CudaGrid.h>

#include <lbm/MacroscopicQuantities.h>
#include <lbm/constants/D3Q27.h>

#include <basics/constants/NumericConstants.h>

#include "Utilities/KernelUtilities.h"
#include "cuda_helper/CudaIndexCalculation.h"
#include "Calculation/Calculation.h"
#include "lbm/ChimeraTransformation.h"
#include "lbm/collision/TurbulentViscosity.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::cuda;

namespace vf::gpu {

__global__ void calculateMacroscopicQuantities_device(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();
   
    //////////////////////////////////////////////////////////////////////////
    if(nodeIndex<numberOfLBnodes)
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
       
        //////////////////////////////////////////////////////////////////////////
        //index
        unsigned int kzero= nodeIndex;
        unsigned int ke   = nodeIndex;
        unsigned int kw   = neighborX[nodeIndex];
        unsigned int kn   = nodeIndex;
        unsigned int ks   = neighborY[nodeIndex];
        unsigned int kt   = nodeIndex;
        unsigned int kb   = neighborZ[nodeIndex];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = nodeIndex;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = nodeIndex;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = nodeIndex;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
       
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            rhoD[nodeIndex] = 
                (dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[d000])[kzero]+ 
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vxD[nodeIndex] =
                (dist.f[dP00])[ke  ]- (dist.f[dM00])[kw  ]+ 
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]- (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]- (dist.f[dM0P])[ktw ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]- (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
           
            vyD[nodeIndex] =
                (dist.f[d0P0])[kn  ]- (dist.f[d0M0])[ks  ]+
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]-
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]- (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]- 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vzD[nodeIndex] =
                (dist.f[d00P])[kt  ]- (dist.f[d00M])[kb  ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]-
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]-
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]- 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
           
            pressD[nodeIndex] =
                ((dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                2.f*(
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ])+
                3.f*(
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+c0o1*rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]; // times zero for incompressible case   
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
       }
    }
}

__global__ void calculateMacroscopicQuantitiesCompressible_device(
    real *vxD,
    real *vyD,
    real *vzD,
    real *rhoD,
    real *pressD,
    const uint *geoD,
    const uint *neighborX,
    const uint *neighborY,
    const uint *neighborZ,
    unsigned long long numberOfLBnodes,
    real *distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    if(nodeIndex >= numberOfLBnodes)
        return;

    if (!isValidFluidNode(geoD[nodeIndex]))
    {
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
        return;
    }

    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
    ListIndices listIndices(nodeIndex, neighborX, neighborY, neighborZ);

    real distribution[27];
    getPreCollisionDistribution(distribution, dist, listIndices);

    rhoD[nodeIndex] = vf::lbm::getDensity(distribution);
    vxD[nodeIndex] = vf::lbm::getCompressibleVelocityX1(distribution, rhoD[nodeIndex]);
    vyD[nodeIndex] = vf::lbm::getCompressibleVelocityX2(distribution, rhoD[nodeIndex]);
    vzD[nodeIndex] = vf::lbm::getCompressibleVelocityX3(distribution, rhoD[nodeIndex]);
    pressD[nodeIndex] = vf::lbm::getPressure(distribution, rhoD[nodeIndex], vxD[nodeIndex], vyD[nodeIndex], vzD[nodeIndex]);
}


__global__ void calculateSubGridScaleFluxesCompressible_device(
    const uint* indices,
    const uint numberOfIndices,
    real *vxvx,
    real *vxvy,
    real *vxvz,
    real *vyvy,
    real *vyvz,
    real *vzvz,
    real *phix,
    real *phiy,
    real *phiz,
    const uint *geoD,
    const real* vx, const real* vy, const real* vz, 
    const real* scalars,
    const real* turbulenceViscosities,
    const real* turbulenceDiffusivities,
    const uint *neighborX,
    const uint *neighborY,
    const uint *neighborZ,
    const unsigned long long numberOfLBnodes,
    real *distributions,
    real* distributionsScalar,
    const real omega,
    const real omegaDiffusivity,
    const bool isEvenTimestep, 
    const bool computeSubGridScaleFluxesScalar)
{
    using namespace vf::lbm;

    const uint nodeIndex = get1DIndexFrom2DBlock();

    if(nodeIndex >= numberOfIndices)
        return;

    const uint k000 = indices[nodeIndex];

    if (!isValidFluidNode(geoD[k000])) {
        vxvx[k000] = c0o1;
        vxvy[k000] = c0o1;
        vxvz[k000] = c0o1;
        vyvy[k000] = c0o1;
        vyvz[k000] = c0o1;
        vzvz[k000] = c0o1;
        return;
    }

    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
    const ListIndices listIndices(k000, neighborX, neighborY, neighborZ);

    real population[27];
    getPreCollisionDistribution(population, dist, listIndices);

    real& m011 = population[dM00];
    real& m101 = population[d0M0];
    real& m110 = population[d00M];
    real& m002 = population[dMMP];
    real& m020 = population[dMPM];
    real& m200 = population[dPMM];
    real& m000 = population[dMMM];

    const real velocityX = vx[k000];
    const real velocityY = vy[k000];
    const real velocityZ = vz[k000];

    const real vx2 = velocityX * velocityX;
    const real vy2 = velocityY * velocityY;
    const real vz2 = velocityZ * velocityZ;

    forwardChimeraWithInverseK(population[dMMM], population[dMM0], population[dMMP], velocityZ, vz2, c36o1, c1o36);
    forwardChimeraWithInverseK(population[dM0M], population[dM00], population[dM0P], velocityZ, vz2, c9o1,  c1o9);
    forwardChimeraWithInverseK(population[dMPM], population[dMP0], population[dMPP], velocityZ, vz2, c36o1, c1o36);
    forwardChimeraWithInverseK(population[d0MM], population[d0M0], population[d0MP], velocityZ, vz2, c9o1,  c1o9);
    forwardChimeraWithInverseK(population[d00M], population[d000], population[d00P], velocityZ, vz2, c9o4,  c4o9);
    forwardChimeraWithInverseK(population[d0PM], population[d0P0], population[d0PP], velocityZ, vz2, c9o1,  c1o9);
    forwardChimeraWithInverseK(population[dPMM], population[dPM0], population[dPMP], velocityZ, vz2, c36o1, c1o36);
    forwardChimeraWithInverseK(population[dP0M], population[dP00], population[dP0P], velocityZ, vz2, c9o1,  c1o9);
    forwardChimeraWithInverseK(population[dPPM], population[dPP0], population[dPPP], velocityZ, vz2, c36o1, c1o36);

    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    forwardChimeraWithInverseK(population[dMMM], population[dM0M], population[dMPM], velocityY, vy2, c6o1,  c1o6); 
    forwardChimera(            population[dMM0], population[dM00], population[dMP0], velocityY, vy2);
    forwardChimeraWithInverseK(population[dMMP], population[dM0P], population[dMPP], velocityY, vy2, c18o1, c1o18);
    forwardChimeraWithInverseK(population[d0MM], population[d00M], population[d0PM], velocityY, vy2, c3o2,  c2o3);
    forwardChimera(            population[d0M0], population[d000], population[d0P0], velocityY, vy2);
    forwardChimeraWithInverseK(population[d0MP], population[d00P], population[d0PP], velocityY, vy2, c9o2,  c2o9);
    forwardChimeraWithInverseK(population[dPMM], population[dP0M], population[dPPM], velocityY, vy2, c6o1,  c1o6);
    forwardChimera(            population[dPM0], population[dP00], population[dPP0], velocityY, vy2);
    forwardChimeraWithInverseK(population[dPMP], population[dP0P], population[dPPP], velocityY, vy2, c18o1, c1o18);

    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    forwardChimeraWithInverseK(population[dMMM], population[d0MM], population[dPMM], velocityX, vx2, c1o1, c1o1); 
    forwardChimera(            population[dM0M], population[d00M], population[dP0M], velocityX, vx2);
    forwardChimeraWithInverseK(population[dMPM], population[d0PM], population[dPPM], velocityX, vx2, c3o1, c1o3);
    forwardChimera(            population[dMM0], population[d0M0], population[dPM0], velocityX, vx2);
    forwardChimera(            population[dM00], population[d000], population[dP00], velocityX, vx2);
    // forwardChimera(            population[dMP0], population[d0P0], population[dPP0], vvx, vx2);
    forwardChimeraWithInverseK(population[dMMP], population[d0MP], population[dPMP], velocityX, vx2, c3o1, c1o3);
    // forwardChimera(            population[dM0P], population[d00P], population[dP0P], vvx, vx2);
    // forwardChimeraWithInverseK(population[dMPP], population[d0PP], population[dPPP], vvx, vx2, c3o1, c1o3);
    
    const real OxxPyyPzz = c1o1;

    const real mxxPyyPzz = m200 + m020 + m002;
    const real mxxMyy    = m200 - m020;
    const real mxxMzz    = m200 - m002;

    const real turbulenceViscosity = turbulenceViscosities[k000];
    const real omegaTurbulent = vf::lbm::calculateOmegaWithTurbulentViscosity(omega, turbulenceViscosity);
    const real turbulenceViscosityOverDensity = turbulenceViscosity / (c1o1 + m000);

    const real Dxy = -c3o1 * omegaTurbulent * m110;
    const real Dxz = -c3o1 * omegaTurbulent * m101;
    const real Dyz = -c3o1 * omegaTurbulent * m011;
    const real dxux = c1o2 * (-omegaTurbulent) * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (m000 - mxxPyyPzz);
    const real dyuy = dxux + omegaTurbulent * c3o2 * mxxMyy;
    const real dzuz = dxux + omegaTurbulent * c3o2 * mxxMzz;

    vxvx[nodeIndex] = -turbulenceViscosityOverDensity * dxux;
    vxvy[nodeIndex] = -turbulenceViscosityOverDensity * Dxy;
    vxvz[nodeIndex] = -turbulenceViscosityOverDensity * Dxz;
    vyvy[nodeIndex] = -turbulenceViscosityOverDensity * dyuy;
    vyvz[nodeIndex] = -turbulenceViscosityOverDensity * Dyz;
    vzvz[nodeIndex] = -turbulenceViscosityOverDensity * dzuz;

    if (!computeSubGridScaleFluxesScalar)
        return;

    getPointersToDistributions(dist, distributionsScalar, numberOfLBnodes, isEvenTimestep);
    getPreCollisionDistribution(population, dist, listIndices);

    const real scalar = scalars[k000];
    const real turbulenceDiffusivity = turbulenceDiffusivities[k000];
    const real omegaDiffusiveTurbulent = vf::lbm::calculateOmegaWithTurbulentViscosity(omegaDiffusivity, turbulenceDiffusivity);
    phix[nodeIndex] = -turbulenceDiffusivity * c3o1 * omegaDiffusiveTurbulent * (vf::lbm::getIncompressibleVelocityX1(population) - scalar * velocityX);
    phiy[nodeIndex] = -turbulenceDiffusivity * c3o1 * omegaDiffusiveTurbulent * (vf::lbm::getIncompressibleVelocityX2(population) - scalar * velocityY);
    phiz[nodeIndex] = -turbulenceDiffusivity * c3o1 * omegaDiffusiveTurbulent * (vf::lbm::getIncompressibleVelocityX3(population) - scalar * velocityZ);
}

__global__ void calculateMean_device(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    //////////////////////////////////////////////////////////////////////////
    if( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
        
        //////////////////////////////////////////////////////////////////////////
        //index
        unsigned int kzero= nodeIndex;
        unsigned int ke   = nodeIndex;
        unsigned int kw   = neighborX[nodeIndex];
        unsigned int kn   = nodeIndex;
        unsigned int ks   = neighborY[nodeIndex];
        unsigned int kt   = nodeIndex;
        unsigned int kb   = neighborZ[nodeIndex];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = nodeIndex;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = nodeIndex;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = nodeIndex;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
        
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            rhoD[nodeIndex] =
                (dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[d000])[kzero]+ 
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw]+
                RHO;
            
            vxD[nodeIndex] =
                (dist.f[dP00])[ke  ]- (dist.f[dM00])[kw  ]+ 
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]- (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]- (dist.f[dM0P])[ktw ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]- (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw]+
                VX;
            
            vyD[nodeIndex] =
                (dist.f[d0P0])[kn  ]- (dist.f[d0M0])[ks  ]+
                (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]-
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]- (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]- 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw]+
                VY;
            
            vzD[nodeIndex] =
                (dist.f[d00P])[kt  ]- (dist.f[d00M])[kb  ]+
                (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]-
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]-
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]- 
                (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw]+
                VZ;
            
            pressD[nodeIndex] =
                ((dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                c2o1*(
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ])+
                c3o1*(
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]+
                PRESS;    
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
        }
    }
}

__global__ void calculateMeanCompressible_device(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    //////////////////////////////////////////////////////////////////////////
    if( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
        
        //////////////////////////////////////////////////////////////////////////
        //index
        //unsigned int kzero= k;
        unsigned int ke   = nodeIndex;
        unsigned int kw   = neighborX[nodeIndex];
        unsigned int kn   = nodeIndex;
        unsigned int ks   = neighborY[nodeIndex];
        unsigned int kt   = nodeIndex;
        unsigned int kb   = neighborZ[nodeIndex];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = nodeIndex;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = nodeIndex;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = nodeIndex;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
        
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            real mfcbb = (dist.f[dP00])[nodeIndex];//[ke   ];
            real mfabb = (dist.f[dM00])[kw];//[kw   ];  
            real mfbcb = (dist.f[d0P0])[nodeIndex];//[kn   ];
            real mfbab = (dist.f[d0M0])[ks];//[ks   ];  
            real mfbbc = (dist.f[d00P])[nodeIndex];//[kt   ];
            real mfbba = (dist.f[d00M])[kb];//[kb   ];  
            real mfccb = (dist.f[dPP0])[nodeIndex];//[kne  ];  
            real mfaab = (dist.f[dMM0])[ksw];//[ksw  ];
            real mfcab = (dist.f[dPM0])[ks];//[kse  ]; 
            real mfacb = (dist.f[dMP0])[kw];//[knw  ]; 
            real mfcbc = (dist.f[dP0P])[nodeIndex];//[kte  ];  
            real mfaba = (dist.f[dM0M])[kbw];//[kbw  ];
            real mfcba = (dist.f[dP0M])[kb];//[kbe  ]; 
            real mfabc = (dist.f[dM0P])[kw];//[ktw  ]; 
            real mfbcc = (dist.f[d0PP])[nodeIndex];//[ktn  ];  
            real mfbaa = (dist.f[d0MM])[kbs];//[kbs  ];
            real mfbca = (dist.f[d0PM])[kb];//[kbn  ]; 
            real mfbac = (dist.f[d0MP])[ks];//[kts  ]; 
            real mfbbb = (dist.f[d000])[nodeIndex];//[kzero];
            real mfccc = (dist.f[dPPP])[nodeIndex];//[ktne ]; 
            real mfaac = (dist.f[dMMP])[ksw];//[ktsw ]; 
            real mfcac = (dist.f[dPMP])[ks];//[ktse ];
            real mfacc = (dist.f[dMPP])[kw];//[ktnw ];
            real mfcca = (dist.f[dPPM])[kb];//[kbne ];
            real mfaaa = (dist.f[dMMM])[kbsw];//[kbsw ];
            real mfcaa = (dist.f[dPMM])[kbs];//[kbse ]; 
            real mfaca = (dist.f[dMPM])[kbw];//[kbnw ]; 
            ////////////////////////////////////////////////////////////////////////////////////
            real drho = 
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

            real rho = c1o1 + drho;

            rhoD[nodeIndex] = drho + RHO;

            vxD[nodeIndex] = 
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                (mfcbb - mfabb)) / rho) + VX;
            vyD[nodeIndex] = 
                (((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                (mfbcb - mfbab)) / rho) + VY;
            vzD[nodeIndex] = 
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                (mfbbc - mfbba)) / rho) + VZ;

            pressD[nodeIndex]  =
                ((dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                c2o1*(
                (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ])+
                c3o1*(
                (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw])-
                rhoD[nodeIndex]-(vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1+rhoD[nodeIndex])) * c1o2+rhoD[nodeIndex]+
                PRESS;    
            //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
        }
    }
}

__global__ void calculateMeanCompressibleAdvectionDiffusion_device(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    real* concD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    real* distributionsAD,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    //////////////////////////////////////////////////////////////////////////
    if ( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist, distAD;
        getPointersToDistributions(dist,   distributions,   numberOfLBnodes, isEvenTimestep);
        getPointersToDistributions(distAD, distributionsAD, numberOfLBnodes, isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////
        //index
        //unsigned int kzero = k;
        unsigned int ke = nodeIndex;
        unsigned int kw = neighborX[nodeIndex];
        unsigned int kn = nodeIndex;
        unsigned int ks = neighborY[nodeIndex];
        unsigned int kt = nodeIndex;
        unsigned int kb = neighborZ[nodeIndex];
        unsigned int ksw = neighborY[kw];
        unsigned int kne = nodeIndex;
        unsigned int kse = ks;
        unsigned int knw = kw;
        unsigned int kbw = neighborZ[kw];
        unsigned int kte = nodeIndex;
        unsigned int kbe = kb;
        unsigned int ktw = kw;
        unsigned int kbs = neighborZ[ks];
        unsigned int ktn = nodeIndex;
        unsigned int kbn = kb;
        unsigned int kts = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = nodeIndex;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        real CONC  = concD[nodeIndex];
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        concD[nodeIndex]  = c0o1;
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
        
        if (geoD[nodeIndex] == GEO_FLUID)
        {
            real mfcbb = (dist.f[dP00])[nodeIndex];//[ke   ];
            real mfabb = (dist.f[dM00])[kw];//[kw   ];  
            real mfbcb = (dist.f[d0P0])[nodeIndex];//[kn   ];
            real mfbab = (dist.f[d0M0])[ks];//[ks   ];  
            real mfbbc = (dist.f[d00P])[nodeIndex];//[kt   ];
            real mfbba = (dist.f[d00M])[kb];//[kb   ];  
            real mfccb = (dist.f[dPP0])[nodeIndex];//[kne  ];  
            real mfaab = (dist.f[dMM0])[ksw];//[ksw  ];
            real mfcab = (dist.f[dPM0])[ks];//[kse  ]; 
            real mfacb = (dist.f[dMP0])[kw];//[knw  ]; 
            real mfcbc = (dist.f[dP0P])[nodeIndex];//[kte  ];  
            real mfaba = (dist.f[dM0M])[kbw];//[kbw  ];
            real mfcba = (dist.f[dP0M])[kb];//[kbe  ]; 
            real mfabc = (dist.f[dM0P])[kw];//[ktw  ]; 
            real mfbcc = (dist.f[d0PP])[nodeIndex];//[ktn  ];  
            real mfbaa = (dist.f[d0MM])[kbs];//[kbs  ];
            real mfbca = (dist.f[d0PM])[kb];//[kbn  ]; 
            real mfbac = (dist.f[d0MP])[ks];//[kts  ]; 
            real mfbbb = (dist.f[d000])[nodeIndex];//[kzero];
            real mfccc = (dist.f[dPPP])[nodeIndex];//[ktne ]; 
            real mfaac = (dist.f[dMMP])[ksw];//[ktsw ]; 
            real mfcac = (dist.f[dPMP])[ks];//[ktse ];
            real mfacc = (dist.f[dMPP])[kw];//[ktnw ];
            real mfcca = (dist.f[dPPM])[kb];//[kbne ];
            real mfaaa = (dist.f[dMMM])[kbsw];//[kbsw ];
            real mfcaa = (dist.f[dPMM])[kbs];//[kbse ]; 
            real mfaca = (dist.f[dMPM])[kbw];//[kbnw ]; 
            ////////////////////////////////////////////////////////////////////////////////////
            real drho =
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) + mfbbb;
            real rho = c1o1 + drho;
            ////////////////////////////////////////////////////////////////////////////////////
            
            rhoD[nodeIndex] = drho + RHO;
            
            vxD[nodeIndex] =
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                    (mfcbb - mfabb)) / rho) + VX;
            
            vyD[nodeIndex] =
                (((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                    (mfbcb - mfbab)) / rho) + VY;
            
            vzD[nodeIndex] =
                (((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                    (mfbbc - mfbba)) / rho) + VZ;
            
            pressD[nodeIndex] = 
                ((dist.f[dP00])[ke] + (dist.f[dM00])[kw] +
                 (dist.f[d0P0])[kn] + (dist.f[d0M0])[ks] +
                 (dist.f[d00P])[kt] + (dist.f[d00M])[kb] +
                 c2o1*(
                 (dist.f[dPP0])[kne] + (dist.f[dMM0])[ksw] +
                 (dist.f[dPM0])[kse] + (dist.f[dMP0])[knw] +
                 (dist.f[dP0P])[kte] + (dist.f[dM0M])[kbw] +
                 (dist.f[dP0M])[kbe] + (dist.f[dM0P])[ktw] +
                 (dist.f[d0PP])[ktn] + (dist.f[d0MM])[kbs] +
                 (dist.f[d0PM])[kbn] + (dist.f[d0MP])[kts]) +
                 c3o1*(
                 (dist.f[dPPP])[ktne] + (dist.f[dMMP])[ktsw] +
                 (dist.f[dPMP])[ktse] + (dist.f[dMPP])[ktnw] +
                 (dist.f[dPPM])[kbne] + (dist.f[dMMM])[kbsw] +
                 (dist.f[dPMM])[kbse] + (dist.f[dMPM])[kbnw]) -
                 rhoD[nodeIndex] - (vxD[nodeIndex] * vxD[nodeIndex] + vyD[nodeIndex] * vyD[nodeIndex] + vzD[nodeIndex] * vzD[nodeIndex]) * (c1o1 + rhoD[nodeIndex])) * c1o2 + rhoD[nodeIndex] +
                 PRESS;
                 //achtung op hart gesetzt Annahme op = 1 ;                                                    ^^^^(1.0/op-0.5)=0.5
            //////////////////////////////////////////////////////////////////////////
            mfcbb = (distAD.f[dP00])[nodeIndex   ];
            mfabb = (distAD.f[dM00])[kw  ];
            mfbcb = (distAD.f[d0P0])[nodeIndex   ];
            mfbab = (distAD.f[d0M0])[ks  ];
            mfbbc = (distAD.f[d00P])[nodeIndex   ];
            mfbba = (distAD.f[d00M])[kb  ];
            mfccb = (distAD.f[dPP0])[nodeIndex   ];
            mfaab = (distAD.f[dMM0])[ksw ];
            mfcab = (distAD.f[dPM0])[ks  ];
            mfacb = (distAD.f[dMP0])[kw  ];
            mfcbc = (distAD.f[dP0P])[nodeIndex   ];
            mfaba = (distAD.f[dM0M])[kbw ];
            mfcba = (distAD.f[dP0M])[kb  ];
            mfabc = (distAD.f[dM0P])[kw  ];
            mfbcc = (distAD.f[d0PP])[nodeIndex   ];
            mfbaa = (distAD.f[d0MM])[kbs ];
            mfbca = (distAD.f[d0PM])[kb  ];
            mfbac = (distAD.f[d0MP])[ks  ];
            mfbbb = (distAD.f[d000])[nodeIndex   ];
            mfccc = (distAD.f[dPPP])[nodeIndex   ];
            mfaac = (distAD.f[dMMP])[ksw ];
            mfcac = (distAD.f[dPMP])[ks  ];
            mfacc = (distAD.f[dMPP])[kw  ];
            mfcca = (distAD.f[dPPM])[kb  ];
            mfaaa = (distAD.f[dMMM])[kbsw];
            mfcaa = (distAD.f[dPMM])[kbs ];
            mfaca = (distAD.f[dMPM])[kbw ];
            //////////////////////////////////////////////////////////////////////////
            concD[nodeIndex] = 
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa)   + (mfaac + mfcca))) +
                 (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba)   + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                  ((mfabb + mfcbb) + (mfbab + mfbcb)  +  (mfbba + mfbbc))) +  mfbbb + CONC;
        }
    }
}

__global__ void calculateMacrosopicMean_device(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned int tdiff,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    //////////////////////////////////////////////////////////////////////////
    if(nodeIndex<numberOfLBnodes)
    {
        //////////////////////////////////////////////////////////////////////////
        real PRESS = pressD[nodeIndex];
        real RHO   = rhoD[nodeIndex];
        real VX    = vxD[nodeIndex];
        real VY    = vyD[nodeIndex];
        real VZ    = vzD[nodeIndex];
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
       
        if(geoD[nodeIndex] == GEO_FLUID)
        {
            rhoD[nodeIndex]    =   RHO   / tdiff;
            vxD[nodeIndex]     =   VX    / tdiff;
            vyD[nodeIndex]     =   VY    / tdiff;
            vzD[nodeIndex]     =   VZ    / tdiff;
            pressD[nodeIndex]  =   PRESS / tdiff;    
        }
    }
}

__global__ void LBResetMeanValuesSP27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    //////////////////////////////////////////////////////////////////////////
    if ( nodeIndex < numberOfLBnodes )
    {
        //////////////////////////////////////////////////////////////////////////
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex] = c0o1;
        vxD[nodeIndex] = c0o1;
        vyD[nodeIndex] = c0o1;
        vzD[nodeIndex] = c0o1;
    }
}

__global__ void LBResetMeanValuesAD27(
    real* vxD,
    real* vyD,
    real* vzD,
    real* rhoD,
    real* pressD,
    real* concD,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    //////////////////////////////////////////////////////////////////////////
    if (nodeIndex < numberOfLBnodes)
    {
        concD[nodeIndex]  = c0o1;
        pressD[nodeIndex] = c0o1;
        rhoD[nodeIndex]   = c0o1;
        vxD[nodeIndex]    = c0o1;
        vyD[nodeIndex]    = c0o1;
        vzD[nodeIndex]    = c0o1;
    }
}

__global__ void LBCalcMeasurePoints(
    real* vxMP,
    real* vyMP,
    real* vzMP,
    real* rhoMP,
    unsigned int* kMP,
    unsigned int numberOfPointskMP,
    unsigned int MPClockCycle,
    unsigned int t,
    unsigned int* geoD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned long long numberOfLBnodes,
    real* distributions,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = get1DIndexFrom2DBlock();

    //////////////////////////////////////////////////////////////////////////
    if( nodeIndex < numberOfPointskMP )
    {
        //////////////////////////////////////////////////////////////////////////
        //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
        //! timestep is based on the esoteric twist algorithm \ref <a
        //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
        //! DOI:10.3390/computation5020019 ]</b></a>
        //!
        Distributions27 dist;
        getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

        //////////////////////////////////////////////////////////////////////////
        //index
        unsigned int kzero= kMP[nodeIndex];//k;
        unsigned int ke   = kzero;
        unsigned int kw   = neighborX[kzero];
        unsigned int kn   = kzero;
        unsigned int ks   = neighborY[kzero];
        unsigned int kt   = kzero;
        unsigned int kb   = neighborZ[kzero];
        unsigned int ksw  = neighborY[kw];
        unsigned int kne  = kzero;
        unsigned int kse  = ks;
        unsigned int knw  = kw;
        unsigned int kbw  = neighborZ[kw];
        unsigned int kte  = kzero;
        unsigned int kbe  = kb;
        unsigned int ktw  = kw;
        unsigned int kbs  = neighborZ[ks];
        unsigned int ktn  = kzero;
        unsigned int kbn  = kb;
        unsigned int kts  = ks;
        unsigned int ktse = ks;
        unsigned int kbnw = kbw;
        unsigned int ktnw = kw;
        unsigned int kbse = kbs;
        unsigned int ktsw = ksw;
        unsigned int kbne = kb;
        unsigned int ktne = kzero;
        unsigned int kbsw = neighborZ[ksw];
        //////////////////////////////////////////////////////////////////////////
        unsigned int kMac = nodeIndex*MPClockCycle + t;
        //////////////////////////////////////////////////////////////////////////
        
        if(geoD[kzero] == GEO_FLUID)
        {
            rhoMP[kMac]= (dist.f[dP00])[ke  ]+ (dist.f[dM00])[kw  ]+ 
                         (dist.f[d0P0])[kn  ]+ (dist.f[d0M0])[ks  ]+
                         (dist.f[d00P])[kt  ]+ (dist.f[d00M])[kb  ]+
                         (dist.f[dPP0])[kne ]+ (dist.f[dMM0])[ksw ]+
                         (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                         (dist.f[dP0P])[kte ]+ (dist.f[dM0M])[kbw ]+
                         (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                         (dist.f[d0PP])[ktn ]+ (dist.f[d0MM])[kbs ]+
                         (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                         (dist.f[d000])[kzero]+ 
                         (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                         (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                         (dist.f[dPPM])[kbne]+ (dist.f[dMMM])[kbsw]+ 
                         (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vxMP[kMac] = (dist.f[dP00])[ke  ]- (dist.f[dM00])[kw  ]+ 
                         (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]+
                         (dist.f[dPM0])[kse ]- (dist.f[dMP0])[knw ]+
                         (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]+
                         (dist.f[dP0M])[kbe ]- (dist.f[dM0P])[ktw ]+
                         (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]+ 
                         (dist.f[dPMP])[ktse]- (dist.f[dMPP])[ktnw]+ 
                         (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]+ 
                         (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
           
            vyMP[kMac] = (dist.f[d0P0])[kn  ]- (dist.f[d0M0])[ks  ]+
                         (dist.f[dPP0])[kne ]- (dist.f[dMM0])[ksw ]-
                         (dist.f[dPM0])[kse ]+ (dist.f[dMP0])[knw ]+
                         (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]+
                         (dist.f[d0PM])[kbn ]- (dist.f[d0MP])[kts ]+
                         (dist.f[dPPP])[ktne]- (dist.f[dMMP])[ktsw]- 
                         (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]+ 
                         (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                         (dist.f[dPMM])[kbse]+ (dist.f[dMPM])[kbnw];
           
            vzMP[kMac] = (dist.f[d00P])[kt  ]- (dist.f[d00M])[kb  ]+
                         (dist.f[dP0P])[kte ]- (dist.f[dM0M])[kbw ]-
                         (dist.f[dP0M])[kbe ]+ (dist.f[dM0P])[ktw ]+
                         (dist.f[d0PP])[ktn ]- (dist.f[d0MM])[kbs ]-
                         (dist.f[d0PM])[kbn ]+ (dist.f[d0MP])[kts ]+
                         (dist.f[dPPP])[ktne]+ (dist.f[dMMP])[ktsw]+ 
                         (dist.f[dPMP])[ktse]+ (dist.f[dMPP])[ktnw]- 
                         (dist.f[dPPM])[kbne]- (dist.f[dMMM])[kbsw]- 
                         (dist.f[dPMM])[kbse]- (dist.f[dMPM])[kbnw];
        }
    }
}

void calculateMacroscopicQuantities(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                 unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                 unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    calculateMacroscopicQuantities_device<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                               numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("calculateMacroscopicQuantities_device execution failed");
}

void calculateMacroscopicQuantitiesCompressible(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, const uint* geoD, const uint* neighborX,
                     const uint* neighborY, const uint* neighborZ, unsigned long long numberOfLBnodes,
                     unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    calculateMacroscopicQuantitiesCompressible_device<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                                   numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("calculateMacroscopicQuantitiesCompressible_device execution failed");
}

void calculateSubGridScaleFluxesCompressible(const uint* indices, const uint numberOfIndices, real* vxvx, real* vxvy, real* vxvz, real* vyvy, real* vyvz, real* vzvz,
                                             real* phix, real* phiy, real* phiz, const uint* geoD, const real* vx,
                                             const real* vy, const real* vz, const real* scalars,
                                             const real* turbulenceViscosities, const real* turbulenceDiffusivities,
                                             const uint* neighborX, const uint* neighborY, const uint* neighborZ,
                                             const unsigned long long numberOfLBnodes, real* distributions,
                                             real* distributionsScalar, const real omega, const real omegaDiffusivity, const uint numberOfThreads,
                                             const bool isEvenTimestep, const bool computeSubGridScaleFluxesScalar)
{
    vf::cuda::CudaGrid grid(numberOfThreads, numberOfIndices);
    calculateSubGridScaleFluxesCompressible_device<<<grid.grid, grid.threads>>>(indices, numberOfIndices,
        vxvx, vxvy, vxvz, vyvy, vyvz, vzvz, phix, phiy, phiz, geoD, vx, vy, vz, scalars, turbulenceViscosities,
        turbulenceDiffusivities, neighborX, neighborY, neighborZ, numberOfLBnodes, distributions, distributionsScalar, omega, omegaDiffusivity,
        isEvenTimestep, computeSubGridScaleFluxesScalar);
    getLastCudaError("calculateSubgridScaleFluxesCompressible_device execution failed");
}

void calculateMean(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                 unsigned int* neighborY, unsigned int* neighborZ, unsigned long long numberOfLBnodes,
                 unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    calculateMean_device<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                               numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("calculateMean_device execution failed");
}

void calculateMeanCompressible(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, uint* geoD, uint* neighborX,
                     uint* neighborY, uint* neighborZ, unsigned long long numberOfLBnodes,
                     unsigned int numberOfThreads, real* DD, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    calculateMeanCompressible_device<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ,
                                                   numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("calculateMeanCompressible_device execution failed");
}

void calculateMeanCompressibleAdvectionDiffusion(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, real* concD, unsigned int* geoD,
                     unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                     unsigned long long numberOfLBnodes, unsigned int numberOfThreads, real* DD, real* DD_AD,
                     bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    calculateMeanCompressibleAdvectionDiffusion_device<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, concD, geoD, neighborX, neighborY, neighborZ,
                                                   numberOfLBnodes, DD, DD_AD, isEvenTimestep);
    getLastCudaError("calculateMeanCompressibleAdvectionDiffusion_device execution failed");
}

void calculateMacrosopicMean(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned int* geoD, unsigned int* neighborX,
                    unsigned int* neighborY, unsigned int* neighborZ, unsigned int tdiff, unsigned long long numberOfLBnodes,
                    unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    calculateMacrosopicMean_device<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, geoD, neighborX, neighborY, neighborZ, tdiff,
                                                  numberOfLBnodes, isEvenTimestep);
    getLastCudaError("calculateMacrosopicMean_device execution failed");
}

void ResetMeanValuesSP27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, unsigned long long numberOfLBnodes,
                         unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBResetMeanValuesSP27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, numberOfLBnodes, isEvenTimestep);
    getLastCudaError("LBResetMeanValuesSP27 execution failed");
}
//////////////////////////////////////////////////////////////////////////
void ResetMeanValuesAD27(real* vxD, real* vyD, real* vzD, real* rhoD, real* pressD, real* concD,
                         unsigned long long numberOfLBnodes, unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfLBnodes);

    LBResetMeanValuesAD27<<<grid.grid, grid.threads>>>(vxD, vyD, vzD, rhoD, pressD, concD, numberOfLBnodes, isEvenTimestep);
    getLastCudaError("LBResetMeanValuesAD27 execution failed");
}

void calculateMeasurePoints(real* vxMP, real* vyMP, real* vzMP, real* rhoMP, unsigned int* kMP,
                           unsigned int numberOfPointskMP, unsigned int MPClockCycle, unsigned int t, unsigned int* geoD,
                           unsigned int* neighborX, unsigned int* neighborY, unsigned int* neighborZ,
                           unsigned long long numberOfLBnodes, real* DD, unsigned int numberOfThreads, bool isEvenTimestep)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(numberOfThreads, numberOfPointskMP);

    LBCalcMeasurePoints<<<grid.grid, grid.threads>>>(vxMP, vyMP, vzMP, rhoMP, kMP, numberOfPointskMP, MPClockCycle, t, geoD,
                                                     neighborX, neighborY, neighborZ, numberOfLBnodes, DD, isEvenTimestep);
    getLastCudaError("LBCalcMeasurePoints execution failed");
}

}

//! \}
