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
//! \file ScalingUtilities.h
//! \ingroup LBM/GPUHelperFunctions
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================
#ifndef SCALING_HELPER_FUNCTIONS_H
#define SCALING_HELPER_FUNCTIONS_H

#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "basics/constants/NumericConstants.h"

using namespace vf::basics::constant;
#include <lbm/KernelParameter.h>

using namespace vf::lbm::dir;

namespace vf::gpu
{



template<bool hasTurbulentViscosity> __device__ void calculate_moment_set(
    vf::lbm::MomentsOnSourceNodeSet& moments_set,
    const unsigned nodeIndex,
    real *distributionsFine,
    unsigned int *neighborXfine,
    unsigned int *neighborYfine,
    unsigned int *neighborZfine,
    unsigned int *indicesFineMMM,
    real* turbulentViscosityFine,
    unsigned long long numberOfLBnodesFine,
    const real omegaFine,
    bool isEvenTimestep
)
{
    real omegaF = omegaFine;
    Distributions27 distFine;
    getPointersToDistributions(distFine, distributionsFine, numberOfLBnodesFine, isEvenTimestep);

    vf::gpu::ListIndices indices;

    //////////////////////////////////////////////////////////////////////////
    //! - Calculate moments for each source node 
    //!
    //////////////////////////////////////////////////////////////////////////
    // source node BSW = MMM
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    unsigned int k_base_000 = indicesFineMMM[nodeIndex];
    unsigned int k_base_M00 = neighborXfine [k_base_000];
    unsigned int k_base_0M0 = neighborYfine [k_base_000];
    unsigned int k_base_00M = neighborZfine [k_base_000];
    unsigned int k_base_MM0 = neighborYfine [k_base_M00];
    unsigned int k_base_M0M = neighborZfine [k_base_M00];
    unsigned int k_base_0MM = neighborZfine [k_base_0M0];
    unsigned int k_base_MMM = neighborZfine [k_base_MM0];
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    real f_fine[27];

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_MMM);

    //////////////////////////////////////////////////////////////////////////
    // source node TSW = MMP
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices - has to be recalculated for the new source node
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZfine[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZfine[indices.k_M0M];
    indices.k_0MM = neighborZfine[indices.k_0MM];
    indices.k_MMM = neighborZfine[indices.k_MMM];

    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_MMP);

    //////////////////////////////////////////////////////////////////////////
    // source node TSE = PMP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXfine[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXfine[indices.k_MM0];
    indices.k_M0M = neighborXfine[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXfine[indices.k_MMM];

    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_PMP);

    //////////////////////////////////////////////////////////////////////////
    // source node BSE = PMM 
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXfine[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXfine[k_base_MM0];

    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_PMM);

    //////////////////////////////////////////////////////////////////////////
    // source node BNW = MPM
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors --> indices of all source nodes
    k_base_000 = k_base_0M0;
    k_base_M00 = k_base_MM0;
    k_base_0M0 = neighborYfine[k_base_0M0];
    k_base_00M = k_base_0MM;
    k_base_MM0 = neighborYfine[k_base_MM0];
    k_base_M0M = k_base_MMM;
    k_base_0MM = neighborYfine[k_base_0MM];
    k_base_MMM = neighborYfine[k_base_MMM];
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_MPM);

    //////////////////////////////////////////////////////////////////////////
    // source node TNW = MPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZfine[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZfine[indices.k_M0M];
    indices.k_0MM = neighborZfine[indices.k_0MM];
    indices.k_MMM = neighborZfine[indices.k_MMM];

    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_MPP);

    //////////////////////////////////////////////////////////////////////////
    // source node TNE = PPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXfine[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXfine[indices.k_MM0];
    indices.k_M0M = neighborXfine[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXfine[indices.k_MMM];

    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_PPP);

    //////////////////////////////////////////////////////////////////////////
    // source node BNE = PPM
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXfine[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXfine[k_base_MM0];
    
    omegaF = hasTurbulentViscosity ? calculateOmega(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    readDistributionFromList(f_fine, distFine, indices);
    vf::lbm::calculateMomentsOnSourceNodes(f_fine, omegaF, moments_set.moments_PPM);
}


__device__ __inline__ void readDistributionFromList(vf::lbm::Distribution27 &distribution, const Distributions27 &dist, unsigned int &k_000,
                                                         unsigned int &k_M00, unsigned int &k_0M0, unsigned int &k_00M,
                                                         unsigned int &k_MM0, unsigned int &k_M0M, unsigned int &k_0MM,
                                                         unsigned int &k_MMM)
{
    distribution.f[DIR_000] = (dist.f[DIR_000])[k_000];
    distribution.f[DIR_P00] = (dist.f[DIR_P00])[k_000];
    distribution.f[DIR_M00] = (dist.f[DIR_M00])[k_M00];
    distribution.f[DIR_0P0] = (dist.f[DIR_0P0])[k_000];
    distribution.f[DIR_0M0] = (dist.f[DIR_0M0])[k_0M0];
    distribution.f[DIR_00P] = (dist.f[DIR_00P])[k_000];
    distribution.f[DIR_00M] = (dist.f[DIR_00M])[k_00M];
    distribution.f[DIR_PP0] = (dist.f[DIR_PP0])[k_000];
    distribution.f[DIR_MM0] = (dist.f[DIR_MM0])[k_MM0];
    distribution.f[DIR_PM0] = (dist.f[DIR_PM0])[k_0M0];
    distribution.f[DIR_MP0] = (dist.f[DIR_MP0])[k_M00];
    distribution.f[DIR_P0P] = (dist.f[DIR_P0P])[k_000];
    distribution.f[DIR_M0M] = (dist.f[DIR_M0M])[k_M0M];
    distribution.f[DIR_P0M] = (dist.f[DIR_P0M])[k_00M];
    distribution.f[DIR_M0P] = (dist.f[DIR_M0P])[k_M00];
    distribution.f[DIR_0PP] = (dist.f[DIR_0PP])[k_000];
    distribution.f[DIR_0MM] = (dist.f[DIR_0MM])[k_0MM];
    distribution.f[DIR_0PM] = (dist.f[DIR_0PM])[k_00M];
    distribution.f[DIR_0MP] = (dist.f[DIR_0MP])[k_0M0];
    distribution.f[DIR_PPP] = (dist.f[DIR_PPP])[k_000];
    distribution.f[DIR_MPP] = (dist.f[DIR_MPP])[k_M00];
    distribution.f[DIR_PMP] = (dist.f[DIR_PMP])[k_0M0];
    distribution.f[DIR_MMP] = (dist.f[DIR_MMP])[k_MM0];
    distribution.f[DIR_PPM] = (dist.f[DIR_PPM])[k_00M];
    distribution.f[DIR_MPM] = (dist.f[DIR_MPM])[k_M0M];
    distribution.f[DIR_PMM] = (dist.f[DIR_PMM])[k_0MM];
    distribution.f[DIR_MMM] = (dist.f[DIR_MMM])[k_MMM];
}

} // namespace vf::gpu

#endif
