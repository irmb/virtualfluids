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
//! \file scaleFC_compressible.cu
//! \ingroup GPU/GridScaling
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================

#include "LBM/GPUHelperFunctions/ChimeraTransformation.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"
#include "LBM/GPUHelperFunctions/ScalingUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

//////////////////////////////////////////////////////////////////////////
//! \brief Interpolate from fine to coarse
//! \details This scaling function is designed for the Cumulant K17 Kernel chimera collision kernel
//! The function is executed in the following steps:
//!

// based on scaleFC_RhoSq_comp_27
template<bool hasTurbulentViscosity> __global__ void scaleFC_compressible(
    real *distributionsCoarse,
    real *distributionsFine,
    unsigned int *neighborXcoarse,
    unsigned int *neighborYcoarse,
    unsigned int *neighborZcoarse,
    unsigned int *neighborXfine,
    unsigned int *neighborYfine,
    unsigned int *neighborZfine,
    unsigned long long numberOfLBnodesCoarse,
    unsigned long long numberOfLBnodesFine,
    bool isEvenTimestep,
    unsigned int *indicesCoarse000,
    unsigned int *indicesFineMMM,
    unsigned int numberOfInterfaceNodes,
    real omegaCoarse,
    real omegaFine,
    real* turbulentViscosityCoarse,
    real* turbulentViscosityFine,
    ICellNeigh neighborFineToCoarse)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get the node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    //////////////////////////////////////////////////////////////////////////
    //! - Return for non-interface node
    if (nodeIndex >= numberOfInterfaceNodes)
        return;

    //////////////////////////////////////////////////////////////////////////
    //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on
    //! timestep is based on the esoteric twist algorithm \ref <a
    //! href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    Distributions27 distFine, distCoarse;
    getPointersToDistributions(distFine, distributionsFine, numberOfLBnodesFine, true);
    getPointersToDistributions(distCoarse, distributionsCoarse, numberOfLBnodesCoarse, isEvenTimestep);

    ////////////////////////////////////////////////////////////////////////////////
    //! - declare local variables for source nodes
    //!
    real eps_new = c2o1; // ratio of grid resolutions
    real omegaF  = omegaFine;
    real omegaC  = omegaCoarse;

    // // zeroth and first order moments at the source nodes
    // real drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP;
    // real drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP;
    // real drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP;
    // real drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP;
    // real drho_PPM, vx1_PPM, vx2_PPM, vx3_PPM;
    // real drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM;
    // real drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM;
    // real drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM;

    // // second order moments at the source nodes
    // real kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP;
    // real kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP;
    // real kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP;
    // real kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP;
    // real kxyFromfcNEQ_PPM, kyzFromfcNEQ_PPM, kxzFromfcNEQ_PPM, kxxMyyFromfcNEQ_PPM, kxxMzzFromfcNEQ_PPM;
    // real kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM;
    // real kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM;
    // real kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM;

    vf::lbm::MomentsOnSourceNode moments_PPP;
    vf::lbm::MomentsOnSourceNode moments_MPP;
    vf::lbm::MomentsOnSourceNode moments_PMP;
    vf::lbm::MomentsOnSourceNode moments_MMP;
    vf::lbm::MomentsOnSourceNode moments_PPM;
    vf::lbm::MomentsOnSourceNode moments_MPM;
    vf::lbm::MomentsOnSourceNode moments_PMM;
    vf::lbm::MomentsOnSourceNode moments_MMM;

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
    unsigned int k_000 = k_base_000;
    unsigned int k_M00 = k_base_M00;
    unsigned int k_0M0 = k_base_0M0;
    unsigned int k_00M = k_base_00M;
    unsigned int k_MM0 = k_base_MM0;
    unsigned int k_M0M = k_base_M0M;
    unsigned int k_0MM = k_base_0MM;
    unsigned int k_MMM = k_base_MMM;

    if(hasTurbulentViscosity) omegaF = omegaFine / (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    vf::lbm::Distribution27 distribution;
    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_MMM);

    // calculateMomentsOnSourceNodes( distFine, omegaF,
    //     k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MMM, vx1_MMM, vx2_MMM, vx3_MMM,
    //     kxyFromfcNEQ_MMM, kyzFromfcNEQ_MMM, kxzFromfcNEQ_MMM, kxxMyyFromfcNEQ_MMM, kxxMzzFromfcNEQ_MMM);


    //////////////////////////////////////////////////////////////////////////
    // source node TSW = MMP
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices - has to be recalculated for the new source node
    k_000 = k_00M;
    k_M00 = k_M0M;
    k_0M0 = k_0MM;
    k_00M = neighborZfine[k_00M];
    k_MM0 = k_MMM;
    k_M0M = neighborZfine[k_M0M];
    k_0MM = neighborZfine[k_0MM];
    k_MMM = neighborZfine[k_MMM];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    
    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_MMP);

    // calculateMomentsOnSourceNodes( distFine, omegaF,
    //     k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MMP, vx1_MMP, vx2_MMP, vx3_MMP,
    //     kxyFromfcNEQ_MMP, kyzFromfcNEQ_MMP, kxzFromfcNEQ_MMP, kxxMyyFromfcNEQ_MMP, kxxMzzFromfcNEQ_MMP);

    //////////////////////////////////////////////////////////////////////////
    // source node TSE = PMP
    //////////////////////////////////////////////////////////////////////////
    // index
    k_000 = k_M00;
    k_M00 = neighborXfine[k_M00];
    k_0M0 = k_MM0;
    k_00M = k_M0M;
    k_MM0 = neighborXfine[k_MM0];
    k_M0M = neighborXfine[k_M0M];
    k_0MM = k_MMM;
    k_MMM = neighborXfine[k_MMM];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_PMP);

    // calculateMomentsOnSourceNodes( distFine, omegaF,
    //     k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_PMP, vx1_PMP, vx2_PMP, vx3_PMP,
    //     kxyFromfcNEQ_PMP, kyzFromfcNEQ_PMP, kxzFromfcNEQ_PMP, kxxMyyFromfcNEQ_PMP, kxxMzzFromfcNEQ_PMP);

    //////////////////////////////////////////////////////////////////////////
    // source node BSE = PMM 
    //////////////////////////////////////////////////////////////////////////
    // index
    k_00M = k_000;
    k_M0M = k_M00;
    k_0MM = k_0M0;
    k_MMM = k_MM0;
    k_000 = k_base_M00;
    k_M00 = neighborXfine[k_base_M00];
    k_0M0 = k_base_MM0;
    k_MM0 = neighborXfine[k_base_MM0];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_PMM);

    // calculateMomentsOnSourceNodes( distFine, omegaF,
    //     k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_PMM, vx1_PMM, vx2_PMM, vx3_PMM,
    //     kxyFromfcNEQ_PMM, kyzFromfcNEQ_PMM, kxzFromfcNEQ_PMM, kxxMyyFromfcNEQ_PMM, kxxMzzFromfcNEQ_PMM);

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
    k_000 = k_base_000;
    k_M00 = k_base_M00;
    k_0M0 = k_base_0M0;
    k_00M = k_base_00M;
    k_MM0 = k_base_MM0;
    k_M0M = k_base_M0M;
    k_0MM = k_base_0MM;
    k_MMM = k_base_MMM;

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_MPM);

    // calculateMomentsOnSourceNodes( distFine, omegaF,
    //     k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MPM, vx1_MPM, vx2_MPM, vx3_MPM,
    //     kxyFromfcNEQ_MPM, kyzFromfcNEQ_MPM, kxzFromfcNEQ_MPM, kxxMyyFromfcNEQ_MPM, kxxMzzFromfcNEQ_MPM);

    //////////////////////////////////////////////////////////////////////////
    // source node TNW = MPP
    //////////////////////////////////////////////////////////////////////////
    // index
    k_000 = k_00M;
    k_M00 = k_M0M;
    k_0M0 = k_0MM;
    k_00M = neighborZfine[k_00M];
    k_MM0 = k_MMM;
    k_M0M = neighborZfine[k_M0M];
    k_0MM = neighborZfine[k_0MM];
    k_MMM = neighborZfine[k_MMM];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_MPP);
    
    // calculateMomentsOnSourceNodes( distFine, omegaF,
    //     k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_MPP, vx1_MPP, vx2_MPP, vx3_MPP,
    //     kxyFromfcNEQ_MPP, kyzFromfcNEQ_MPP, kxzFromfcNEQ_MPP, kxxMyyFromfcNEQ_MPP, kxxMzzFromfcNEQ_MPP);

    //////////////////////////////////////////////////////////////////////////
    // source node TNE = PPP
    //////////////////////////////////////////////////////////////////////////
    // index
    k_000 = k_M00;
    k_M00 = neighborXfine[k_M00];
    k_0M0 = k_MM0;
    k_00M = k_M0M;
    k_MM0 = neighborXfine[k_MM0];
    k_M0M = neighborXfine[k_M0M];
    k_0MM = k_MMM;
    k_MMM = neighborXfine[k_MMM];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);


    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_PPP);

    // calculateMomentsOnSourceNodes( distFine, omegaF,
    //     k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM, drho_PPP, vx1_PPP, vx2_PPP, vx3_PPP,
    //     kxyFromfcNEQ_PPP, kyzFromfcNEQ_PPP, kxzFromfcNEQ_PPP, kxxMyyFromfcNEQ_PPP, kxxMzzFromfcNEQ_PPP);

    //////////////////////////////////////////////////////////////////////////
    // source node BNE = PPM
    //////////////////////////////////////////////////////////////////////////
    // index
    k_00M = k_000;
    k_M0M = k_M00;
    k_0MM = k_0M0;
    k_MMM = k_MM0;
    k_000 = k_base_M00;
    k_M00 = neighborXfine[k_base_M00];
    k_0M0 = k_base_MM0;
    k_MM0 = neighborXfine[k_base_MM0];
    
    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    readDistributionFromList(distribution, distFine, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaF, moments_PPM);

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set the relative position of the offset cell {-1, 0, 1}
    //!
    const real xoff    = neighborFineToCoarse.x[nodeIndex];
    const real yoff    = neighborFineToCoarse.y[nodeIndex];
    const real zoff    = neighborFineToCoarse.z[nodeIndex];
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position Coarse 0., 0., 0.
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // x = 0.;
    // y = 0.;
    // z = 0.;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // index of the destination node and its neighbors
    k_000 = indicesCoarse000[nodeIndex];
    k_M00 = neighborXcoarse [k_000];
    k_0M0 = neighborYcoarse [k_000];
    k_00M = neighborZcoarse [k_000];
    k_MM0 = neighborYcoarse [k_M00];
    k_M0M = neighborZcoarse [k_M00];
    k_0MM = neighborZcoarse [k_0M0];
    k_MMM = neighborZcoarse [k_MM0];
    ////////////////////////////////////////////////////////////////////////////////////
    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse * turbulentViscosityCoarse[k_000]);

    vf::lbm::Coefficients coefficients;
    calculateCoefficients(xoff, yoff, zoff, coefficients, 
        moments_PPP,
        moments_MPP,
        moments_PMP,
        moments_MMP,
        moments_PPM,
        moments_MPM,
        moments_PMM,
        moments_MMM);
    
    real f[27];
    vf::lbm::interpolate_fc(f, eps_new, omegaC, coefficients);

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distCoarse.f[DIR_000])[k_000] = f[DIR_000];
    (distCoarse.f[DIR_P00])[k_000] = f[DIR_P00];
    (distCoarse.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distCoarse.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distCoarse.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distCoarse.f[DIR_00P])[k_000] = f[DIR_00P];
    (distCoarse.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distCoarse.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distCoarse.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distCoarse.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distCoarse.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distCoarse.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distCoarse.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distCoarse.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distCoarse.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distCoarse.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distCoarse.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distCoarse.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distCoarse.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distCoarse.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distCoarse.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distCoarse.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distCoarse.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distCoarse.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distCoarse.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distCoarse.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distCoarse.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
    ////////////////////////////////////////////////////////////////////////////////////
}

template __global__ void scaleFC_compressible<true>( real *distributionsCoarse, real *distributionsFine, unsigned int *neighborXcoarse, unsigned int *neighborYcoarse, unsigned int *neighborZcoarse, unsigned int *neighborXfine, unsigned int *neighborYfine, unsigned int *neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int *indicesCoarse000, unsigned int *indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh neighborFineToCoarse);

template __global__ void scaleFC_compressible<false>( real *distributionsCoarse, real *distributionsFine, unsigned int *neighborXcoarse, unsigned int *neighborYcoarse, unsigned int *neighborZcoarse, unsigned int *neighborXfine, unsigned int *neighborYfine, unsigned int *neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int *indicesCoarse000, unsigned int *indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh neighborFineToCoarse);