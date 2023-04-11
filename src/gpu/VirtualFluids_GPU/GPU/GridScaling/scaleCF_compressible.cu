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

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

//////////////////////////////////////////////////////////////////////////
//! \brief Interpolate from coarse to fine nodes
//! \details This scaling function is designed for the Cumulant K17 Kernel chimera collision kernel.
//!
//! The function is executed in the following steps:
//!
// based on scaleCF_RhoSq_comp_27
template<bool hasTurbulentViscosity> __global__ void scaleCF_compressible(
    real* distributionsCoarse, 
    real* distributionsFine, 
    unsigned int* neighborXcoarse,
    unsigned int* neighborYcoarse,
    unsigned int* neighborZcoarse,
    unsigned int* neighborXfine,
    unsigned int* neighborYfine,
    unsigned int* neighborZfine,
    unsigned long long numberOfLBnodesCoarse, 
    unsigned long long numberOfLBnodesFine, 
    bool isEvenTimestep,
    unsigned int* indicesCoarseMMM, 
    unsigned int* indicesFineMMM, 
    unsigned int numberOfInterfaceNodes, 
    real omegaCoarse, 
    real omegaFine, 
    real* turbulentViscosityCoarse,
    real* turbulentViscosityFine,
    ICellNeigh neighborCoarseToFine)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get the node index coordinates from threadId_100, blockId_100, blockDim and gridDim.
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
    real eps_new = c1o2; // ratio of grid resolutions
    real omegaC  = omegaCoarse;
    real omegaF  = omegaFine;

    vf::lbm::MomentsOnSourceNode moments_PPP;
    vf::lbm::MomentsOnSourceNode moments_MPP;
    vf::lbm::MomentsOnSourceNode moments_PMP;
    vf::lbm::MomentsOnSourceNode moments_MMP;
    vf::lbm::MomentsOnSourceNode moments_PPM;
    vf::lbm::MomentsOnSourceNode moments_MPM;
    vf::lbm::MomentsOnSourceNode moments_PMM;
    vf::lbm::MomentsOnSourceNode moments_MMM;

    ////////////////////////////////////////////////////////////////////////////////
    //! - Calculate moments for each source node 
    //!
    ////////////////////////////////////////////////////////////////////////////////
    // source node BSW = MMM
    ////////////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    unsigned int k_base_000 = indicesCoarseMMM[nodeIndex];
    unsigned int k_base_M00 = neighborXcoarse [k_base_000];
    unsigned int k_base_0M0 = neighborYcoarse [k_base_000];
    unsigned int k_base_00M = neighborZcoarse [k_base_000];
    unsigned int k_base_MM0 = neighborYcoarse [k_base_M00];
    unsigned int k_base_M0M = neighborZcoarse [k_base_M00];
    unsigned int k_base_0MM = neighborZcoarse [k_base_0M0];
    unsigned int k_base_MMM = neighborZcoarse [k_base_MM0];
    ////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    unsigned int k_000 = k_base_000;
    unsigned int k_M00 = k_base_M00;
    unsigned int k_0M0 = k_base_0M0;
    unsigned int k_00M = k_base_00M;
    unsigned int k_MM0 = k_base_MM0;
    unsigned int k_M0M = k_base_M0M;
    unsigned int k_0MM = k_base_0MM;
    unsigned int k_MMM = k_base_MMM;

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);

    vf::lbm::Distribution27 distribution;

    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_MMM);

    //////////////////////////////////////////////////////////////////////////
    // source node TSW = MMP
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices - has to be recalculated for the new source node
    k_000 = k_00M;
    k_M00 = k_M0M;
    k_0M0 = k_0MM;
    k_00M = neighborZcoarse[k_00M];
    k_MM0 = k_MMM;
    k_M0M = neighborZcoarse[k_M0M];
    k_0MM = neighborZcoarse[k_0MM];
    k_MMM = neighborZcoarse[k_MMM];

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);

    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_MMP);

    //////////////////////////////////////////////////////////////////////////
    // source node TSE = PMP
    //////////////////////////////////////////////////////////////////////////
    // index
    k_000 = k_M00;
    k_M00 = neighborXcoarse[k_M00];
    k_0M0 = k_MM0;
    k_00M = k_M0M;
    k_MM0 = neighborXcoarse[k_MM0];
    k_M0M = neighborXcoarse[k_M0M];
    k_0MM = k_MMM;
    k_MMM = neighborXcoarse[k_MMM];

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);

    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_PMP);

    //////////////////////////////////////////////////////////////////////////
    // source node BSE = PMM 
    //////////////////////////////////////////////////////////////////////////
    // index
    k_00M = k_000;
    k_M0M = k_M00;
    k_0MM = k_0M0;
    k_MMM = k_MM0;
    k_000 = k_base_M00;
    k_M00 = neighborXcoarse[k_base_M00];
    k_0M0 = k_base_MM0;
    k_MM0 = neighborXcoarse[k_base_MM0];

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);

    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_PMM);

    //////////////////////////////////////////////////////////////////////////
    // source node BNW = MPM
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors --> indices of all source nodes
    k_base_000 = k_base_0M0;
    k_base_M00 = k_base_MM0;
    k_base_0M0 = neighborYcoarse[k_base_0M0];
    k_base_00M = k_base_0MM;
    k_base_MM0 = neighborYcoarse[k_base_MM0];
    k_base_M0M = k_base_MMM;
    k_base_0MM = neighborYcoarse[k_base_0MM];
    k_base_MMM = neighborYcoarse[k_base_MMM];
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

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);

    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_MPM);

    //////////////////////////////////////////////////////////////////////////
    // source node TNW = MPP
    //////////////////////////////////////////////////////////////////////////
    // index
    k_000 = k_00M;
    k_M00 = k_M0M;
    k_0M0 = k_0MM;
    k_00M = neighborZcoarse[k_00M];
    k_MM0 = k_MMM;
    k_M0M = neighborZcoarse[k_M0M];
    k_0MM = neighborZcoarse[k_0MM];
    k_MMM = neighborZcoarse[k_MMM];

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);
    
    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_MPP);
    //////////////////////////////////////////////////////////////////////////
    // source node TNE = PPP
    //////////////////////////////////////////////////////////////////////////
    // index
    // index
    k_000 = k_M00;
    k_M00 = neighborXcoarse[k_M00];
    k_0M0 = k_MM0;
    k_00M = k_M0M;
    k_MM0 = neighborXcoarse[k_MM0];
    k_M0M = neighborXcoarse[k_M0M];
    k_0MM = k_MMM;
    k_MMM = neighborXcoarse[k_MMM];

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);

    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_PPP);
    //////////////////////////////////////////////////////////////////////////
    // source node BNE = PPM
    //////////////////////////////////////////////////////////////////////////
    // index
    k_00M = k_000;
    k_M0M = k_M00;
    k_0MM = k_0M0;
    k_MMM = k_MM0;
    k_000 = k_base_M00;
    k_M00 = neighborXcoarse[k_base_M00];
    k_0M0 = k_base_MM0;
    k_MM0 = neighborXcoarse[k_base_MM0];

    if(hasTurbulentViscosity) omegaC = omegaCoarse / (c1o1 + c3o1*omegaCoarse*turbulentViscosityCoarse[k_000]);

    readDistributionFromList(distribution, distCoarse, k_000, k_M00, k_0M0, k_00M, k_MM0, k_M0M, k_0MM, k_MMM);
    vf::lbm::calculateMomentsOnSourceNodes(distribution.f, omegaC, moments_PPM);


    ////////////////////////////////////////////////////////////////////////////////
    //! - Set the relative position of the offset cell {-1, 0, 1}
    //!
    const real xoff    = neighborCoarseToFine.x[nodeIndex];
    const real yoff    = neighborCoarseToFine.y[nodeIndex];
    const real zoff    = neighborCoarseToFine.z[nodeIndex];

    vf::lbm::Coefficients coefficients;
    vf::lbm::calculateCoefficients(xoff, yoff, zoff, coefficients, 
        moments_PPP,
        moments_MPP,
        moments_PMP,
        moments_MMP,
        moments_PPM,
        moments_MPM,
        moments_PMM,
        moments_MMM);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position BSW = MMM: -0.25, -0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real x = -c1o4;
    real y = -c1o4;
    real z = -c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    k_base_000 = indicesFineMMM[nodeIndex];
    k_base_M00 = neighborXfine [k_base_000];
    k_base_0M0 = neighborYfine [k_base_000];
    k_base_00M = neighborZfine [k_base_000];
    k_base_MM0 = neighborYfine [k_base_M00];
    k_base_M0M = neighborZfine [k_base_M00];
    k_base_0MM = neighborZfine [k_base_0M0];
    k_base_MMM = neighborZfine [k_base_MM0];
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_000 = k_base_000;
    k_M00 = k_base_M00;
    k_0M0 = k_base_0M0;
    k_00M = k_base_00M;
    k_MM0 = k_base_MM0;
    k_M0M = k_base_M0M;
    k_0MM = k_base_0MM;
    k_MMM = k_base_MMM;
    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (zeroth to sixth order) on destination node
    //!

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    real f[27];
    vf::lbm::interpolate_cf(
        x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
    //////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TSW = MMP: -0.25, -0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = -c1o4;
    y = -c1o4;
    z =  c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_000 = k_00M;
    k_M00 = k_M0M;
    k_0M0 = k_0MM;
    k_00M = neighborZfine[k_00M];
    k_MM0 = k_MMM;
    k_M0M = neighborZfine[k_M0M];
    k_0MM = neighborZfine[k_0MM];
    k_MMM = neighborZfine[k_MMM];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    vf::lbm::interpolate_cf(
        x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TSE = PMP: 0.25, -0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x =  c1o4;
    y = -c1o4;
    z =  c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_000 = k_M00;
    k_M00 = neighborXfine[k_M00];
    k_0M0 = k_MM0;
    k_00M = k_M0M;
    k_MM0 = neighborXfine[k_MM0];
    k_M0M = neighborXfine[k_M0M];
    k_0MM = k_MMM;
    k_MMM = neighborXfine[k_MMM];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    vf::lbm::interpolate_cf(
        x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position BSE = PMM: 0.25, -0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x =  c1o4;
    y = -c1o4;
    z = -c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_00M = k_000;
    k_M0M = k_M00;
    k_0MM = k_0M0;
    k_MMM = k_MM0;
    k_000 = k_base_M00;
    k_M00 = neighborXfine[k_base_M00];
    k_0M0 = k_base_MM0;
    k_MM0 = neighborXfine[k_base_MM0];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);


    // real f[27];
    vf::lbm::interpolate_cf(
        x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position BNW = MPM: -0.25, 0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = -c1o4;
    y =  c1o4;
    z = -c1o4;
    
    //////////////////////////////////////////////////////////////////////////
    // index of the base node and its neighbors
    k_base_000 = k_base_0M0;
    k_base_M00 = k_base_MM0;
    k_base_0M0 = neighborYfine[k_base_0M0];
    k_base_00M = k_base_0MM;
    k_base_MM0 = neighborYfine[k_base_MM0];
    k_base_M0M = k_base_MMM;
    k_base_0MM = neighborYfine[k_base_0MM];
    k_base_MMM = neighborYfine[k_base_MMM];

    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_000 = k_base_000;
    k_M00 = k_base_M00;
    k_0M0 = k_base_0M0;
    k_00M = k_base_00M;
    k_MM0 = k_base_MM0;
    k_M0M = k_base_M0M;
    k_0MM = k_base_0MM;
    k_MMM = k_base_MMM;

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    vf::lbm::interpolate_cf(
        x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TNW = MPP: -0.25, 0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = -c1o4;
    y =  c1o4;
    z =  c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_000 = k_00M;
    k_M00 = k_M0M;
    k_0M0 = k_0MM;
    k_00M = neighborZfine[k_00M];
    k_MM0 = k_MMM;
    k_M0M = neighborZfine[k_M0M];
    k_0MM = neighborZfine[k_0MM];
    k_MMM = neighborZfine[k_MMM];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    vf::lbm::interpolate_cf(
        x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Position TNE = PPP: 0.25, 0.25, 0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x = c1o4;
    y = c1o4;
    z = c1o4;
    ////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_000 = k_M00;
    k_M00 = neighborXfine[k_M00];
    k_0M0 = k_MM0;
    k_00M = k_M0M;
    k_MM0 = neighborXfine[k_MM0];
    k_M0M = neighborXfine[k_M0M];
    k_0MM = k_MMM;
    k_MMM = neighborXfine[k_MMM];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);

    vf::lbm::interpolate_cf(
        x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //Position BNE = PPM: 0.25, 0.25, -0.25
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x =  c1o4;
    y =  c1o4;
    z = -c1o4;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set neighbor indices
    k_00M = k_000;
    k_M0M = k_M00;
    k_0MM = k_0M0;
    k_MMM = k_MM0;
    k_000 = k_base_M00;
    k_M00 = neighborXfine[k_base_M00];
    k_0M0 = k_base_MM0;
    k_MM0 = neighborXfine[k_base_MM0];

    if(hasTurbulentViscosity) omegaF = omegaFine/ (c1o1 + c3o1*omegaFine*turbulentViscosityFine[k_000]);
    vf::lbm::interpolate_cf(x, y, z, f, coefficients, eps_new, omegaF
    );

    //////////////////////////////////////////////////////////////////////////
    //! - Write distributions: style of reading and writing the distributions from/to
    //! stored arrays dependent on timestep is based on the esoteric twist algorithm
    //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017),
    //! DOI:10.3390/computation5020019 ]</b></a>
    //!
    (distFine.f[DIR_000])[k_000] = f[DIR_000];
    (distFine.f[DIR_P00])[k_000] = f[DIR_P00];
    (distFine.f[DIR_M00])[k_M00] = f[DIR_M00];
    (distFine.f[DIR_0P0])[k_000] = f[DIR_0P0];
    (distFine.f[DIR_0M0])[k_0M0] = f[DIR_0M0];
    (distFine.f[DIR_00P])[k_000] = f[DIR_00P];
    (distFine.f[DIR_00M])[k_00M] = f[DIR_00M];
    (distFine.f[DIR_PP0])[k_000] = f[DIR_PP0];
    (distFine.f[DIR_MM0])[k_MM0] = f[DIR_MM0];
    (distFine.f[DIR_PM0])[k_0M0] = f[DIR_PM0];
    (distFine.f[DIR_MP0])[k_M00] = f[DIR_MP0];
    (distFine.f[DIR_P0P])[k_000] = f[DIR_P0P];
    (distFine.f[DIR_M0M])[k_M0M] = f[DIR_M0M];
    (distFine.f[DIR_P0M])[k_00M] = f[DIR_P0M];
    (distFine.f[DIR_M0P])[k_M00] = f[DIR_M0P];
    (distFine.f[DIR_0PP])[k_000] = f[DIR_0PP];
    (distFine.f[DIR_0MM])[k_0MM] = f[DIR_0MM];
    (distFine.f[DIR_0PM])[k_00M] = f[DIR_0PM];
    (distFine.f[DIR_0MP])[k_0M0] = f[DIR_0MP];
    (distFine.f[DIR_PPP])[k_000] = f[DIR_PPP];
    (distFine.f[DIR_MPP])[k_M00] = f[DIR_MPP];
    (distFine.f[DIR_PMP])[k_0M0] = f[DIR_PMP];
    (distFine.f[DIR_MMP])[k_MM0] = f[DIR_MMP];
    (distFine.f[DIR_PPM])[k_00M] = f[DIR_PPM];
    (distFine.f[DIR_MPM])[k_M0M] = f[DIR_MPM];
    (distFine.f[DIR_PMM])[k_0MM] = f[DIR_PMM];
    (distFine.f[DIR_MMM])[k_MMM] = f[DIR_MMM];
}

template __global__ void scaleCF_compressible<true>( real* distributionsCoarse, real* distributionsFine, unsigned int* neighborXcoarse, unsigned int* neighborYcoarse, unsigned int* neighborZcoarse, unsigned int* neighborXfine, unsigned int* neighborYfine, unsigned int* neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int* indicesCoarseMMM, unsigned int* indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh offsetCF);

template __global__ void scaleCF_compressible<false>( real* distributionsCoarse, real* distributionsFine, unsigned int* neighborXcoarse, unsigned int* neighborYcoarse, unsigned int* neighborZcoarse, unsigned int* neighborXfine, unsigned int* neighborYfine, unsigned int* neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int* indicesCoarseMMM, unsigned int* indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh offsetCF);