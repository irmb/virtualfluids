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
//! \addtogroup gpu_GridScaling GridScaling
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================

#include <basics/DataTypes.h>
#include "Utilities/KernelUtilities.h"
#include "Utilities/ScalingUtilities.h"

#include <lbm/refinement/InterpolationCF.h>
#include <lbm/interpolation/InterpolationCoefficients.h>
#include <lbm/collision/TurbulentViscosity.h>



template <bool hasTurbulentViscosity> __device__ void interpolate(
    vf::lbm::InterpolationCoefficients& coefficients,
    const unsigned int nodeIndex,
    real* distributionsFine, 
    unsigned int* neighborXfine,
    unsigned int* neighborYfine,
    unsigned int* neighborZfine,
    unsigned long long numberOfLBnodesFine,
    unsigned int* indicesFineMMM,
    real omegaFine,
    real* turbulentViscosityFine)
{
    Distributions27 distFine;
    vf::gpu::getPointersToDistributions(distFine, distributionsFine, numberOfLBnodesFine, true);

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
    vf::gpu::ListIndices indices;
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;
    ////////////////////////////////////////////////////////////////////////////////
    //! - Set moments (zeroth to sixth order) on destination node
    //!
    real omegaF  = omegaFine;
    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    const real epsilon_new = c1o2; // ratio of grid resolutions
    real f[27];
    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);


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
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZfine[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZfine[indices.k_M0M];
    indices.k_0MM = neighborZfine[indices.k_0MM];
    indices.k_MMM = neighborZfine[indices.k_MMM];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);

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
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXfine[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXfine[indices.k_MM0];
    indices.k_M0M = neighborXfine[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXfine[indices.k_MMM];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);

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
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXfine[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXfine[k_base_MM0];

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);
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
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    ////////////////////////////////////////////////////////////////////////////////
    // Set moments (zeroth to sixth orders) on destination node

    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);
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
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZfine[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZfine[indices.k_M0M];
    indices.k_0MM = neighborZfine[indices.k_0MM];
    indices.k_MMM = neighborZfine[indices.k_MMM];

    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);
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
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXfine[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXfine[indices.k_MM0];
    indices.k_M0M = neighborXfine[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXfine[indices.k_MMM];

    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);
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
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXfine[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXfine[k_base_MM0];

    omegaF = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaFine, turbulentViscosityFine[indices.k_000]) : omegaFine;

    vf::lbm::interpolateCF(f, omegaF, epsilon_new, coefficients, x, y, z);

    setPreCollisionDistribution(distFine, indices, f);
}

//////////////////////////////////////////////////////////////////////////
//! \brief Interpolate from coarse to fine nodes
//! \details This scaling function is designed for the Cumulant K17 Kernel chimera collision kernel.
//!
//! The function is executed in the following steps:
//!
template<bool hasTurbulentViscosity> __global__ void scaleCoarseToFineCompressible_Device(
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
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if (nodeIndex >= numberOfInterfaceNodes)
        return;

    // 1.calculate moments
    vf::lbm::MomentsOnSourceNodeSet momentsSet;

    // - ###########################################
    // TODO: The following part starting from here
    //       should be exchanged by function call of calculate_moment_set<hasTurbulentViscosity>()
    //       However, this leads to a very small deviation in the results (~10e-6).
    real omegaC  = omegaCoarse;
    Distributions27 distCoarse;
    vf::gpu::getPointersToDistributions(distCoarse, distributionsCoarse, numberOfLBnodesCoarse, isEvenTimestep);

    vf::gpu::ListIndices indices;
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
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    real f_coarse[27];

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculateMMM(f_coarse, omegaC);

    //////////////////////////////////////////////////////////////////////////
    // source node TSW = MMP
    //////////////////////////////////////////////////////////////////////////
    // Set neighbor indices - has to be recalculated for the new source node
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZcoarse[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZcoarse[indices.k_M0M];
    indices.k_0MM = neighborZcoarse[indices.k_0MM];
    indices.k_MMM = neighborZcoarse[indices.k_MMM];

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculateMMP(f_coarse, omegaC);

    //////////////////////////////////////////////////////////////////////////
    // source node TSE = PMP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXcoarse[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXcoarse[indices.k_MM0];
    indices.k_M0M = neighborXcoarse[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXcoarse[indices.k_MMM];

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculatePMP(f_coarse, omegaC);

    //////////////////////////////////////////////////////////////////////////
    // source node BSE = PMM 
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXcoarse[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXcoarse[k_base_MM0];

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculatePMM(f_coarse, omegaC);

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
    indices.k_000 = k_base_000;
    indices.k_M00 = k_base_M00;
    indices.k_0M0 = k_base_0M0;
    indices.k_00M = k_base_00M;
    indices.k_MM0 = k_base_MM0;
    indices.k_M0M = k_base_M0M;
    indices.k_0MM = k_base_0MM;
    indices.k_MMM = k_base_MMM;

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculateMPM(f_coarse, omegaC);

    //////////////////////////////////////////////////////////////////////////
    // source node TNW = MPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_00M;
    indices.k_M00 = indices.k_M0M;
    indices.k_0M0 = indices.k_0MM;
    indices.k_00M = neighborZcoarse[indices.k_00M];
    indices.k_MM0 = indices.k_MMM;
    indices.k_M0M = neighborZcoarse[indices.k_M0M];
    indices.k_0MM = neighborZcoarse[indices.k_0MM];
    indices.k_MMM = neighborZcoarse[indices.k_MMM];

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculateMPP(f_coarse, omegaC);
    //////////////////////////////////////////////////////////////////////////
    // source node TNE = PPP
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_000 = indices.k_M00;
    indices.k_M00 = neighborXcoarse[indices.k_M00];
    indices.k_0M0 = indices.k_MM0;
    indices.k_00M = indices.k_M0M;
    indices.k_MM0 = neighborXcoarse[indices.k_MM0];
    indices.k_M0M = neighborXcoarse[indices.k_M0M];
    indices.k_0MM = indices.k_MMM;
    indices.k_MMM = neighborXcoarse[indices.k_MMM];

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculatePPP(f_coarse, omegaC);
    //////////////////////////////////////////////////////////////////////////
    // source node BNE = PPM
    //////////////////////////////////////////////////////////////////////////
    // index
    indices.k_00M = indices.k_000;
    indices.k_M0M = indices.k_M00;
    indices.k_0MM = indices.k_0M0;
    indices.k_MMM = indices.k_MM0;
    indices.k_000 = k_base_M00;
    indices.k_M00 = neighborXcoarse[k_base_M00];
    indices.k_0M0 = k_base_MM0;
    indices.k_MM0 = neighborXcoarse[k_base_MM0];

    omegaC = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;

    getPreCollisionDistribution(f_coarse, distCoarse, indices);
    momentsSet.calculatePPM(f_coarse, omegaC);
    // should be extractel until this line
    // - ###################################################################

    // 2.calculate coefficients
    vf::lbm::InterpolationCoefficients coefficients;
    momentsSet.calculateCoefficients(coefficients, neighborCoarseToFine.x[nodeIndex], neighborCoarseToFine.y[nodeIndex], neighborCoarseToFine.z[nodeIndex]);

    // 3. interpolate coarse to fine
    interpolate<hasTurbulentViscosity>(
        coefficients,
        nodeIndex,
        distributionsFine, 
        neighborXfine,
        neighborYfine,
        neighborZfine,
        numberOfLBnodesFine,
        indicesFineMMM,
        omegaFine,
        turbulentViscosityFine);
}

template __global__ void scaleCoarseToFineCompressible_Device<true>( real* distributionsCoarse, real* distributionsFine, unsigned int* neighborXcoarse, unsigned int* neighborYcoarse, unsigned int* neighborZcoarse, unsigned int* neighborXfine, unsigned int* neighborYfine, unsigned int* neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int* indicesCoarseMMM, unsigned int* indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh offsetCF);

template __global__ void scaleCoarseToFineCompressible_Device<false>( real* distributionsCoarse, real* distributionsFine, unsigned int* neighborXcoarse, unsigned int* neighborYcoarse, unsigned int* neighborZcoarse, unsigned int* neighborXfine, unsigned int* neighborYfine, unsigned int* neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int* indicesCoarseMMM, unsigned int* indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh offsetCF);
//! \}
