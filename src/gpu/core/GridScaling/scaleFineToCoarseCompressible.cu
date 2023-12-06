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
//! \author Martin Schoenherr, Anna Wellmann
//=======================================================================================

#include "Utilities/KernelUtilities.h"
#include "Utilities/ScalingUtilities.h"

#include <lbm/refinement/InterpolationFC.h>
#include <lbm/interpolation/InterpolationCoefficients.h>


template <bool hasTurbulentViscosity> __device__ void interpolate(
    vf::lbm::InterpolationCoefficients& coefficients,
    const unsigned int nodeIndex,
    real* distributionsCoarse, 
    unsigned int* neighborXcoarse,
    unsigned int* neighborYcoarse,
    unsigned int* neighborZcoarse,
    unsigned long long numberOfLBnodesCoarse,
    unsigned int* indicesCoarse000,
    real omegaCoarse,
    real* turbulentViscosityCoarse,
    bool isEvenTimestep
)
{
    // Position Coarse 0., 0., 0.
    Distributions27 distCoarse;
    vf::gpu::getPointersToDistributions(distCoarse, distributionsCoarse, numberOfLBnodesCoarse, isEvenTimestep);

    vf::gpu::ListIndices indices(indicesCoarse000[nodeIndex], neighborXcoarse, neighborYcoarse, neighborZcoarse);

    const real epsilonNew = c2o1; // ratio of grid resolutions
    const real omegaCoarseNew = hasTurbulentViscosity ? vf::lbm::calculateOmegaWithturbulentViscosity(omegaCoarse, turbulentViscosityCoarse[indices.k_000]) : omegaCoarse;
    real fCoarse[27];
    vf::lbm::interpolateFC(fCoarse, epsilonNew, omegaCoarseNew, coefficients);

    vf::gpu::setPreCollisionDistribution(distCoarse, indices, fCoarse);
}



//////////////////////////////////////////////////////////////////////////
//! \brief Interpolate from fine to coarse
//! \details This scaling function is designed for the Cumulant K17 Kernel chimera collision kernel
template<bool hasTurbulentViscosity> __global__ void scaleFineToCoarseCompressible_Device(
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
    const real omega_coarse,
    const real omegaFine,
    real* turbulentViscosityCoarse,
    real* turbulentViscosityFine,
    ICellNeigh neighborFineToCoarse)
{
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if (nodeIndex >= numberOfInterfaceNodes)
        return;

    // 1.calculate moments
    vf::lbm::MomentsOnSourceNodeSet momentsSet;
    vf::gpu::calculateMomentSet<hasTurbulentViscosity>(
        momentsSet, nodeIndex, distributionsFine, neighborXfine, neighborYfine, neighborZfine, indicesFineMMM, turbulentViscosityFine, numberOfLBnodesFine, omegaFine, true);

    // 2.calculate coefficients
    vf::lbm::InterpolationCoefficients coefficients;
    momentsSet.calculateCoefficients(coefficients, neighborFineToCoarse.x[nodeIndex], neighborFineToCoarse.y[nodeIndex], neighborFineToCoarse.z[nodeIndex]);

    // 3. interpolate fine to coarse
    interpolate<hasTurbulentViscosity>(
        coefficients,
        nodeIndex,
        distributionsCoarse, 
        neighborXcoarse,
        neighborYcoarse,
        neighborZcoarse,
        numberOfLBnodesCoarse,
        indicesCoarse000,
        omega_coarse,
        turbulentViscosityCoarse,
        isEvenTimestep);
}

template __global__ void scaleFineToCoarseCompressible_Device<true>(real *distributionsCoarse, real *distributionsFine, unsigned int *neighborXcoarse, unsigned int *neighborYcoarse, unsigned int *neighborZcoarse, unsigned int *neighborXfine, unsigned int *neighborYfine, unsigned int *neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int *indicesCoarse000, unsigned int *indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh neighborFineToCoarse);

template __global__ void scaleFineToCoarseCompressible_Device<false>(real *distributionsCoarse, real *distributionsFine, unsigned int *neighborXcoarse, unsigned int *neighborYcoarse, unsigned int *neighborZcoarse, unsigned int *neighborXfine, unsigned int *neighborYfine, unsigned int *neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int *indicesCoarse000, unsigned int *indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh neighborFineToCoarse);