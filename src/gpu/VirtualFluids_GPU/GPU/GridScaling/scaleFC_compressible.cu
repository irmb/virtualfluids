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

#include <lbm/KernelParameter.h>
#include <lbm/refinement/Interpolation_FC.h>
#include <lbm/refinement/Coefficients.h>


template <bool hasTurbulentViscosity> __device__ void interpolate(
    vf::lbm::Coefficients& coefficients,
    const unsigned int nodeIndex,
    real* distributionsCoarse, 
    unsigned int* neighborXcoarse,
    unsigned int* neighborYcoarse,
    unsigned int* neighborZcoarse,
    unsigned long long numberOfLBnodesCoarse,
    unsigned int* indicesCoarse000,
    real omega_coarse,
    real* turbulentViscosityCoarse,
    bool isEvenTimestep
)
{
    // Position Coarse 0., 0., 0.
    Distributions27 distCoarse;
    vf::gpu::getPointersToDistributions(distCoarse, distributionsCoarse, numberOfLBnodesCoarse, isEvenTimestep);

    vf::gpu::ListIndices indices;
    indices.k_000 = indicesCoarse000[nodeIndex];
    indices.k_M00 = neighborXcoarse [indices.k_000];
    indices.k_0M0 = neighborYcoarse [indices.k_000];
    indices.k_00M = neighborZcoarse [indices.k_000];
    indices.k_MM0 = neighborYcoarse [indices.k_M00];
    indices.k_M0M = neighborZcoarse [indices.k_M00];
    indices.k_0MM = neighborZcoarse [indices.k_0M0];
    indices.k_MMM = neighborZcoarse [indices.k_MM0];

    const real epsilon_new = c2o1; // ratio of grid resolutions
    const real omega_coarse_new = hasTurbulentViscosity ? vf::gpu::calculateOmega(omega_coarse, turbulentViscosityCoarse[indices.k_000]) : omega_coarse;
    real f_coarse[27];
    vf::lbm::interpolate_fc(f_coarse, epsilon_new, omega_coarse_new, coefficients);

    vf::gpu::write(distCoarse, indices, f_coarse);
}



//////////////////////////////////////////////////////////////////////////
//! \brief Interpolate from fine to coarse
//! \details This scaling function is designed for the Cumulant K17 Kernel chimera collision kernel
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
    vf::lbm::MomentsOnSourceNodeSet moments_set;
    vf::gpu::calculate_moment_set<hasTurbulentViscosity>(
        moments_set, nodeIndex, distributionsFine, neighborXfine, neighborYfine, neighborZfine, indicesFineMMM, turbulentViscosityFine, numberOfLBnodesFine, omegaFine, true);

    // 2.calculate coefficients
    vf::lbm::Coefficients coefficients;
    moments_set.calculateCoefficients(coefficients, neighborFineToCoarse.x[nodeIndex], neighborFineToCoarse.y[nodeIndex], neighborFineToCoarse.z[nodeIndex]);

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

template __global__ void scaleFC_compressible<true>(real *distributionsCoarse, real *distributionsFine, unsigned int *neighborXcoarse, unsigned int *neighborYcoarse, unsigned int *neighborZcoarse, unsigned int *neighborXfine, unsigned int *neighborYfine, unsigned int *neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int *indicesCoarse000, unsigned int *indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh neighborFineToCoarse);

template __global__ void scaleFC_compressible<false>(real *distributionsCoarse, real *distributionsFine, unsigned int *neighborXcoarse, unsigned int *neighborYcoarse, unsigned int *neighborZcoarse, unsigned int *neighborXfine, unsigned int *neighborYfine, unsigned int *neighborZfine, unsigned long long numberOfLBnodesCoarse, unsigned long long numberOfLBnodesFine, bool isEvenTimestep, unsigned int *indicesCoarse000, unsigned int *indicesFineMMM, unsigned int numberOfInterfaceNodes, real omegaCoarse, real omegaFine, real* turbulentViscosityCoarse, real* turbulentViscosityFine, ICellNeigh neighborFineToCoarse);