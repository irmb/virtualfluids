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
//! \file ExpansionCoefficients.cuh
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#ifndef ExpansionCoefficients_CUH
#define ExpansionCoefficients_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

//! \brief computes spatial and temporal expansion coefficients of the Taylor 
//! expansion of the equilibrium distribution based on corresponding derivatives
//!
//! The equations for the computation of the expansion coefficints can be found in
//! Appendix C in 
//! <a href="https://doi.org/10.1142/9324"><b>[ Kun Xu, (2015), DOI: 10.1142/9324 ]</b></a>.
//!
//! \param[in]  facePrim               flow state on the interface as \ref PrimitiveVariables
//! \param[in]  gradient               derivative of the \ref ConservedVariables
//! \param[in]  K                      number of internal degrees of freedom
//! \param[out] expansionCoefficient   array of expansion coefficients with respect to the provided derivatives
__host__ __device__ inline void computeExpansionCoefficients(const PrimitiveVariables & facePrim, 
                                                             const ConservedVariables & gradient,
                                                             const real K, 
                                                             real expansionCoefficient[LENGTH_CELL_DATA])
{
    real two_E, 
             rho_dU_dx, 
             rho_dV_dx, 
             rho_dW_dx, 
         two_rho_dE_dx;

    two_E = facePrim.U * facePrim.U 
          + facePrim.V * facePrim.V 
          + facePrim.W * facePrim.W 
          + c1o2 * ( K + c3o1 ) / facePrim.lambda;

    rho_dU_dx     =       gradient.rhoU - facePrim.U  * gradient.rho;
    rho_dV_dx     =       gradient.rhoV - facePrim.V  * gradient.rho;
    rho_dW_dx     =       gradient.rhoW - facePrim.W  * gradient.rho;
    two_rho_dE_dx = c2o1 * gradient.rhoE -      two_E  * gradient.rho;

    expansionCoefficient[4] = ( c4o1 * facePrim.lambda * facePrim.lambda ) / ( K + c3o1 )
                            * ( two_rho_dE_dx - c2o1 * facePrim.U * rho_dU_dx 
                                              - c2o1 * facePrim.V * rho_dV_dx 
                                              - c2o1 * facePrim.W * rho_dW_dx );

    expansionCoefficient[3] = c2o1 * facePrim.lambda * rho_dW_dx - facePrim.W * expansionCoefficient[4];

    expansionCoefficient[2] = c2o1 * facePrim.lambda * rho_dV_dx - facePrim.V * expansionCoefficient[4];

    expansionCoefficient[1] = c2o1 * facePrim.lambda * rho_dU_dx - facePrim.U * expansionCoefficient[4];

    expansionCoefficient[0] = gradient.rho -   facePrim.U * expansionCoefficient[1] 
                                           -   facePrim.V * expansionCoefficient[2] 
                                           -   facePrim.W * expansionCoefficient[3] 
                                           - c1o2 * two_E * expansionCoefficient[4];

#ifdef USE_PASSIVE_SCALAR
	expansionCoefficient[5] = c2o1 * facePrim.lambda * (gradient.rhoS_1 - facePrim.S_1 * gradient.rho);
	expansionCoefficient[6] = c2o1 * facePrim.lambda * (gradient.rhoS_2 - facePrim.S_2 * gradient.rho);
#endif // USE_PASSIVE_SCALAR
}

#endif