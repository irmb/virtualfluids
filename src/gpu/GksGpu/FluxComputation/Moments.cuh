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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file Moments.cuh
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#ifndef Moments_CUH
#define Moments_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#define NUMBER_OF_MOMENTS 7

//! \brief Computes the moments of the equilibrium distribution
//!
//! The equations for the computation of the moments can be found in
//! Appendix C in 
//! <a href="https://doi.org/10.1142/9324"><b>[ Kun Xu, (2015), DOI: 10.1142/9324 ]</b></a>.
//! 
//! The moments do not contain 
//! 
//! \param[in]  facePrim   flow state on the interface as \ref PrimitiveVariables
//! \param[in]  K          number of internal degrees of freedom
//! \param[out] momentU    array of moments of the equilibrium function with respect to x-velocity
//! \param[out] momentV    array of moments of the equilibrium function with respect to y-velocity
//! \param[out] momentW    array of moments of the equilibrium function with respect to z-velocity
//! \param[out] momentXi   array of moments of the equilibrium function with respect to internal degrees of freedom
__host__ __device__ inline void computeMoments( const PrimitiveVariables & facePrim,
                                                const real K,
                                                real momentU [NUMBER_OF_MOMENTS], 
                                                real momentV [NUMBER_OF_MOMENTS], 
                                                real momentW [NUMBER_OF_MOMENTS], 
                                                real momentXi[NUMBER_OF_MOMENTS] )
{
    momentU[0] = c1o1;
    momentU[1] = facePrim.U;
#pragma unroll
    for ( uint i = 2; i < NUMBER_OF_MOMENTS; i++ )
        momentU[i] = facePrim.U * momentU[i - 1] + ( real(i - 1) * momentU[i - 2] )/( c2o1 * facePrim.lambda );

    momentV[0] = c1o1;
    momentV[1] = facePrim.V;
#pragma unroll
    for ( uint i = 2; i < NUMBER_OF_MOMENTS; i++ )
        momentV[i] = facePrim.V * momentV[i - 1] + ( real(i - 1) * momentV[i - 2] )/( c2o1 * facePrim.lambda );

    momentW[0] = c1o1;
    momentW[1] = facePrim.W;
#pragma unroll
    for ( uint i = 2; i < NUMBER_OF_MOMENTS; i++ )
        momentW[i] = facePrim.W * momentW[i - 1] + ( real(i - 1) * momentW[i - 2] )/( c2o1 * facePrim.lambda );

    momentXi[0] = c1o1;
    momentXi[1] = c0o1;
    momentXi[2] = K / ( c2o1 * facePrim.lambda );
    momentXi[3] = c0o1;
    momentXi[4] = K * ( c2o1 + K ) / ( c4o1 * facePrim.lambda * facePrim.lambda );
    momentXi[5] = c0o1;
    momentXi[6] = ( K + c4o1 ) / ( c2o1 * facePrim.lambda ) * momentXi[4];
}



#endif