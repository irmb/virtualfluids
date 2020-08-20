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
//! \file SutherlandsLaw.cuh
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#ifndef SutherlandsLaw_CUH
#define SutherlandsLaw_CUH

#include <cmath>

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

//! \brief Computes temperature dependent viscosity based on Sutherlands law
//! 
//! \param parameters       \ref Parameters struct
//! \param r                ratio of temperature and reference temperature at which viscosity is defined
//! 
//! \return viscosity for given temperature ratio
inline __host__ __device__ real sutherlandsLaw(const Parameters & parameters, const real r)
{
    real S  = real( 110.5 );

    real T0 = real( 600.0 );

    real C = S / T0;

    return parameters.mu * sqrt( r * r * r ) * ( C  + c1o1 ) / ( r  + C );
}

#endif
