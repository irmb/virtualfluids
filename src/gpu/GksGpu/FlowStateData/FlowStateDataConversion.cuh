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
//! \file FlowStateDataConversion.cuh
//! \ingroup FlowStateData
//! \author Stephan Lenz
//=======================================================================================
#ifndef FlowStateDataConversion_H
#define FlowStateDataConversion_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//! \brief transforms PrimitiveVariables objects to ConservedVariables objects
//! \param prim        PrimitiveVariables object
//! \param K           number of internal degrees of freedom
//! \param overrideK   deprecated argument
//! \return the corresponding ConservedVariables of the PrimitiveVariables prim
__host__ __device__ inline ConservedVariables toConservedVariables( const PrimitiveVariables& prim, real K, bool overrideK = true )
{
#ifdef USE_PASSIVE_SCALAR
    return ConservedVariables(prim.rho
                             ,prim.U * prim.rho
                             ,prim.V * prim.rho
                             ,prim.W * prim.rho
                             ,( K + c3o1 ) / ( c4o1 * prim.lambda ) * prim.rho + c1o2 * prim.rho * ( prim.U * prim.U + prim.V * prim.V + prim.W * prim.W )
                             ,prim.S_1 * prim.rho
                             ,prim.S_2 * prim.rho
    );
#else
    return ConservedVariables(prim.rho
                             ,prim.U * prim.rho
                             ,prim.V * prim.rho
                             ,prim.W * prim.rho
                             ,( K + c3o1 ) / ( c4o1 * prim.lambda ) * prim.rho + c1o2 * prim.rho * ( prim.U * prim.U + prim.V * prim.V + prim.W * prim.W )
    );
#endif
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//! \brief transforms ConservedVariables objects to PrimitiveVariables objects
//! \param cons        ConservedVariables object
//! \param K           number of internal degrees of freedom
//! \param overrideK   deprecated argument
//! \return the corresponding PrimitiveVariables of the ConservedVariables prim
__host__ __device__ inline PrimitiveVariables toPrimitiveVariables( const ConservedVariables& cons, real K, bool overrideK = true )
{

#ifdef USE_PASSIVE_SCALAR
	return PrimitiveVariables(cons.rho
						     ,cons.rhoU / cons.rho
						     ,cons.rhoV / cons.rho
						     ,cons.rhoW / cons.rho
						     ,( K + c3o1 ) * cons.rho / ( c4o1 * ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) )
                             ,cons.rhoS_1 / cons.rho
                             ,cons.rhoS_2 / cons.rho
	);
#else
	return PrimitiveVariables(cons.rho
						     ,cons.rhoU / cons.rho
						     ,cons.rhoV / cons.rho
						     ,cons.rhoW / cons.rho
						     ,( K + c3o1 ) * cons.rho / ( c4o1 * ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) )
	);
#endif
}

#endif

