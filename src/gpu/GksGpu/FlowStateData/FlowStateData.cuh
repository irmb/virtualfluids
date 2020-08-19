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
//! \file FlowStateData.cuh
//! \ingroup FlowStateData
//! \author Stephan Lenz
//=======================================================================================
#ifndef FlowStateData_H
#define FlowStateData_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Definitions/PassiveScalar.h"

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//! \brief Holds flow state data in terms of primitive variables
struct PrimitiveVariables
{
    real rho;
    real U;
    real V;
    real W;
    real lambda;
    #ifdef USE_PASSIVE_SCALAR
    real S_1;
    real S_2;
    #endif

    //////////////////////////////////////////////////////////////////////////

    //! default constructor that sets all values to zero
    __host__ __device__ PrimitiveVariables()
		: rho   (c0o1)
         ,U     (c0o1)
         ,V     (c0o1)
         ,W     (c0o1)
         ,lambda(c0o1)
    #ifdef USE_PASSIVE_SCALAR
         ,S_1   (c0o1)
         ,S_2   (c0o1)
    #endif
    {}

    //////////////////////////////////////////////////////////////////////////

    //! constructor that initializes the variables according to the arguments
    __host__ __device__ PrimitiveVariables(real rho
                                          ,real U
                                          ,real V
                                          ,real W
                                          ,real lambda
    #ifdef USE_PASSIVE_SCALAR
                                          ,real S_1 = c0o1
                                          ,real S_2 = c0o1
    #endif
    )
        : rho   (rho   )
         ,U     (U     )
         ,V     (V     )
         ,W     (W     )
         ,lambda(lambda)
    #ifdef USE_PASSIVE_SCALAR
         ,S_1   (S_1   )
         ,S_2   (S_2   )
    #endif
    {}

    //////////////////////////////////////////////////////////////////////////
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//! \brief Holds flow state data in terms of conserved variables
struct ConservedVariables
{
    real rho;
    real rhoU;
    real rhoV;
    real rhoW;
    real rhoE;
    #ifdef USE_PASSIVE_SCALAR
    real rhoS_1;
    real rhoS_2;
    #endif

    //////////////////////////////////////////////////////////////////////////

    //! default constructor that sets all values to zero
    __host__ __device__ ConservedVariables()
        : rho (c0o1)
         ,rhoU(c0o1)
         ,rhoV(c0o1)
         ,rhoW(c0o1)
         ,rhoE(c0o1)
    #ifdef USE_PASSIVE_SCALAR
         ,rhoS_1(c0o1)
         ,rhoS_2(c0o1)
    #endif
    {}

    //////////////////////////////////////////////////////////////////////////
		  
    //! constructor that initializes the variables according to the arguments
    __host__ __device__ ConservedVariables(real rho
                                          ,real rhoU
                                          ,real rhoV
                                          ,real rhoW
                                          ,real rhoE
    #ifdef USE_PASSIVE_SCALAR
                                          ,real rhoS_1 = c0o1
                                          ,real rhoS_2 = c0o1
    #endif
    )
        : rho (rho )
         ,rhoU(rhoU)
         ,rhoV(rhoV)
         ,rhoW(rhoW)
         ,rhoE(rhoE)
    #ifdef USE_PASSIVE_SCALAR
         ,rhoS_1(rhoS_1)
         ,rhoS_2(rhoS_2)
    #endif
    {}

    //////////////////////////////////////////////////////////////////////////
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \return the sum of all components of the two PrimitiveVariables left and right as PrimitiveVariables object
__host__ __device__ inline PrimitiveVariables operator+ ( const PrimitiveVariables& left, const PrimitiveVariables& right )
{
    PrimitiveVariables result;

    result.rho    = left.rho    + right.rho   ;
    result.U      = left.U      + right.U     ;
    result.V      = left.V      + right.V     ;
    result.W      = left.W      + right.W     ;
    result.lambda = left.lambda + right.lambda;

#ifdef USE_PASSIVE_SCALAR
    result.S_1    = left.S_1    + right.S_1   ;
    result.S_2    = left.S_2    + right.S_2   ;
#endif

    return result;
}

//! \return the sum of all components of two ConservedVariables left and right as ConservedVariables object
__host__ __device__ inline ConservedVariables operator+ ( const ConservedVariables& left, const ConservedVariables& right )
{
    ConservedVariables result;

    result.rho    = left.rho    + right.rho   ;
    result.rhoU   = left.rhoU   + right.rhoU  ;
    result.rhoV   = left.rhoV   + right.rhoV  ;
    result.rhoW   = left.rhoW   + right.rhoW  ;
    result.rhoE   = left.rhoE   + right.rhoE  ;

#ifdef USE_PASSIVE_SCALAR
    result.rhoS_1 = left.rhoS_1 + right.rhoS_1;
    result.rhoS_2 = left.rhoS_2 + right.rhoS_2;
#endif

    return result;
}

//////////////////////////////////////////////////////////////////////////

//! \return the difference of all components of the two PrimitiveVariables left and right as PrimitiveVariables object
__host__ __device__ inline PrimitiveVariables operator- ( const PrimitiveVariables& left, const PrimitiveVariables& right )
{
    PrimitiveVariables result;

    result.rho    = left.rho    - right.rho   ;
    result.U      = left.U      - right.U     ;
    result.V      = left.V      - right.V     ;
    result.W      = left.W      - right.W     ;
    result.lambda = left.lambda - right.lambda;

#ifdef USE_PASSIVE_SCALAR
    result.S_1    = left.S_1    - right.S_1   ;
    result.S_2    = left.S_2    - right.S_2   ;
#endif

    return result;
}

//! \return the difference of all components of two ConservedVariables left and right as ConservedVariables object
__host__ __device__ inline ConservedVariables operator- ( const ConservedVariables& left, const ConservedVariables& right )
{
    ConservedVariables result;

    result.rho    = left.rho    - right.rho   ;
    result.rhoU   = left.rhoU   - right.rhoU  ;
    result.rhoV   = left.rhoV   - right.rhoV  ;
    result.rhoW   = left.rhoW   - right.rhoW  ;
    result.rhoE   = left.rhoE   - right.rhoE  ;

#ifdef USE_PASSIVE_SCALAR
    result.rhoS_1 = left.rhoS_1 - right.rhoS_1;
    result.rhoS_2 = left.rhoS_2 - right.rhoS_2;
#endif

    return result;
}

//////////////////////////////////////////////////////////////////////////

//! \return the product of the scalar left and all components of the PrimitiveVariables right as PrimitiveVariables object
__host__ __device__ inline PrimitiveVariables operator* ( const real left, const PrimitiveVariables& right )
{
    PrimitiveVariables result;

    result.rho    = left * right.rho   ;
    result.U      = left * right.U     ;
    result.V      = left * right.V     ;
    result.W      = left * right.W     ;
    result.lambda = left * right.lambda;

#ifdef USE_PASSIVE_SCALAR
    result.S_1    = left * right.S_1   ;
    result.S_2    = left * right.S_2   ;
#endif

    return result;
}

//! \return the product of the scalar left and all components of the ConservedVariables right as ConservedVariables object
__host__ __device__ inline ConservedVariables operator* ( const real left, const ConservedVariables& right )
{
    ConservedVariables result;

    result.rho    = left * right.rho   ;
    result.rhoU   = left * right.rhoU  ;
    result.rhoV   = left * right.rhoV  ;
    result.rhoW   = left * right.rhoW  ;
    result.rhoE   = left * right.rhoE  ;

#ifdef USE_PASSIVE_SCALAR
    result.rhoS_1 = left * right.rhoS_1;
    result.rhoS_2 = left * right.rhoS_2;
#endif

    return result;
}

#endif

