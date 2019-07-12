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

