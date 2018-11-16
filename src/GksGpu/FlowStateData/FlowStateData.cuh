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
    real S;
    #endif

    //////////////////////////////////////////////////////////////////////////

    __host__ __device__ PrimitiveVariables()
		: rho   (zero)
         ,U     (zero)
         ,V     (zero)
         ,W     (zero)
         ,lambda(zero)
    #ifdef USE_PASSIVE_SCALAR
         ,S     (zero)
    #endif
    {}

    //////////////////////////////////////////////////////////////////////////

    __host__ __device__ PrimitiveVariables(real rho
                                          ,real U
                                          ,real V
                                          ,real W
                                          ,real lambda
    #ifdef USE_PASSIVE_SCALAR
                                          ,real S
    #endif
    )
        : rho   (rho   )
         ,U     (U     )
         ,V     (V     )
         ,W     (W     )
         ,lambda(lambda)
    #ifdef USE_PASSIVE_SCALAR
         ,S     (S     )
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
    real rhoS;
    #endif

    //////////////////////////////////////////////////////////////////////////

    __host__ __device__ ConservedVariables()
        : rho (zero)
         ,rhoU(zero)
         ,rhoV(zero)
         ,rhoW(zero)
         ,rhoE(zero)
    #ifdef USE_PASSIVE_SCALAR
         ,rhoS(zero)
    #endif
    {}

    //////////////////////////////////////////////////////////////////////////
		  
    __host__ __device__ ConservedVariables(real rho
                                          ,real rhoU
                                          ,real rhoV
                                          ,real rhoW
                                          ,real rhoE
    #ifdef USE_PASSIVE_SCALAR
                                          ,real rhoS
    #endif
    )
        : rho (rho )
         ,rhoU(rhoU)
         ,rhoV(rhoV)
         ,rhoW(rhoW)
         ,rhoE(rhoE)
    #ifdef USE_PASSIVE_SCALAR
         ,rhoS(rhoS)
    #endif
    {}

    //////////////////////////////////////////////////////////////////////////
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline ConservedVariables toConservedVariables( const PrimitiveVariables& prim, real K )
{
    return ConservedVariables(prim.rho
                             ,prim.U * prim.rho
                             ,prim.V * prim.rho
                             ,prim.W * prim.rho
                             ,( K + two ) / ( four * prim.lambda ) * prim.rho + c1o2 * prim.rho * ( prim.U * prim.U + prim.V * prim.V + prim.W * prim.W )
    #ifdef USE_PASSIVE_SCALAR
                             ,prim.S * prim.rho
    #endif
    );
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline PrimitiveVariables toPrimitiveVariables( const ConservedVariables& cons, real K )
{
	return PrimitiveVariables(cons.rho
						     ,cons.rhoU / cons.rho
						     ,cons.rhoV / cons.rho
						     ,cons.rhoW / cons.rho
						     ,( K + two ) * cons.rho / ( four * ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) )
    #ifdef USE_PASSIVE_SCALAR
                             ,cons.rhoS / cons.rho
    #endif
	);
}

#endif

