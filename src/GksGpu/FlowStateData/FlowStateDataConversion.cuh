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
#include "FlowStateData/ThermalDependencies.cuh"

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline ConservedVariables toConservedVariables( const PrimitiveVariables& prim, real K, bool overrideK = true )
{
    #ifdef USE_PASSIVE_SCALAR
    //if( overrideK ) K = getK(prim);
    #endif

    return ConservedVariables(prim.rho
                             ,prim.U * prim.rho
                             ,prim.V * prim.rho
                             ,prim.W * prim.rho
                             ,( K + three ) / ( four * prim.lambda ) * prim.rho + c1o2 * prim.rho * ( prim.U * prim.U + prim.V * prim.V + prim.W * prim.W )
    #ifdef USE_PASSIVE_SCALAR
                             ,prim.S_1 * prim.rho
                             ,prim.S_2 * prim.rho
    #endif
    );
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline PrimitiveVariables toPrimitiveVariables( const ConservedVariables& cons, real K, bool overrideK = true )
{
    #ifdef USE_PASSIVE_SCALAR
    //if( overrideK ) K = getK(cons);
    #endif

	return PrimitiveVariables(cons.rho
						     ,cons.rhoU / cons.rho
						     ,cons.rhoV / cons.rho
						     ,cons.rhoW / cons.rho
						     ,( K + three ) * cons.rho / ( four * ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) )
    #ifdef USE_PASSIVE_SCALAR
                             ,cons.rhoS_1 / cons.rho
                             ,cons.rhoS_2 / cons.rho
    #endif
	);
}

#endif

