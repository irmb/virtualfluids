#ifndef FlowStateDataConversion_H
#define FlowStateDataConversion_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#endif

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/ThermalDependencies.cuh"

namespace GksGpu {

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline ConservedVariables toConservedVariables( const PrimitiveVariables& prim, real K, bool overrideK = true )
{
    //#ifdef USE_PASSIVE_SCALAR
    //if( overrideK ) K = getK(prim);
    //#endif

#ifdef USE_PASSIVE_SCALAR
    return ConservedVariables(prim.rho
                             ,prim.U * prim.rho
                             ,prim.V * prim.rho
                             ,prim.W * prim.rho
                             //,getEint(prim) * prim.rho + c1o2 * prim.rho * ( prim.U * prim.U + prim.V * prim.V + prim.W * prim.W )
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

__host__ __device__ inline PrimitiveVariables toPrimitiveVariables( const ConservedVariables& cons, real K, bool overrideK = true )
{
    //#ifdef USE_PASSIVE_SCALAR
    //if( overrideK ) K = getK(cons);
    //#endif

#ifdef USE_PASSIVE_SCALAR
	return PrimitiveVariables(cons.rho
						     ,cons.rhoU / cons.rho
						     ,cons.rhoV / cons.rho
						     ,cons.rhoW / cons.rho
						     //,getlambda(cons)
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

} // namespace GksGpu

#endif

