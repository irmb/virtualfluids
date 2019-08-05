#include "CellUpdate.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/ThermalDependencies.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void chemicalReaction(DataBaseStruct dataBase, Parameters parameters, uint cellIndex, ConservedVariables& cons)
{
    // see FDS 5 Technical reference guide, section 6.1.4 for combustion model
#ifdef USE_PASSIVE_SCALAR
    if (parameters.enableReaction)
    {
        CellProperties cellProperties = dataBase.cellProperties[ cellIndex ];

        if( isCellProperties( cellProperties, CELL_PROPERTIES_FINE_GHOST ) ) return;

        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

        //////////////////////////////////////////////////////////////////////////

        real diffusivity = dataBase.diffusivity[ cellIndex ] / ( c6o1 * parameters.dx * parameters.dx * parameters.dt );
        dataBase.diffusivity[ cellIndex ] = c0o1;

        //////////////////////////////////////////////////////////////////////////

        real mixingTimeScale = real(0.1) * parameters.dx * parameters.dx / diffusivity;

        //real mixingTimeScale = parameters.dt;

        //if( mixingTimeScale < one )
        //    mixingTimeScale = one;

        //////////////////////////////////////////////////////////////////////////

        real Y_F = prim.S_1;
        real Y_P = prim.S_2;

        real Y_A = c1o1 - Y_F - Y_P;

        ///////////////////////////////////////////////////////////////////////////////

        real Y_O2 = rX * ( M_O2 / M_A ) * Y_A;

        ///////////////////////////////////////////////////////////////////////////////

        real s = M_F / ( c2o1 * M_O2 );

        real heatReleaseRate = cons.rho * fminf(Y_F, s * Y_O2) / mixingTimeScale * ( parameters.heatOfReaction / M_F );

        //////////////////////////////////////////////////////////////////////////

        if( heatReleaseRate < c0o1 )
            heatReleaseRate = c0o1;

        //////////////////////////////////////////////////////////////////////////

        if( parameters.useHeatReleaseRateLimiter )
        if( heatReleaseRate > parameters.heatReleaseRateLimiter )
            heatReleaseRate = parameters.heatReleaseRateLimiter;

        //////////////////////////////////////////////////////////////////////////

        real drhoY_F = heatReleaseRate * parameters.dt / ( parameters.heatOfReaction / M_F );

        //real r = c1o1 + ( c1o2 / rX ) * ( M_A / M_F );
        real r = c1o1 + ( c2o1 / rX ) * ( M_A / M_F );

        cons.rhoS_1 -=     drhoY_F;
        cons.rhoS_2 += r * drhoY_F;
        cons.rhoE   += heatReleaseRate * parameters.dt;

        //////////////////////////////////////////////////////////////////////////
    }

#endif // USE_PASSIVE_SCALAR
}
