#ifndef Moments_CUH
#define Moments_CUH


#include "GksGpu_export.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#define NUMBER_OF_MOMENTS 7

namespace GksGpu {

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

} // namespace GksGpu



#endif