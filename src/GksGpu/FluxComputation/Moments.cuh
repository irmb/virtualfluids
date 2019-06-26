#ifndef Moments_CUH
#define Moments_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#define NUMBER_OF_MOMENTS 7

__host__ __device__ inline void computeMoments( const PrimitiveVariables & facePrim,
                                                const real K,
                                                real momentU [NUMBER_OF_MOMENTS], 
                                                real momentV [NUMBER_OF_MOMENTS], 
                                                real momentW [NUMBER_OF_MOMENTS], 
                                                real momentXi[NUMBER_OF_MOMENTS] )
{
    momentU[0] = one;
    momentU[1] = facePrim.U;
#pragma unroll
    for ( uint i = 2; i < NUMBER_OF_MOMENTS; i++ )
        momentU[i] = facePrim.U * momentU[i - 1] + ( real(i - 1) * momentU[i - 2] )/( two * facePrim.lambda );

    momentV[0] = one;
    momentV[1] = facePrim.V;
#pragma unroll
    for ( uint i = 2; i < NUMBER_OF_MOMENTS; i++ )
        momentV[i] = facePrim.V * momentV[i - 1] + ( real(i - 1) * momentV[i - 2] )/( two * facePrim.lambda );

    momentW[0] = one;
    momentW[1] = facePrim.W;
#pragma unroll
    for ( uint i = 2; i < NUMBER_OF_MOMENTS; i++ )
        momentW[i] = facePrim.W * momentW[i - 1] + ( real(i - 1) * momentW[i - 2] )/( two * facePrim.lambda );

    momentXi[0] = one;
    momentXi[1] = zero;
    momentXi[2] = K / ( two * facePrim.lambda );
    momentXi[3] = zero;
    momentXi[4] = K * ( two + K ) / ( four * facePrim.lambda * facePrim.lambda );
    momentXi[5] = zero;
    momentXi[6] = ( K + four ) / ( two * facePrim.lambda ) * momentXi[4];
}



#endif