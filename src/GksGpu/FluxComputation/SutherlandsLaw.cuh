#ifndef SutherlandsLaw_CUH
#define SutherlandsLaw_CUH

#include <cmath>

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

inline __host__ __device__ real sutherlandsLaw(const Parameters & parameters, const real r)
{
    real S  = real( 110.5 );

    real T0 = real( 600.0 );

    real C = S / T0;

    return parameters.mu * sqrt( r * r * r ) * ( C  + one ) / ( r  + C );
}

#endif
