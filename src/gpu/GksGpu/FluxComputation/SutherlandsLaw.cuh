#ifndef SutherlandsLaw_CUH
#define SutherlandsLaw_CUH

#include <cmath>


#include "GksGpu_export.h"

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

namespace GksGpu {

inline __host__ __device__ real sutherlandsLaw(Parameters & parameters, const real r)
{
    real S  = real( 110.5 );

    real T0 = real( 600.0 );

    real C = S / T0;

    return parameters.mu * sqrt( r * r * r ) * ( C  + c1o1 ) / ( r  + C );
}

inline __host__ __device__ real sutherlandsLaw2(Parameters & parameters, const real r)
{
    real Smu = real( 0.648 );

    real Sk  = real( 0.368 );

    parameters.Pr *= ( ( Smu  + c1o1 ) / ( Sk  + c1o1 ) ) * ( ( r  + Sk ) / ( r  + Smu ) );

    return parameters.mu * sqrt( r * r * r ) * ( Smu  + c1o1 ) / ( r  + Smu );
}

inline __host__ __device__ real getViscosity(Parameters & parameters, const real r)
{
    if ( parameters.viscosityModel == ViscosityModel::sutherlandsLaw ){
        return sutherlandsLaw( parameters, r );
    }
    else if ( parameters.viscosityModel == ViscosityModel::sutherlandsLaw2 ){
        return sutherlandsLaw2( parameters, r );
    }

    return parameters.mu;
}

} // namespace GksGpu

#endif
