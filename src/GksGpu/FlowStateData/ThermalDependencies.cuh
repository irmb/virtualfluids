#ifndef ThermalDependencies_H
#define ThermalDependencies_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

#include <Math.h>

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/HeatCapacities.cuh"


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

#ifdef USE_PASSIVE_SCALAR

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void getMassFractions( real Z1, real Z2, 
                                                  real& Y_O2, 
                                                  real& Y_N2, 
                                                  real& Y_CH4, 
                                                  real& Y_H2O, 
                                                  real& Y_CO2,
                                                  real& M_O2, 
                                                  real& M_N2, 
                                                  real& M_CH4, 
                                                  real& M_H2O, 
                                                  real& M_CO2 )
{
    M_O2  = 32.00e-3;
    M_N2  = 28.00e-3;
    M_CH4 = 16.00e-3;
    M_H2O = 18.00e-3;
    M_CO2 = 44.00e-3;

    const real Y_CH4_Inflow = 1.0;
    const real Y_N2_ambient = 0.767;
    const real Y_O2_ambient = 0.233;

    //////////////////////////////////////////////////////////////////////////

    const real Z = Z1 + Z2;
    
    Y_CH4 = Y_CH4_Inflow * Z1;
    Y_N2  = (1.0 - Z) * Y_N2_ambient;
    Y_O2  = (1.0 - Z) * Y_O2_ambient - 2.0  * ( M_O2  / M_CH4 ) * Y_CH4_Inflow * Z2;
    Y_H2O =                            2.0  * ( M_H2O / M_CH4 ) * Y_CH4_Inflow * Z2;
    Y_CO2 =                                   ( M_CO2 / M_CH4 ) * Y_CH4_Inflow * Z2;
}

__host__ __device__ inline void getMassFractions( const ConservedVariables& cons, 
                                                  real& Y_O2, 
                                                  real& Y_N2, 
                                                  real& Y_CH4, 
                                                  real& Y_H2O, 
                                                  real& Y_CO2,
                                                  real& M_O2, 
                                                  real& M_N2, 
                                                  real& M_CH4, 
                                                  real& M_H2O, 
                                                  real& M_CO2 )
{
    getMassFractions(cons.rhoS_1 / cons.rho,
                     cons.rhoS_2 / cons.rho,
                     Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                     M_O2, M_N2, M_CH4, M_H2O, M_CO2);
}

__host__ __device__ inline void getMassFractions( const PrimitiveVariables& prim, 
                                                  real& Y_O2, 
                                                  real& Y_N2, 
                                                  real& Y_CH4, 
                                                  real& Y_H2O, 
                                                  real& Y_CO2,
                                                  real& M_O2, 
                                                  real& M_N2, 
                                                  real& M_CH4, 
                                                  real& M_H2O, 
                                                  real& M_CO2 )
{
    getMassFractions(prim.S_1,
                     prim.S_2,
                     Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                     M_O2, M_N2, M_CH4, M_H2O, M_CO2);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getCp( const real T,
                                       const real Y_O2, 
                                       const real Y_N2, 
                                       const real Y_CH4, 
                                       const real Y_H2O, 
                                       const real Y_CO2 )
{
    return  Y_O2  * getCpO2 (T)
          + Y_N2  * getCpN2 (T)
          + Y_CH4 * getCpCH4(T)
          + Y_H2O * getCpH2O(T)
          + Y_CO2 * getCpCO2(T);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getK( const real T,
                                      const real Y_O2, 
                                      const real Y_N2, 
                                      const real Y_CH4, 
                                      const real Y_H2O, 
                                      const real Y_CO2 )
{
    real Ru = 8.31445984848;
    return 2.0 * getCp( T, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2 ) / Ru - 5.0;
}

__host__ __device__ inline real getK( const ConservedVariables& cons )
{
    real Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2;
    real M_O2, M_N2, M_CH4, M_H2O, M_CO2;

    getMassFractions(cons, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                           M_O2, M_N2, M_CH4, M_H2O, M_CO2);
    
    real M = Y_O2  * M_O2
           + Y_N2  * M_N2
           + Y_CH4 * M_CH4
           + Y_H2O * M_H2O
           + Y_CO2 * M_CO2;

    real Eint = ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) / cons.rho;

    real TAitken[3];

    real T0 = 600.0;
    real T  = T0;

    real Ru = 8.31445984848;

    //////////////////////////////////////////////////////////////////////////

    uint nIter = 5;
    uint iAitkenStart = nIter - 3;

#pragma unroll
    for( uint i = 0; i < nIter; i++ )
    {
        real Cp = getCp( T, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2 );
    
        T = 0.5 * ( T + ( Eint * M ) / ( Cp - Ru) );
        //T = ( Eint * M ) / ( Cp - Ru);

        if( i >= iAitkenStart ) TAitken[i - iAitkenStart] = T;
    }

    if( fabs( T - T0 ) / T0 > real(0.01) )
    {
        T = ( TAitken[1] * TAitken[1] - TAitken[0] * TAitken[2] )
          / ( TAitken[1] + TAitken[1] - TAitken[0] - TAitken[2] );
    }

    //////////////////////////////////////////////////////////////////////////

    return getK(T, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2);
}

__host__ __device__ inline real getK( const PrimitiveVariables& prim )
{
    real Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2;
    real M_O2, M_N2, M_CH4, M_H2O, M_CO2;

    getMassFractions(prim, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                           M_O2, M_N2, M_CH4, M_H2O, M_CO2);
    
    real M = Y_O2  * M_O2
           + Y_N2  * M_N2
           + Y_CH4 * M_CH4
           + Y_H2O * M_H2O
           + Y_CO2 * M_CO2;

    real Ru = 8.31445984848;

    //////////////////////////////////////////////////////////////////////////

    real T = M / ( two * prim.lambda * Ru );

    //////////////////////////////////////////////////////////////////////////

    return getK(T, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getT( const PrimitiveVariables& prim )
{
    real Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2;
    real M_O2, M_N2, M_CH4, M_H2O, M_CO2;

    getMassFractions(prim, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                           M_O2, M_N2, M_CH4, M_H2O, M_CO2);
    
    real M = Y_O2  * M_O2
           + Y_N2  * M_N2
           + Y_CH4 * M_CH4
           + Y_H2O * M_H2O
           + Y_CO2 * M_CO2;

    real Ru = 8.31445984848;

    return M / ( two * prim.lambda * Ru );
}

__host__ __device__ inline void setLambdaFromT( PrimitiveVariables& prim, real T )
{
    real Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2;
    real M_O2, M_N2, M_CH4, M_H2O, M_CO2;

    getMassFractions(prim, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                           M_O2, M_N2, M_CH4, M_H2O, M_CO2);
    
    real M = Y_O2  * M_O2
           + Y_N2  * M_N2
           + Y_CH4 * M_CH4
           + Y_H2O * M_H2O
           + Y_CO2 * M_CO2;

    real Ru = 8.31445984848;

    prim.lambda =  M / ( two * T * Ru );
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



#endif // USE_PASSIVE_SCALAR



#endif

