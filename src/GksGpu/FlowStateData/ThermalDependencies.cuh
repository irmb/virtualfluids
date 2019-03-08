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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_PASSIVE_SCALAR

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void getMolarMasses( real& M_O2, 
                                                real& M_N2, 
                                                real& M_CH4, 
                                                real& M_H2O, 
                                                real& M_CO2 )
{
    M_O2  = real(32.00e-3);    // kg / mol
    M_N2  = real(28.00e-3);    // kg / mol
    M_CH4 = real(16.00e-3);    // kg / mol
    M_H2O = real(18.00e-3);    // kg / mol
    M_CO2 = real(44.00e-3);    // kg / mol
}

__host__ __device__ inline real getRu()
{
    return real(8.31445984848);    // J / (mol K)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void getMassFractions( real Z1, real Z2, 
                                                  real& Y_O2, 
                                                  real& Y_N2, 
                                                  real& Y_CH4, 
                                                  real& Y_H2O, 
                                                  real& Y_CO2,
                                                  real& M )
{
    real M_O2, M_N2, M_CH4, M_H2O, M_CO2;

    getMolarMasses(M_O2, M_N2, M_CH4, M_H2O, M_CO2);

    const real Y_CH4_Inflow = real(1.0  );
    const real Y_N2_ambient = real(0.767);
    const real Y_O2_ambient = real(0.233);

    //////////////////////////////////////////////////////////////////////////

    const real Z = Z1 + Z2;
    
    Y_CH4 = Y_CH4_Inflow * Z1;
    Y_N2  = (one - Z) * Y_N2_ambient;
    Y_O2  = (one - Z) * Y_O2_ambient - two  * ( M_O2  / M_CH4 ) * Y_CH4_Inflow * Z2;
    Y_H2O =                            two  * ( M_H2O / M_CH4 ) * Y_CH4_Inflow * Z2;
    Y_CO2 =                                   ( M_CO2 / M_CH4 ) * Y_CH4_Inflow * Z2;
    
    M = one / ( Y_O2  / M_O2
              + Y_N2  / M_N2
              + Y_CH4 / M_CH4
              + Y_H2O / M_H2O
              + Y_CO2 / M_CO2 );
}

__host__ __device__ inline void getMassFractions( const ConservedVariables& cons, 
                                                  real& Y_O2, 
                                                  real& Y_N2, 
                                                  real& Y_CH4, 
                                                  real& Y_H2O, 
                                                  real& Y_CO2,
                                                  real& M )
{
    getMassFractions(cons.rhoS_1 / cons.rho,
                     cons.rhoS_2 / cons.rho,
                     Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                     M);
}

__host__ __device__ inline void getMassFractions( const PrimitiveVariables& prim, 
                                                  real& Y_O2, 
                                                  real& Y_N2, 
                                                  real& Y_CH4, 
                                                  real& Y_H2O, 
                                                  real& Y_CO2,
                                                  real& M )
{
    getMassFractions(prim.S_1,
                     prim.S_2,
                     Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                     M);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getMolarMass( const ConservedVariables& cons )
{
    real Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2;
    real M;

    getMassFractions(cons,
                     Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                     M);

    return M;
}

__host__ __device__ inline real getMolarMass( const PrimitiveVariables& prim )
{
    real Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2;
    real M;

    getMassFractions(prim,
                     Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, 
                     M);

    return M;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void getMoleFractions( real Z1, real Z2, 
                                                  real& X_O2, 
                                                  real& X_N2, 
                                                  real& X_CH4, 
                                                  real& X_H2O, 
                                                  real& X_CO2,
                                                  real& M )
{
    real Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2;
    real M_O2, M_N2, M_CH4, M_H2O, M_CO2;

    getMassFractions(Z1, Z2, Y_O2, Y_N2, Y_CH4, Y_H2O, Y_CO2, M);
    getMolarMasses  (        M_O2, M_N2, M_CH4, M_H2O, M_CO2 );

    X_O2  = Y_O2  * M / M_O2 ;
    X_N2  = Y_N2  * M / M_N2 ;
    X_CH4 = Y_CH4 * M / M_CH4;
    X_H2O = Y_H2O * M / M_H2O;
    X_CO2 = Y_CO2 * M / M_CO2;
}

__host__ __device__ inline void getMoleFractions( const ConservedVariables& cons, 
                                                  real& X_O2, 
                                                  real& X_N2, 
                                                  real& X_CH4, 
                                                  real& X_H2O, 
                                                  real& X_CO2,
                                                  real& M )
{
    getMoleFractions(cons.rhoS_1 / cons.rho,
                     cons.rhoS_2 / cons.rho,
                     X_O2, X_N2, X_CH4, X_H2O, X_CO2, 
                     M);
}

__host__ __device__ inline void getMoleFractions( const PrimitiveVariables& prim, 
                                                  real& X_O2, 
                                                  real& X_N2, 
                                                  real& X_CH4, 
                                                  real& X_H2O, 
                                                  real& X_CO2,
                                                  real& M )
{
    getMoleFractions(prim.S_1,
                     prim.S_2,
                     X_O2, X_N2, X_CH4, X_H2O, X_CO2, 
                     M);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getCp( const real T,
                                       const real X_O2, 
                                       const real X_N2, 
                                       const real X_CH4, 
                                       const real X_H2O, 
                                       const real X_CO2 )
{
    return X_O2  * getCpO2 (T)
         + X_N2  * getCpN2 (T)
         + X_CH4 * getCpCH4(T)
         + X_H2O * getCpH2O(T)
         + X_CO2 * getCpCO2(T);
}

__host__ __device__ inline real getIntegratedCv( const real T,
                                       const real X_O2, 
                                       const real X_N2, 
                                       const real X_CH4, 
                                       const real X_H2O, 
                                       const real X_CO2 )
{
    return X_O2  * getIntegratedCvO2 (T)
         + X_N2  * getIntegratedCvN2 (T)
         + X_CH4 * getIntegratedCvCH4(T)
         + X_H2O * getIntegratedCvH2O(T)
         + X_CO2 * getIntegratedCvCO2(T);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getK( const real T,
                                      const real X_O2, 
                                      const real X_N2, 
                                      const real X_CH4, 
                                      const real X_H2O, 
                                      const real X_CO2 )
{
    return two * getCp( T, X_O2, X_N2, X_CH4, X_H2O, X_CO2 ) / getRu() - five;
}

__host__ __device__ inline real getK( const ConservedVariables& cons )
{
    real X_O2, X_N2, X_CH4, X_H2O, X_CO2;
    real M;

    getMoleFractions(cons, X_O2, X_N2, X_CH4, X_H2O, X_CO2, M);

    real Eint = ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) / cons.rho; // J / kg

    //////////////////////////////////////////////////////////////////////////

    real TAitken[3];

    real T0 = real(600.0);
    real T  = T0;

    //////////////////////////////////////////////////////////////////////////

    uint nIter = 9;
    uint iAitkenStart = nIter - 3;

#pragma unroll
    for( uint i = 0; i < nIter; i++ )
    {
        real Cp = getCp( T, X_O2, X_N2, X_CH4, X_H2O, X_CO2 );
    
        T = c1o2 * ( T + ( Eint * M ) / ( Cp - getRu() ) );
        //T = ( Eint * M ) / ( Cp - Ru);

        if( i >= iAitkenStart ) TAitken[i - iAitkenStart] = T;
    }

    //if( fabs( T - T0 ) / T0 > real(0.01) )
    if( fabs( TAitken[2] - TAitken[1] ) / TAitken[1] > real(0.0001) )
    {
        T = ( TAitken[1] * TAitken[1] - TAitken[0] * TAitken[2] )
          / ( TAitken[1] + TAitken[1] - TAitken[0] - TAitken[2] );
    }

    //////////////////////////////////////////////////////////////////////////

    return getK(T, X_O2, X_N2, X_CH4, X_H2O, X_CO2);
}

//__host__ __device__ inline real getK( const PrimitiveVariables& prim )
//{
//    real X_O2, X_N2, X_CH4, X_H2O, X_CO2;
//    real M;
//
//    getMoleFractions(prim, X_O2, X_N2, X_CH4, X_H2O, X_CO2, M);
//
//    //////////////////////////////////////////////////////////////////////////
//
//    real T = M / ( two * prim.lambda * getRu() );
//
//    //////////////////////////////////////////////////////////////////////////
//
//    return getK(T, X_O2, X_N2, X_CH4, X_H2O, X_CO2);
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getT( const PrimitiveVariables& prim )
{
    real M = getMolarMass(prim);

    return M / ( two * prim.lambda * getRu() );
}

__host__ __device__ inline real getT( const ConservedVariables& cons )
{
    real X_O2, X_N2, X_CH4, X_H2O, X_CO2;
    real M;

    getMoleFractions(cons, X_O2, X_N2, X_CH4, X_H2O, X_CO2, M);

    real Eint = ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) / cons.rho; // J / kg

    //////////////////////////////////////////////////////////////////////////

    real T0 = real(1000.0);
    real T  = T0;

    //////////////////////////////////////////////////////////////////////////

    uint nIter = 10;

#pragma unroll
    for( uint i = 0; i < nIter; i++ )
    {
        real Cv = getCp( T, X_O2, X_N2, X_CH4, X_H2O, X_CO2 ) - getRu();
    
        real f = getIntegratedCv( T, X_O2, X_N2, X_CH4, X_H2O, X_CO2 ) - Eint * M;

        T = T - real(0.9) * f / Cv;
    }

    //////////////////////////////////////////////////////////////////////////

    return T;
}

__host__ __device__ inline real getEint( const PrimitiveVariables& prim )
{
    real X_O2, X_N2, X_CH4, X_H2O, X_CO2;
    real M;

    getMoleFractions(prim, X_O2, X_N2, X_CH4, X_H2O, X_CO2, M);

    //////////////////////////////////////////////////////////////////////////
    
    real T =  getT(prim);

    real Eint = getIntegratedCv( T, X_O2, X_N2, X_CH4, X_H2O, X_CO2 ) / M;

    //////////////////////////////////////////////////////////////////////////

    return Eint;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getlambda( const ConservedVariables& cons )
{
    real M = getMolarMass(cons);

    real T = getT(cons);

    return M / ( two * T * getRu() );
}

__host__ __device__ inline void setLambdaFromT( PrimitiveVariables& prim, real T )
{
    real M = getMolarMass(prim);

    prim.lambda =  M / ( two * T * getRu() );
}




#endif // USE_PASSIVE_SCALAR



#endif

