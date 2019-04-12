#ifndef ThermalDependencies_H
#define ThermalDependencies_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

#include <math.h>

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

#define R_U real(8.31445984848)

#define CP_INT_FIT_T0  real(300)
#define CP_INT_FIT_A_A real(0.00161062666357469)
#define CP_INT_FIT_A_B real(-18613.4961810207)
#define CP_INT_FIT_F_A real(0.011624194441801)
#define CP_INT_FIT_F_B real(-4131.10340570215)
#define CP_INT_FIT_P_A real(0.00232864732198785)
#define CP_INT_FIT_P_B real(-13716.8715588326)

#define CP_FIT_A_A real(3.5368e-10)
#define CP_FIT_A_B real(-3.1248e-06)
#define CP_FIT_A_C real(0.010139)
#define CP_FIT_A_D real(25.7644)
#define CP_FIT_F_A real(3.7886e-09)
#define CP_FIT_F_B real(-3.098e-05)
#define CP_FIT_F_C real(0.089563)
#define CP_FIT_F_D real(9.1489)
#define CP_FIT_P_A real(4.601e-10)
#define CP_FIT_P_B real(-4.2602e-06)
#define CP_FIT_P_C real(0.014334)
#define CP_FIT_P_D real(25.8835)

#define M_A real(0.02884)
#define M_P real(0.0276199095022624)
#define M_F real(0.016)

#define T_FAKTOR real(1.0)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getMolarMass( const ConservedVariables& cons )
{
    real Y_F = cons.rhoS_1 / cons.rho;
    real Y_P = cons.rhoS_2 / cons.rho;

    real Y_A = one - Y_F - Y_P;

    real M = one / ( Y_A / M_A
                   + Y_F / M_F
                   + Y_P / M_P );

    return M;
}

__host__ __device__ inline real getMolarMass( const PrimitiveVariables& prim )
{
    
    real Y_F = prim.S_1;
    real Y_P = prim.S_2;

    real Y_A = one - Y_F - Y_P;

    real M = one / ( Y_A / M_A
                   + Y_F / M_F
                   + Y_P / M_P );

    return M;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void getMoleFractions( real Z1, real Z2, 
                                                  real& X_A, 
                                                  real& X_F, 
                                                  real& X_P,
                                                  real& M )
{
    real Y_F = Z1;
    real Y_P = Z2;

    real Y_A = one - Y_F - Y_P;

    M = one / ( Y_A / M_A
              + Y_F / M_F
              + Y_P / M_P );

    X_A = Y_A * M / M_A;
    X_F = Y_F * M / M_F;
    X_P = Y_P * M / M_P;
}

__host__ __device__ inline void getMoleFractions( const ConservedVariables& cons, 
                                                  real& X_A, 
                                                  real& X_F, 
                                                  real& X_P,
                                                  real& M )
{
    getMoleFractions(cons.rhoS_1 / cons.rho,
                     cons.rhoS_2 / cons.rho,
                     X_A, X_F, X_P, 
                     M);
}

__host__ __device__ inline void getMoleFractions( const PrimitiveVariables& prim, 
                                                  real& X_A, 
                                                  real& X_F, 
                                                  real& X_P,
                                                  real& M )
{
    getMoleFractions(prim.S_1,
                     prim.S_2,
                     X_A, X_F, X_P, 
                     M);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getCp( real T,
                                       const real X_A, 
                                       const real X_F, 
                                       const real X_P )
{
    T *= T_FAKTOR;

    return X_A * ( CP_FIT_A_A * T * T * T   +   CP_FIT_A_B * T * T   +   CP_FIT_A_C * T   +   CP_FIT_A_D )
         + X_F * ( CP_FIT_F_A * T * T * T   +   CP_FIT_F_B * T * T   +   CP_FIT_F_C * T   +   CP_FIT_F_D )
         + X_P * ( CP_FIT_P_A * T * T * T   +   CP_FIT_P_B * T * T   +   CP_FIT_P_C * T   +   CP_FIT_P_D );
}

__host__ __device__ inline real getIntegratedCv( real T,
                                                 const real X_A, 
                                                 const real X_F, 
                                                 const real X_P )
{
    T *= T_FAKTOR;

    return X_A * ( CP_INT_FIT_A_A * ( T - CP_INT_FIT_T0 ) * ( T - CP_INT_FIT_A_B ) )
         + X_F * ( CP_INT_FIT_F_A * ( T - CP_INT_FIT_T0 ) * ( T - CP_INT_FIT_F_B ) )
         + X_P * ( CP_INT_FIT_P_A * ( T - CP_INT_FIT_T0 ) * ( T - CP_INT_FIT_P_B ) )
         - ( T - CP_INT_FIT_T0 ) * R_U;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getT( const PrimitiveVariables& prim )
{
    real M = getMolarMass(prim);

    real T = M / ( two * prim.lambda * R_U );

    return T;
}

__host__ __device__ inline real getT( const ConservedVariables& cons )
{
    real X_A, X_F, X_P;
    real M;

    getMoleFractions(cons, X_A, X_F, X_P, M);

    real Eint = ( cons.rhoE - c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho ) / cons.rho; // J / kg

    Eint -= real(100000.0);

    real a = X_A * CP_INT_FIT_A_A
           + X_P * CP_INT_FIT_P_A
           + X_F * CP_INT_FIT_F_A;

    real mp = ( X_A * CP_INT_FIT_A_A * ( CP_INT_FIT_T0 + CP_INT_FIT_A_B )
              + X_P * CP_INT_FIT_P_A * ( CP_INT_FIT_T0 + CP_INT_FIT_P_B )
              + X_F * CP_INT_FIT_F_A * ( CP_INT_FIT_T0 + CP_INT_FIT_F_B ) 
              + R_U ) / a;

    real q  = ( X_A * CP_INT_FIT_T0 * CP_INT_FIT_A_A * CP_INT_FIT_A_B
              + X_P * CP_INT_FIT_T0 * CP_INT_FIT_P_A * CP_INT_FIT_P_B
              + X_F * CP_INT_FIT_T0 * CP_INT_FIT_F_A * CP_INT_FIT_F_B
              - Eint * M + CP_INT_FIT_T0 * R_U ) / a;

    real T = c1o2 * mp + sqrt( c1o4 * mp * mp - q );

    T /= T_FAKTOR;

    return T;
}

__host__ __device__ inline real getEint( const PrimitiveVariables& prim )
{
    real X_A, X_F, X_P;
    real M;

    getMoleFractions(prim, X_A, X_F, X_P, M);

    real T = getT(prim);

    real Eint = getIntegratedCv( T, X_A, X_F, X_P ) / M;

    Eint += real(100000.0);

    return Eint;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getK( const real T,
                                      const real X_A, 
                                      const real X_F, 
                                      const real X_P )
{
    return two * getCp( T, X_A, X_F, X_P ) / R_U - five;
}

__host__ __device__ inline real getK( const ConservedVariables& cons )
{
    real T = getT(cons);

    real X_A, X_F, X_P;
    real M;

    getMoleFractions(cons, X_A, X_F, X_P, M);

    return getK(T, X_A, X_F, X_P);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline real getlambda( const ConservedVariables& cons )
{
    real M = getMolarMass(cons);

    real T = getT(cons);

    return M / ( two * T * R_U );
}

__host__ __device__ inline void setLambdaFromT( PrimitiveVariables& prim, real T )
{
    real M = getMolarMass(prim);

    prim.lambda =  M / ( two * T * R_U );
}

#endif // USE_PASSIVE_SCALAR



#endif

