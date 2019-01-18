#ifndef ExpansionCoefficients_CUH
#define ExpansionCoefficients_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

__host__ __device__ inline void computeExpansionCoefficients(const PrimitiveVariables & facePrim, 
                                                             const ConservedVariables & gradient,
                                                             const real K, 
                                                             real expansionCoefficient[LENGTH_CELL_DATA])
{
    real two_E, 
             rho_dU_dx, 
             rho_dV_dx, 
             rho_dW_dx, 
         two_rho_dE_dx;

    two_E = facePrim.U * facePrim.U 
          + facePrim.V * facePrim.V 
          + facePrim.W * facePrim.W 
          + c1o2 * ( K + three ) / facePrim.lambda;

    rho_dU_dx     =       gradient.rhoU - facePrim.U  * gradient.rho;
    rho_dV_dx     =       gradient.rhoV - facePrim.V  * gradient.rho;
    rho_dW_dx     =       gradient.rhoW - facePrim.W  * gradient.rho;
    two_rho_dE_dx = two * gradient.rhoE -      two_E  * gradient.rho;

    expansionCoefficient[4] = ( four * facePrim.lambda * facePrim.lambda ) / ( K + three )
                            * ( two_rho_dE_dx - two * facePrim.U * rho_dU_dx 
                                              - two * facePrim.V * rho_dV_dx 
                                              - two * facePrim.W * rho_dW_dx );

    expansionCoefficient[3] = two * facePrim.lambda * rho_dW_dx - facePrim.W * expansionCoefficient[4];

    expansionCoefficient[2] = two * facePrim.lambda * rho_dV_dx - facePrim.V * expansionCoefficient[4];

    expansionCoefficient[1] = two * facePrim.lambda * rho_dU_dx - facePrim.U * expansionCoefficient[4];

    expansionCoefficient[0] = gradient.rho -   facePrim.U * expansionCoefficient[1] 
                                           -   facePrim.V * expansionCoefficient[2] 
                                           -   facePrim.W * expansionCoefficient[3] 
                                           - c1o2 * two_E * expansionCoefficient[4];

#ifdef USE_PASSIVE_SCALAR
	expansionCoefficient[5] = two * facePrim.lambda * (gradient.rhoS_1 - facePrim.S_1 * gradient.rho);
	expansionCoefficient[6] = two * facePrim.lambda * (gradient.rhoS_2 - facePrim.S_2 * gradient.rho);
#endif // USE_PASSIVE_SCALAR
}

#endif