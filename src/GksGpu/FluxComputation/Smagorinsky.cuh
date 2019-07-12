#ifndef Smagorinsky_CUH
#define Smagorinsky_CUH

#include <cmath>

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

inline __host__ __device__ real getTurbulentViscositySmagorinsky(const Parameters & parameters, 
                                                                 const PrimitiveVariables& facePrim, 
                                                                 const ConservedVariables gradX1, 
                                                                 const ConservedVariables gradX2, 
                                                                 const ConservedVariables gradX3 )
{
    // See FDS 6 Technical Reference Guide, Section 4.2.8

    real dUdx1 = ( gradX1.rhoU - facePrim.U * gradX1.rho )/* / facePrim.rho*/;
    real dUdx2 = ( gradX2.rhoU - facePrim.U * gradX2.rho )/* / facePrim.rho*/;
    real dUdx3 = ( gradX3.rhoU - facePrim.U * gradX3.rho )/* / facePrim.rho*/;
    real dVdx1 = ( gradX1.rhoV - facePrim.V * gradX1.rho )/* / facePrim.rho*/;
    real dVdx2 = ( gradX2.rhoV - facePrim.V * gradX2.rho )/* / facePrim.rho*/;
    real dVdx3 = ( gradX3.rhoV - facePrim.V * gradX3.rho )/* / facePrim.rho*/;
    real dWdx1 = ( gradX1.rhoW - facePrim.W * gradX1.rho )/* / facePrim.rho*/;
    real dWdx2 = ( gradX2.rhoW - facePrim.W * gradX2.rho )/* / facePrim.rho*/;
    real dWdx3 = ( gradX3.rhoW - facePrim.W * gradX3.rho )/* / facePrim.rho*/;

    real S11sq = dUdx1*dUdx1;
    real S22sq = dVdx2*dVdx2;
    real S33sq = dWdx3*dWdx3;

    real S12sq = c1o4 * ( dUdx2 + dVdx1 ) * ( dUdx2 + dVdx1 );
    real S13sq = c1o4 * ( dUdx3 + dWdx1 ) * ( dUdx3 + dWdx1 );
    real S23sq = c1o4 * ( dVdx3 + dWdx2 ) * ( dVdx3 + dWdx2 );

    real divergence = dUdx1 + dVdx2 + dWdx3;

    real S = sqrt( c2o1 * ( S11sq + S22sq + S33sq + c2o1 * ( S12sq + S13sq + S23sq ) ) - c2o3 * divergence * divergence );

    real Cs = parameters.smagorinskyConstant;

    return facePrim.rho * Cs*Cs * parameters.dx*parameters.dx * S;
}

#endif
