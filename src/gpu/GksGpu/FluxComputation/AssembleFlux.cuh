//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file AssembleFlux.cuh
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#ifndef AssembleFlux_CUH
#define AssembleFlux_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "FluxComputation/SutherlandsLaw.cuh"
#include "FluxComputation/Moments.cuh"

extern __device__ real atomicAdd(real* address, real val);

#define NUMBER_OF_MOMENTS 7

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//! \brief Computes the time derivative of the conserved variables based on 
//!  expansion coefficients and moments of the equilibrium distribution
//!
//! This function implements the expanded equations found in Eq. (17) in 
//! <a href="https://doi.org/10.1016/j.ijthermalsci.2018.10.004"><b>[ Lenz et al. (2019), DOI: 10.1016/j.ijthermalsci.2018.10.004 ]</b></a>.
//! which is based on 
//! <a href="https://doi.org/10.1006/jcph.2001.6790"><b>[ Kun Xu (2001), DOI: 10.1006/jcph.2001.6790 ]</b></a>
//! and 
//! <a href="https://doi.org/10.1016/j.jcp.2007.06.024"><b>[ Tian et al. (2007), DOI: 10.1016/j.jcp.2007.06.024 ]</b></a>.
//!
//! Part of the equations was expanded using a computer algebra system 
//! and are, hence, not well structured.
//! 
//! \param[in]  facePrim       flow state on the interface as \ref PrimitiveVariables
//! \param[in]  momentU        array of moments of the equilibrium function with respect to x-velocity
//! \param[in]  momentV        array of moments of the equilibrium function with respect to y-velocity
//! \param[in]  momentW        array of moments of the equilibrium function with respect to z-velocity
//! \param[in]  momentXi       array of moments of the equilibrium function with respect to internal degrees of freedom
//! \param[in]  ax             array of spatial expansion coefficients of the Taylor expansion of the equilibrium distribution
//! \param[in]  ay             array of spatial expansion coefficients of the Taylor expansion of the equilibrium distribution 
//! \param[in]  az             array of spatial expansion coefficients of the Taylor expansion of the equilibrium distribution
//! \param[in]  force          volume force vector
//! \param[out] timeGrad       time derivatives of the \ref ConservedVariables
__host__ __device__ inline void computeTimeDerivative( const PrimitiveVariables& facePrim, 
                                                       const real momentU [ NUMBER_OF_MOMENTS ], 
                                                       const real momentV [ NUMBER_OF_MOMENTS ], 
                                                       const real momentW [ NUMBER_OF_MOMENTS ], 
                                                       const real momentXi[ NUMBER_OF_MOMENTS ],
                                                       const real ax[LENGTH_CELL_DATA],
                                                       const real ay[LENGTH_CELL_DATA],
                                                       const real az[LENGTH_CELL_DATA],
                                                       const Vec3& force,
                                                       ConservedVariables& timeGrad )
{
    timeGrad.rho = ax[0]*momentU[1] + ax[1]*momentU[2] + c1o2*ax[4]*momentU[3] + ay[0]*momentV[1] + 
   ax[2]*momentU[1]*momentV[1] + ay[1]*momentU[1]*momentV[1] + 
   c1o2*ay[4]*momentU[2]*momentV[1] + ay[2]*momentV[2] + c1o2*ax[4]*momentU[1]*momentV[2] + 
   c1o2*ay[4]*momentV[3] + az[0]*momentW[1] + ax[3]*momentU[1]*momentW[1] + 
   az[1]*momentU[1]*momentW[1] + c1o2*az[4]*momentU[2]*momentW[1] + 
   ay[3]*momentV[1]*momentW[1] + az[2]*momentV[1]*momentW[1] + 
   c1o2*az[4]*momentV[2]*momentW[1] + az[3]*momentW[2] + c1o2*ax[4]*momentU[1]*momentW[2] + 
   c1o2*ay[4]*momentV[1]*momentW[2] + c1o2*az[4]*momentW[3] + 
   c1o2*ax[4]*momentU[1]*momentXi[2] + c1o2*ay[4]*momentV[1]*momentXi[2] + 
   c1o2*az[4]*momentW[1]*momentXi[2];

    timeGrad.rhoU = ax[0]*momentU[2] + ax[1]*momentU[3] + c1o2*ax[4]*momentU[4] + 
   ay[0]*momentU[1]*momentV[1] + ax[2]*momentU[2]*momentV[1] + 
   ay[1]*momentU[2]*momentV[1] + c1o2*ay[4]*momentU[3]*momentV[1] + 
   ay[2]*momentU[1]*momentV[2] + c1o2*ax[4]*momentU[2]*momentV[2] + 
   c1o2*ay[4]*momentU[1]*momentV[3] + az[0]*momentU[1]*momentW[1] + 
   ax[3]*momentU[2]*momentW[1] + az[1]*momentU[2]*momentW[1] + 
   c1o2*az[4]*momentU[3]*momentW[1] + ay[3]*momentU[1]*momentV[1]*momentW[1] + 
   az[2]*momentU[1]*momentV[1]*momentW[1] + c1o2*az[4]*momentU[1]*momentV[2]*momentW[1] + 
   az[3]*momentU[1]*momentW[2] + c1o2*ax[4]*momentU[2]*momentW[2] + 
   c1o2*ay[4]*momentU[1]*momentV[1]*momentW[2] + c1o2*az[4]*momentU[1]*momentW[3] + 
   c1o2*ax[4]*momentU[2]*momentXi[2] + c1o2*ay[4]*momentU[1]*momentV[1]*momentXi[2] + 
   c1o2*az[4]*momentU[1]*momentW[1]*momentXi[2];

    timeGrad.rhoV = ax[0]*momentU[1]*momentV[1] + ax[1]*momentU[2]*momentV[1] + 
   c1o2*ax[4]*momentU[3]*momentV[1] + ay[0]*momentV[2] + ax[2]*momentU[1]*momentV[2] + 
   ay[1]*momentU[1]*momentV[2] + c1o2*ay[4]*momentU[2]*momentV[2] + ay[2]*momentV[3] + 
   c1o2*ax[4]*momentU[1]*momentV[3] + c1o2*ay[4]*momentV[4] + az[0]*momentV[1]*momentW[1] + 
   ax[3]*momentU[1]*momentV[1]*momentW[1] + az[1]*momentU[1]*momentV[1]*momentW[1] + 
   c1o2*az[4]*momentU[2]*momentV[1]*momentW[1] + ay[3]*momentV[2]*momentW[1] + 
   az[2]*momentV[2]*momentW[1] + c1o2*az[4]*momentV[3]*momentW[1] + 
   az[3]*momentV[1]*momentW[2] + c1o2*ax[4]*momentU[1]*momentV[1]*momentW[2] + 
   c1o2*ay[4]*momentV[2]*momentW[2] + c1o2*az[4]*momentV[1]*momentW[3] + 
   c1o2*ax[4]*momentU[1]*momentV[1]*momentXi[2] + c1o2*ay[4]*momentV[2]*momentXi[2] + 
   c1o2*az[4]*momentV[1]*momentW[1]*momentXi[2];

    timeGrad.rhoW = ax[0]*momentU[1]*momentW[1] + ax[1]*momentU[2]*momentW[1] + 
   c1o2*ax[4]*momentU[3]*momentW[1] + ay[0]*momentV[1]*momentW[1] + 
   ax[2]*momentU[1]*momentV[1]*momentW[1] + ay[1]*momentU[1]*momentV[1]*momentW[1] + 
   c1o2*ay[4]*momentU[2]*momentV[1]*momentW[1] + ay[2]*momentV[2]*momentW[1] + 
   c1o2*ax[4]*momentU[1]*momentV[2]*momentW[1] + c1o2*ay[4]*momentV[3]*momentW[1] + 
   az[0]*momentW[2] + ax[3]*momentU[1]*momentW[2] + az[1]*momentU[1]*momentW[2] + 
   c1o2*az[4]*momentU[2]*momentW[2] + ay[3]*momentV[1]*momentW[2] + 
   az[2]*momentV[1]*momentW[2] + c1o2*az[4]*momentV[2]*momentW[2] + az[3]*momentW[3] + 
   c1o2*ax[4]*momentU[1]*momentW[3] + c1o2*ay[4]*momentV[1]*momentW[3] + 
   c1o2*az[4]*momentW[4] + c1o2*ax[4]*momentU[1]*momentW[1]*momentXi[2] + 
   c1o2*ay[4]*momentV[1]*momentW[1]*momentXi[2] + c1o2*az[4]*momentW[2]*momentXi[2];

    timeGrad.rhoE = c1o2*ax[0]*momentU[3] + c1o2*ax[1]*momentU[4] + c1o4*ax[4]*momentU[5] + 
   c1o2*ay[0]*momentU[2]*momentV[1] + c1o2*ax[2]*momentU[3]*momentV[1] + 
   c1o2*ay[1]*momentU[3]*momentV[1] + c1o4*ay[4]*momentU[4]*momentV[1] + 
   c1o2*ax[0]*momentU[1]*momentV[2] + c1o2*ax[1]*momentU[2]*momentV[2] + 
   c1o2*ay[2]*momentU[2]*momentV[2] + c1o2*ax[4]*momentU[3]*momentV[2] + 
   c1o2*ay[0]*momentV[3] + c1o2*ax[2]*momentU[1]*momentV[3] + 
   c1o2*ay[1]*momentU[1]*momentV[3] + c1o2*ay[4]*momentU[2]*momentV[3] + 
   c1o2*ay[2]*momentV[4] + c1o4*ax[4]*momentU[1]*momentV[4] + c1o4*ay[4]*momentV[5] + 
   c1o2*az[0]*momentU[2]*momentW[1] + c1o2*ax[3]*momentU[3]*momentW[1] + 
   c1o2*az[1]*momentU[3]*momentW[1] + c1o4*az[4]*momentU[4]*momentW[1] + 
   c1o2*ay[3]*momentU[2]*momentV[1]*momentW[1] + c1o2*az[2]*momentU[2]*momentV[1]*momentW[1] + 
   c1o2*az[0]*momentV[2]*momentW[1] + c1o2*ax[3]*momentU[1]*momentV[2]*momentW[1] + 
   c1o2*az[1]*momentU[1]*momentV[2]*momentW[1] + c1o2*az[4]*momentU[2]*momentV[2]*momentW[1] + 
   c1o2*ay[3]*momentV[3]*momentW[1] + c1o2*az[2]*momentV[3]*momentW[1] + 
   c1o4*az[4]*momentV[4]*momentW[1] + c1o2*ax[0]*momentU[1]*momentW[2] + 
   c1o2*ax[1]*momentU[2]*momentW[2] + c1o2*az[3]*momentU[2]*momentW[2] + 
   c1o2*ax[4]*momentU[3]*momentW[2] + c1o2*ay[0]*momentV[1]*momentW[2] + 
   c1o2*ax[2]*momentU[1]*momentV[1]*momentW[2] + c1o2*ay[1]*momentU[1]*momentV[1]*momentW[2] + 
   c1o2*ay[4]*momentU[2]*momentV[1]*momentW[2] + c1o2*ay[2]*momentV[2]*momentW[2] + 
   c1o2*az[3]*momentV[2]*momentW[2] + c1o2*ax[4]*momentU[1]*momentV[2]*momentW[2] + 
   c1o2*ay[4]*momentV[3]*momentW[2] + c1o2*az[0]*momentW[3] + 
   c1o2*ax[3]*momentU[1]*momentW[3] + c1o2*az[1]*momentU[1]*momentW[3] + 
   c1o2*az[4]*momentU[2]*momentW[3] + c1o2*ay[3]*momentV[1]*momentW[3] + 
   c1o2*az[2]*momentV[1]*momentW[3] + c1o2*az[4]*momentV[2]*momentW[3] + 
   c1o2*az[3]*momentW[4] + c1o4*ax[4]*momentU[1]*momentW[4] + 
   c1o4*ay[4]*momentV[1]*momentW[4] + c1o4*az[4]*momentW[5] + 
   c1o2*ax[0]*momentU[1]*momentXi[2] + c1o2*ax[1]*momentU[2]*momentXi[2] + 
   c1o2*ax[4]*momentU[3]*momentXi[2] + c1o2*ay[0]*momentV[1]*momentXi[2] + 
   c1o2*ax[2]*momentU[1]*momentV[1]*momentXi[2] + 
   c1o2*ay[1]*momentU[1]*momentV[1]*momentXi[2] + 
   c1o2*ay[4]*momentU[2]*momentV[1]*momentXi[2] + c1o2*ay[2]*momentV[2]*momentXi[2] + 
   c1o2*ax[4]*momentU[1]*momentV[2]*momentXi[2] + c1o2*ay[4]*momentV[3]*momentXi[2] + 
   c1o2*az[0]*momentW[1]*momentXi[2] + c1o2*ax[3]*momentU[1]*momentW[1]*momentXi[2] + 
   c1o2*az[1]*momentU[1]*momentW[1]*momentXi[2] + 
   c1o2*az[4]*momentU[2]*momentW[1]*momentXi[2] + 
   c1o2*ay[3]*momentV[1]*momentW[1]*momentXi[2] + 
   c1o2*az[2]*momentV[1]*momentW[1]*momentXi[2] + 
   c1o2*az[4]*momentV[2]*momentW[1]*momentXi[2] + c1o2*az[3]*momentW[2]*momentXi[2] + 
   c1o2*ax[4]*momentU[1]*momentW[2]*momentXi[2] + 
   c1o2*ay[4]*momentV[1]*momentW[2]*momentXi[2] + c1o2*az[4]*momentW[3]*momentXi[2] + 
   c1o4*ax[4]*momentU[1]*momentXi[4] + c1o4*ay[4]*momentV[1]*momentXi[4] + 
   c1o4*az[4]*momentW[1]*momentXi[4];

    //////////////////////////////////////////////////////////////////////////

    timeGrad.rho  += c2o1 * facePrim.lambda * (                                                          facePrim.U - momentU[1]                           ) * force.x
                   + c2o1 * facePrim.lambda * (                                                          facePrim.V -              momentV[1]              ) * force.y
                   + c2o1 * facePrim.lambda * (                                                          facePrim.W -                           momentW[1] ) * force.z ;
                                                                                                         
    timeGrad.rhoU += c2o1 * facePrim.lambda * (   momentU[1] *                                           facePrim.U - momentU[2]                           ) * force.x
                   + c2o1 * facePrim.lambda * (   momentU[1] *                                           facePrim.V - momentU[1] * momentV[1]              ) * force.y
                   + c2o1 * facePrim.lambda * (   momentU[1] *                                           facePrim.W - momentU[1] *              momentW[1] ) * force.z ;
                                                                                                         
    timeGrad.rhoV += c2o1 * facePrim.lambda * (                momentV[1] *                              facePrim.U - momentU[1] * momentV[1]              ) * force.x
                   + c2o1 * facePrim.lambda * (                momentV[1] *                              facePrim.V -              momentV[2]              ) * force.y
                   + c2o1 * facePrim.lambda * (                momentV[1] *                              facePrim.W -              momentV[1] * momentW[1] ) * force.z ;
                                                                                                         
    timeGrad.rhoW += c2o1 * facePrim.lambda * (                             momentW[1] *                 facePrim.U - momentU[1] *              momentW[1] ) * force.x
                   + c2o1 * facePrim.lambda * (                             momentW[1] *                 facePrim.V -              momentV[1] * momentW[1] ) * force.y
                   + c2o1 * facePrim.lambda * (                             momentW[1] *                 facePrim.W -                           momentW[2] ) * force.z ;

    timeGrad.rhoE +=       facePrim.lambda * ( ( momentU[2] + momentV[2] + momentW[2] + momentXi[2] ) * facePrim.U

                                             - ( momentU[3]
                                               + momentU[1] * momentV[2] 
                                               + momentU[1] *              momentW[2] 
                                               + momentU[1] *                           momentXi[2] )
                                             ) * force.x

                   +       facePrim.lambda * ( ( momentU[2] + momentV[2] + momentW[2] + momentXi[2] ) * facePrim.V

                                             - ( momentU[2] * momentV[1]
                                               +              momentV[3]
                                               +              momentV[1] * momentW[2]
                                               +              momentV[1] *              momentXi[2] )
                                             ) * force.y

                   +       facePrim.lambda * ( ( momentU[2] + momentV[2] + momentW[2] + momentXi[2] ) * facePrim.W

                                             - ( momentU[2] *              momentW[1]
                                               +              momentV[2] * momentW[1]
                                               +                           momentW[3]
                                               +                           momentW[1] * momentXi[2] )
                                             ) * force.z ;

    //////////////////////////////////////////////////////////////////////////

    timeGrad.rho  *= - c1o1;
    timeGrad.rhoU *= - c1o1;
    timeGrad.rhoV *= - c1o1;
    timeGrad.rhoW *= - c1o1;
    timeGrad.rhoE *= - c1o1;

    //////////////////////////////////////////////////////////////////////////

#ifdef USE_PASSIVE_SCALAR
	timeGrad.rhoS_1 = timeGrad.rho * facePrim.S_1
		            + ( ax[5] * momentU[1]  
		              + ay[5] *              momentV[1]
		              + az[5] *                           momentW[1] )
		            / (c1o2 * facePrim.lambda);

	timeGrad.rhoS_2 = timeGrad.rho * facePrim.S_2
		            + ( ax[6] * momentU[1]  
		              + ay[6] *              momentV[1]
		              + az[6] *                           momentW[1] )
		            / (c1o2 * facePrim.lambda);

    timeGrad.rhoS_1 *= - c1o1;
    timeGrad.rhoS_2 *= - c1o1;
#endif // USE_PASSIVE_SCALAR

    //////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Computes the products of relaxation time \f$\tau\f$ and time step 
//! \f$\Delta t\f$ that stem from the symbolic time integration of the flux
//!
//! It further considers temperature dependent viscisity based on 
//! \ref sutherlandsLaw, if enabled.
//!
//! \param[in]  facePrim          flow state on the interface as \ref PrimitiveVariables
//! \param[in]  parameters        \ref Parameters struct
//! \param[out] timeCoefficients  various products of relaxation time \f$\tau\f$ and time step \f$\Delta t\f$ that stem from the symbolic time integration of the flux
__host__ __device__ inline void computeTimeCoefficients(const PrimitiveVariables & facePrim, const Parameters& parameters, real timeCoefficients[4])
{
    real r   = parameters.lambdaRef / facePrim.lambda;

    real mu;
    if ( parameters.viscosityModel == ViscosityModel::constant ){
        mu = parameters.mu;
    }
    else if ( parameters.viscosityModel == ViscosityModel::sutherlandsLaw ){
        mu = sutherlandsLaw( parameters, r );
    }

    real tau = c2o1 * facePrim.lambda * mu / facePrim.rho;

    real dt = parameters.dt;

    timeCoefficients[0] =                         dt;
    timeCoefficients[1] =                 - tau * dt;
    timeCoefficients[2] =  c1o2 * dt * dt - tau * dt;

    timeCoefficients[3] =                   tau     ;
}

//! \brief Computes the relaxation time
//!
//! It further considers temperature dependent viscisity based on 
//! \ref sutherlandsLaw, if enabled.
//!
//! \param[in]  facePrim          flow state on the interface as \ref PrimitiveVariables
//! \param[in]  parameters        \ref Parameters struct
//! \param[out] tau               relaxation time
__host__ __device__ inline void getTau(const PrimitiveVariables & facePrim, const Parameters& parameters, real& tau)
{
    real r   = parameters.lambdaRef / facePrim.lambda;

    real mu;
    if ( parameters.viscosityModel == ViscosityModel::constant ){
        mu = parameters.mu;
    }
    else if ( parameters.viscosityModel == ViscosityModel::sutherlandsLaw ){
        mu = sutherlandsLaw( parameters, r );
    }

    tau = c2o1 * facePrim.lambda * mu / facePrim.rho;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief Computes the flux of the conserved variables based on expansion
//! coefficients and moments of the equilibrium distribution
//!
//! 
//! This function implements the expanded equations found in Eq. (19) in 
//! <a href="https://doi.org/10.1016/j.ijthermalsci.2018.10.004"><b>[ Lenz et al. (2019), DOI: 10.1016/j.ijthermalsci.2018.10.004 ]</b></a>.
//! which is based on 
//! <a href="https://doi.org/10.1006/jcph.2001.6790"><b>[ Kun Xu (2001), DOI: 10.1006/jcph.2001.6790 ]</b></a>
//! and 
//! <a href="https://doi.org/10.1016/j.jcp.2007.06.024"><b>[ Tian et al. (2007), DOI: 10.1016/j.jcp.2007.06.024 ]</b></a>.
//!
//! This function incorporates the Pandtl number fix from Eq. Eq. (27) in 
//! <a href="https://doi.org/10.1006/jcph.2001.6790"><b>[ Kun Xu (2001), DOI: 10.1006/jcph.2001.6790 ]</b></a>.
//!
//! Part of the equations was expanded using a computer algebra system 
//! and are, hence, not well structured.
//!
//! \param[in]  facePrim          flow state on the interface as \ref PrimitiveVariables
//! \param[in]  momentU           array of moments of the equilibrium function with respect to x-velocity
//! \param[in]  momentV           array of moments of the equilibrium function with respect to y-velocity
//! \param[in]  momentW           array of moments of the equilibrium function with respect to z-velocity
//! \param[in]  momentXi          array of moments of the equilibrium function with respect to internal degrees of freedom
//! \param[in]  ax                array of spatial expansion coefficients of the Taylor expansion of the equilibrium distribution
//! \param[in]  ay                array of spatial expansion coefficients of the Taylor expansion of the equilibrium distribution 
//! \param[in]  az                array of spatial expansion coefficients of the Taylor expansion of the equilibrium distribution
//! \param[in]  at                array of temporal expansion coefficients of the Taylor expansion of the equilibrium distribution
//! \param[in]  timeCoefficients  various products of relaxation time \f$\tau\f$ and time step \f$\Delta t\f$ that stem from the symbolic time integration of the flux
//! \param[in]  parameters        \ref Parameters struct
//! \param[in]  force             volume force vector
//! \param[out] flux              time integrated absolute amount of \ref ConservedVariables transported over the face in this time step
//! \param[out] heatFlux          amount of heat transported over the face
__host__ __device__ inline void assembleFlux( const PrimitiveVariables& facePrim, 
                                              const real momentU [ NUMBER_OF_MOMENTS ], 
                                              const real momentV [ NUMBER_OF_MOMENTS ], 
                                              const real momentW [ NUMBER_OF_MOMENTS ], 
                                              const real momentXi[ NUMBER_OF_MOMENTS ],
                                              const real ax[LENGTH_CELL_DATA],
                                              const real ay[LENGTH_CELL_DATA],
                                              const real az[LENGTH_CELL_DATA],
                                              const real at[LENGTH_CELL_DATA],
                                              const real timeCoefficients[4],
                                              const Parameters& parameters,
                                              const Vec3 force,
                                              ConservedVariables& flux,
                                              real& heatFlux )
{
    ConservedVariables flux_1, flux_2, flux_3;

    flux_1.rho  =           momentU[0+1]                          ;
    flux_1.rhoU =           momentU[1+1]                          ;
    flux_1.rhoV =           momentU[0+1] * momentV[1]             ;
    flux_1.rhoW =           momentU[0+1] *              momentW[1];

    flux_1.rhoE =  c1o2 * ( momentU[2+1]             
                          + momentU[0+1] * momentV[2]
                          + momentU[0+1] *              momentW[2]
                          + momentU[0+1] *                           momentXi[2] );

    //////////////////////////////////////////////////////////////////////////

    flux_2.rho  = ax[0]*momentU[2] + ax[1]*momentU[3] + c1o2*ax[4]*momentU[4] + 
   ay[0]*momentU[1]*momentV[1] + ax[2]*momentU[2]*momentV[1] + 
   ay[1]*momentU[2]*momentV[1] + c1o2*ay[4]*momentU[3]*momentV[1] + 
   ay[2]*momentU[1]*momentV[2] + c1o2*ax[4]*momentU[2]*momentV[2] + 
   c1o2*ay[4]*momentU[1]*momentV[3] + az[0]*momentU[1]*momentW[1] + 
   ax[3]*momentU[2]*momentW[1] + az[1]*momentU[2]*momentW[1] + 
   c1o2*az[4]*momentU[3]*momentW[1] + ay[3]*momentU[1]*momentV[1]*momentW[1] + 
   az[2]*momentU[1]*momentV[1]*momentW[1] + c1o2*az[4]*momentU[1]*momentV[2]*momentW[1] + 
   az[3]*momentU[1]*momentW[2] + c1o2*ax[4]*momentU[2]*momentW[2] + 
   c1o2*ay[4]*momentU[1]*momentV[1]*momentW[2] + c1o2*az[4]*momentU[1]*momentW[3] + 
   c1o2*ax[4]*momentU[2]*momentXi[2] + c1o2*ay[4]*momentU[1]*momentV[1]*momentXi[2] + 
   c1o2*az[4]*momentU[1]*momentW[1]*momentXi[2];

    flux_2.rhoU = ax[0]*momentU[3] + ax[1]*momentU[4] + c1o2*ax[4]*momentU[5] + 
   ay[0]*momentU[2]*momentV[1] + ax[2]*momentU[3]*momentV[1] + 
   ay[1]*momentU[3]*momentV[1] + c1o2*ay[4]*momentU[4]*momentV[1] + 
   ay[2]*momentU[2]*momentV[2] + c1o2*ax[4]*momentU[3]*momentV[2] + 
   c1o2*ay[4]*momentU[2]*momentV[3] + az[0]*momentU[2]*momentW[1] + 
   ax[3]*momentU[3]*momentW[1] + az[1]*momentU[3]*momentW[1] + 
   c1o2*az[4]*momentU[4]*momentW[1] + ay[3]*momentU[2]*momentV[1]*momentW[1] + 
   az[2]*momentU[2]*momentV[1]*momentW[1] + c1o2*az[4]*momentU[2]*momentV[2]*momentW[1] + 
   az[3]*momentU[2]*momentW[2] + c1o2*ax[4]*momentU[3]*momentW[2] + 
   c1o2*ay[4]*momentU[2]*momentV[1]*momentW[2] + c1o2*az[4]*momentU[2]*momentW[3] + 
   c1o2*ax[4]*momentU[3]*momentXi[2] + c1o2*ay[4]*momentU[2]*momentV[1]*momentXi[2] + 
   c1o2*az[4]*momentU[2]*momentW[1]*momentXi[2];

    flux_2.rhoV = ax[0]*momentU[2]*momentV[1] + ax[1]*momentU[3]*momentV[1] + 
   c1o2*ax[4]*momentU[4]*momentV[1] + ay[0]*momentU[1]*momentV[2] + 
   ax[2]*momentU[2]*momentV[2] + ay[1]*momentU[2]*momentV[2] + 
   c1o2*ay[4]*momentU[3]*momentV[2] + ay[2]*momentU[1]*momentV[3] + 
   c1o2*ax[4]*momentU[2]*momentV[3] + c1o2*ay[4]*momentU[1]*momentV[4] + 
   az[0]*momentU[1]*momentV[1]*momentW[1] + ax[3]*momentU[2]*momentV[1]*momentW[1] + 
   az[1]*momentU[2]*momentV[1]*momentW[1] + c1o2*az[4]*momentU[3]*momentV[1]*momentW[1] + 
   ay[3]*momentU[1]*momentV[2]*momentW[1] + az[2]*momentU[1]*momentV[2]*momentW[1] + 
   c1o2*az[4]*momentU[1]*momentV[3]*momentW[1] + az[3]*momentU[1]*momentV[1]*momentW[2] + 
   c1o2*ax[4]*momentU[2]*momentV[1]*momentW[2] + c1o2*ay[4]*momentU[1]*momentV[2]*momentW[2] + 
   c1o2*az[4]*momentU[1]*momentV[1]*momentW[3] + 
   c1o2*ax[4]*momentU[2]*momentV[1]*momentXi[2] + 
   c1o2*ay[4]*momentU[1]*momentV[2]*momentXi[2] + 
   c1o2*az[4]*momentU[1]*momentV[1]*momentW[1]*momentXi[2];
    
    flux_2.rhoW = ax[0]*momentU[2]*momentW[1] + ax[1]*momentU[3]*momentW[1] + 
   c1o2*ax[4]*momentU[4]*momentW[1] + ay[0]*momentU[1]*momentV[1]*momentW[1] + 
   ax[2]*momentU[2]*momentV[1]*momentW[1] + ay[1]*momentU[2]*momentV[1]*momentW[1] + 
   c1o2*ay[4]*momentU[3]*momentV[1]*momentW[1] + ay[2]*momentU[1]*momentV[2]*momentW[1] + 
   c1o2*ax[4]*momentU[2]*momentV[2]*momentW[1] + c1o2*ay[4]*momentU[1]*momentV[3]*momentW[1] + 
   az[0]*momentU[1]*momentW[2] + ax[3]*momentU[2]*momentW[2] + 
   az[1]*momentU[2]*momentW[2] + c1o2*az[4]*momentU[3]*momentW[2] + 
   ay[3]*momentU[1]*momentV[1]*momentW[2] + az[2]*momentU[1]*momentV[1]*momentW[2] + 
   c1o2*az[4]*momentU[1]*momentV[2]*momentW[2] + az[3]*momentU[1]*momentW[3] + 
   c1o2*ax[4]*momentU[2]*momentW[3] + c1o2*ay[4]*momentU[1]*momentV[1]*momentW[3] + 
   c1o2*az[4]*momentU[1]*momentW[4] + c1o2*ax[4]*momentU[2]*momentW[1]*momentXi[2] + 
   c1o2*ay[4]*momentU[1]*momentV[1]*momentW[1]*momentXi[2] + 
   c1o2*az[4]*momentU[1]*momentW[2]*momentXi[2];

    flux_2.rhoE = c1o2*ax[0]*momentU[4] + c1o2*ax[1]*momentU[5] + c1o4*ax[4]*momentU[6] + 
   c1o2*ay[0]*momentU[3]*momentV[1] + c1o2*ax[2]*momentU[4]*momentV[1] + 
   c1o2*ay[1]*momentU[4]*momentV[1] + c1o4*ay[4]*momentU[5]*momentV[1] + 
   c1o2*ax[0]*momentU[2]*momentV[2] + c1o2*ax[1]*momentU[3]*momentV[2] + 
   c1o2*ay[2]*momentU[3]*momentV[2] + c1o2*ax[4]*momentU[4]*momentV[2] + 
   c1o2*ay[0]*momentU[1]*momentV[3] + c1o2*ax[2]*momentU[2]*momentV[3] + 
   c1o2*ay[1]*momentU[2]*momentV[3] + c1o2*ay[4]*momentU[3]*momentV[3] + 
   c1o2*ay[2]*momentU[1]*momentV[4] + c1o4*ax[4]*momentU[2]*momentV[4] + 
   c1o4*ay[4]*momentU[1]*momentV[5] + c1o2*az[0]*momentU[3]*momentW[1] + 
   c1o2*ax[3]*momentU[4]*momentW[1] + c1o2*az[1]*momentU[4]*momentW[1] + 
   c1o4*az[4]*momentU[5]*momentW[1] + c1o2*ay[3]*momentU[3]*momentV[1]*momentW[1] + 
   c1o2*az[2]*momentU[3]*momentV[1]*momentW[1] + 
   c1o2*az[0]*momentU[1]*momentV[2]*momentW[1] + 
   c1o2*ax[3]*momentU[2]*momentV[2]*momentW[1] + c1o2*az[1]*momentU[2]*momentV[2]*momentW[1] + 
   c1o2*az[4]*momentU[3]*momentV[2]*momentW[1] + c1o2*ay[3]*momentU[1]*momentV[3]*momentW[1] + 
   c1o2*az[2]*momentU[1]*momentV[3]*momentW[1] + c1o4*az[4]*momentU[1]*momentV[4]*momentW[1] + 
   c1o2*ax[0]*momentU[2]*momentW[2] + c1o2*ax[1]*momentU[3]*momentW[2] + 
   c1o2*az[3]*momentU[3]*momentW[2] + c1o2*ax[4]*momentU[4]*momentW[2] + 
   c1o2*ay[0]*momentU[1]*momentV[1]*momentW[2] + 
   c1o2*ax[2]*momentU[2]*momentV[1]*momentW[2] + c1o2*ay[1]*momentU[2]*momentV[1]*momentW[2] + 
   c1o2*ay[4]*momentU[3]*momentV[1]*momentW[2] + c1o2*ay[2]*momentU[1]*momentV[2]*momentW[2] + 
   c1o2*az[3]*momentU[1]*momentV[2]*momentW[2] + c1o2*ax[4]*momentU[2]*momentV[2]*momentW[2] + 
   c1o2*ay[4]*momentU[1]*momentV[3]*momentW[2] + c1o2*az[0]*momentU[1]*momentW[3] + 
   c1o2*ax[3]*momentU[2]*momentW[3] + c1o2*az[1]*momentU[2]*momentW[3] + 
   c1o2*az[4]*momentU[3]*momentW[3] + c1o2*ay[3]*momentU[1]*momentV[1]*momentW[3] + 
   c1o2*az[2]*momentU[1]*momentV[1]*momentW[3] + c1o2*az[4]*momentU[1]*momentV[2]*momentW[3] + 
   c1o2*az[3]*momentU[1]*momentW[4] + c1o4*ax[4]*momentU[2]*momentW[4] + 
   c1o4*ay[4]*momentU[1]*momentV[1]*momentW[4] + c1o4*az[4]*momentU[1]*momentW[5] + 
   c1o2*ax[0]*momentU[2]*momentXi[2] + c1o2*ax[1]*momentU[3]*momentXi[2] + 
   c1o2*ax[4]*momentU[4]*momentXi[2] + c1o2*ay[0]*momentU[1]*momentV[1]*momentXi[2] + 
   c1o2*ax[2]*momentU[2]*momentV[1]*momentXi[2] + 
   c1o2*ay[1]*momentU[2]*momentV[1]*momentXi[2] + 
   c1o2*ay[4]*momentU[3]*momentV[1]*momentXi[2] + 
   c1o2*ay[2]*momentU[1]*momentV[2]*momentXi[2] + 
   c1o2*ax[4]*momentU[2]*momentV[2]*momentXi[2] + 
   c1o2*ay[4]*momentU[1]*momentV[3]*momentXi[2] + 
   c1o2*az[0]*momentU[1]*momentW[1]*momentXi[2] + 
   c1o2*ax[3]*momentU[2]*momentW[1]*momentXi[2] + 
   c1o2*az[1]*momentU[2]*momentW[1]*momentXi[2] + 
   c1o2*az[4]*momentU[3]*momentW[1]*momentXi[2] + 
   c1o2*ay[3]*momentU[1]*momentV[1]*momentW[1]*momentXi[2] + 
   c1o2*az[2]*momentU[1]*momentV[1]*momentW[1]*momentXi[2] + 
   c1o2*az[4]*momentU[1]*momentV[2]*momentW[1]*momentXi[2] + 
   c1o2*az[3]*momentU[1]*momentW[2]*momentXi[2] + 
   c1o2*ax[4]*momentU[2]*momentW[2]*momentXi[2] + 
   c1o2*ay[4]*momentU[1]*momentV[1]*momentW[2]*momentXi[2] + 
   c1o2*az[4]*momentU[1]*momentW[3]*momentXi[2] + c1o4*ax[4]*momentU[2]*momentXi[4] + 
   c1o4*ay[4]*momentU[1]*momentV[1]*momentXi[4] + c1o4*az[4]*momentU[1]*momentW[1]*momentXi[4];

    //////////////////////////////////////////////////////////////////////////

    flux_2.rho  += c2o1 * facePrim.lambda * (   momentU[0+1] *                                           facePrim.U - momentU[1+1]                           ) * force.x
                 + c2o1 * facePrim.lambda * (   momentU[0+1] *                                           facePrim.V - momentU[0+1] * momentV[1]              ) * force.y
                 + c2o1 * facePrim.lambda * (   momentU[0+1] *                                           facePrim.W - momentU[0+1] *              momentW[1] ) * force.z ;
                                                                                                         
    flux_2.rhoU += c2o1 * facePrim.lambda * (   momentU[1+1] *                                           facePrim.U - momentU[2+1]                           ) * force.x
                 + c2o1 * facePrim.lambda * (   momentU[1+1] *                                           facePrim.V - momentU[1+1] * momentV[1]              ) * force.y
                 + c2o1 * facePrim.lambda * (   momentU[1+1] *                                           facePrim.W - momentU[1+1] *              momentW[1] ) * force.z ;
                                                                                                         
    flux_2.rhoV += c2o1 * facePrim.lambda * (   momentU[0+1] * momentV[1] *                              facePrim.U - momentU[1+1] * momentV[1]              ) * force.x
                 + c2o1 * facePrim.lambda * (   momentU[0+1] * momentV[1] *                              facePrim.V - momentU[0+1] * momentV[2]              ) * force.y
                 + c2o1 * facePrim.lambda * (   momentU[0+1] * momentV[1] *                              facePrim.W - momentU[0+1] * momentV[1] * momentW[1] ) * force.z ;
                                                                                                         
    flux_2.rhoW += c2o1 * facePrim.lambda * (   momentU[0+1] *              momentW[1] *                 facePrim.U - momentU[1+1] *              momentW[1] ) * force.x
                 + c2o1 * facePrim.lambda * (   momentU[0+1] *              momentW[1] *                 facePrim.V - momentU[0+1] * momentV[1] * momentW[1] ) * force.y
                 + c2o1 * facePrim.lambda * (   momentU[0+1] *              momentW[1] *                 facePrim.W - momentU[0+1] *              momentW[2] ) * force.z ;

    flux_2.rhoE +=       facePrim.lambda * ( ( momentU[2+1] + momentV[2] + momentW[2] + momentXi[2] ) * facePrim.U

                                           - ( momentU[3+1]
                                             + momentU[1+1] * momentV[2] 
                                             + momentU[1+1] *              momentW[2] 
                                             + momentU[1+1] *                           momentXi[2] )
                                           ) * force.x

                 +       facePrim.lambda * ( ( momentU[2+1] + momentV[2] + momentW[2] + momentXi[2] ) * facePrim.V

                                           - ( momentU[2+1] * momentV[1]
                                             + momentU[0+1] * momentV[3]
                                             + momentU[0+1] * momentV[1] * momentW[2]
                                             + momentU[0+1] * momentV[1] *              momentXi[2] )
                                           ) * force.y

                 +       facePrim.lambda * ( ( momentU[2+1] + momentV[2] + momentW[2] + momentXi[2] ) * facePrim.W

                                           - ( momentU[2+1] *              momentW[1]
                                             + momentU[0+1] * momentV[2] * momentW[1]
                                             + momentU[0+1] *              momentW[3]
                                             + momentU[0+1] *              momentW[1] * momentXi[2] )
                                           ) * force.z ;

    //////////////////////////////////////////////////////////////////////////

    flux_3.rho  = at[0]*momentU[1] + at[1]*momentU[2] + c1o2*at[4]*momentU[3] + at[2]*momentU[1]*momentV[1] + 
   c1o2*at[4]*momentU[1]*momentV[2] + at[3]*momentU[1]*momentW[1] + 
   c1o2*at[4]*momentU[1]*momentW[2] + c1o2*at[4]*momentU[1]*momentXi[2];

    flux_3.rhoU = at[0]*momentU[2] + at[1]*momentU[3] + c1o2*at[4]*momentU[4] + at[2]*momentU[2]*momentV[1] + 
   c1o2*at[4]*momentU[2]*momentV[2] + at[3]*momentU[2]*momentW[1] + 
   c1o2*at[4]*momentU[2]*momentW[2] + c1o2*at[4]*momentU[2]*momentXi[2];

    flux_3.rhoV = at[0]*momentU[1]*momentV[1] + at[1]*momentU[2]*momentV[1] + 
   c1o2*at[4]*momentU[3]*momentV[1] + at[2]*momentU[1]*momentV[2] + 
   c1o2*at[4]*momentU[1]*momentV[3] + at[3]*momentU[1]*momentV[1]*momentW[1] + 
   c1o2*at[4]*momentU[1]*momentV[1]*momentW[2] + c1o2*at[4]*momentU[1]*momentV[1]*momentXi[2];

    flux_3.rhoW = at[0]*momentU[1]*momentW[1] + at[1]*momentU[2]*momentW[1] + 
   c1o2*at[4]*momentU[3]*momentW[1] + at[2]*momentU[1]*momentV[1]*momentW[1] + 
   c1o2*at[4]*momentU[1]*momentV[2]*momentW[1] + at[3]*momentU[1]*momentW[2] + 
   c1o2*at[4]*momentU[1]*momentW[3] + c1o2*at[4]*momentU[1]*momentW[1]*momentXi[2];

    flux_3.rhoE = c1o2*at[0]*momentU[3] + c1o2*at[1]*momentU[4] + c1o4*at[4]*momentU[5] + 
   c1o2*at[2]*momentU[3]*momentV[1] + c1o2*at[0]*momentU[1]*momentV[2] + 
   c1o2*at[1]*momentU[2]*momentV[2] + c1o2*at[4]*momentU[3]*momentV[2] + 
   c1o2*at[2]*momentU[1]*momentV[3] + c1o4*at[4]*momentU[1]*momentV[4] + 
   c1o2*at[3]*momentU[3]*momentW[1] + c1o2*at[3]*momentU[1]*momentV[2]*momentW[1] + 
   c1o2*at[0]*momentU[1]*momentW[2] + c1o2*at[1]*momentU[2]*momentW[2] + 
   c1o2*at[4]*momentU[3]*momentW[2] + c1o2*at[2]*momentU[1]*momentV[1]*momentW[2] + 
   c1o2*at[4]*momentU[1]*momentV[2]*momentW[2] + c1o2*at[3]*momentU[1]*momentW[3] + 
   c1o4*at[4]*momentU[1]*momentW[4] + c1o2*at[0]*momentU[1]*momentXi[2] + 
   c1o2*at[1]*momentU[2]*momentXi[2] + c1o2*at[4]*momentU[3]*momentXi[2] + 
   c1o2*at[2]*momentU[1]*momentV[1]*momentXi[2] + 
   c1o2*at[4]*momentU[1]*momentV[2]*momentXi[2] + 
   c1o2*at[3]*momentU[1]*momentW[1]*momentXi[2] + 
   c1o2*at[4]*momentU[1]*momentW[2]*momentXi[2] + c1o4*at[4]*momentU[1]*momentXi[4];

    //////////////////////////////////////////////////////////////////////////

    flux.rho  = ( timeCoefficients[0] * flux_1.rho  + timeCoefficients[1] * flux_2.rho  + timeCoefficients[2] * flux_3.rho  ) * parameters.dx * parameters.dx * facePrim.rho;
    flux.rhoU = ( timeCoefficients[0] * flux_1.rhoU + timeCoefficients[1] * flux_2.rhoU + timeCoefficients[2] * flux_3.rhoU ) * parameters.dx * parameters.dx * facePrim.rho;
    flux.rhoV = ( timeCoefficients[0] * flux_1.rhoV + timeCoefficients[1] * flux_2.rhoV + timeCoefficients[2] * flux_3.rhoV ) * parameters.dx * parameters.dx * facePrim.rho;
    flux.rhoW = ( timeCoefficients[0] * flux_1.rhoW + timeCoefficients[1] * flux_2.rhoW + timeCoefficients[2] * flux_3.rhoW ) * parameters.dx * parameters.dx * facePrim.rho;
    flux.rhoE = ( timeCoefficients[0] * flux_1.rhoE + timeCoefficients[1] * flux_2.rhoE + timeCoefficients[2] * flux_3.rhoE ) * parameters.dx * parameters.dx * facePrim.rho;

    //////////////////////////////////////////////////////////////////////////

    real qConstPr = timeCoefficients[1] * ( ( flux_2.rhoE + flux_3.rhoE ) 
                                          - facePrim.U * ( flux_2.rhoU + flux_3.rhoU ) 
                                          - facePrim.V * ( flux_2.rhoV + flux_3.rhoV ) 
                                          - facePrim.W * ( flux_2.rhoW + flux_3.rhoW )
                                          ) * parameters.dx * parameters.dx * facePrim.rho;

    flux.rhoE += ( c1o1 / parameters.Pr - c1o1 ) * qConstPr;

    heatFlux   = qConstPr / parameters.Pr;

    //////////////////////////////////////////////////////////////////////////

#ifdef USE_PASSIVE_SCALAR

	flux_2.rhoS_1 = ( ax[5] * momentU[1+1]                          
				    + ay[5] * momentU[0+1] * momentV[1]             
				    + az[5] * momentU[0+1] *              momentW[1]
				    ) / (c2o1 * facePrim.lambda);

	flux_2.rhoS_2 = ( ax[6] * momentU[1+1]                          
				    + ay[6] * momentU[0+1] * momentV[1]             
				    + az[6] * momentU[0+1] *              momentW[1]
				    ) / (c2o1 * facePrim.lambda);

	flux_3.rhoS_1 = at[5] * momentU[0 + 1]
				  / ( c2o1 * facePrim.lambda );

	flux_3.rhoS_2 = at[6] * momentU[0 + 1]
				  / ( c2o1 * facePrim.lambda );

    //////////////////////////////////////////////////////////////////////////

	real dt = parameters.dt;
    real timeCoefficientsScalar [3];

    {
        real tauS = parameters.D1 * c2o1 * facePrim.lambda;
        timeCoefficientsScalar[0] = dt;
        timeCoefficientsScalar[1] = -tauS * dt;
        timeCoefficientsScalar[2] = c1o2 * dt * dt - tauS * dt;

        flux.rhoS_1 = flux.rho * facePrim.S_1 + ( timeCoefficientsScalar[1] * flux_2.rhoS_1 + timeCoefficientsScalar[2] * flux_3.rhoS_1 ) * parameters.dx * parameters.dx * facePrim.rho;
    }

    {
        real tauS = parameters.D2 * c2o1 * facePrim.lambda;
        timeCoefficientsScalar[0] = dt;
        timeCoefficientsScalar[1] = -tauS * dt;
        timeCoefficientsScalar[2] = c1o2 * dt * dt - tauS * dt;

        flux.rhoS_2 = flux.rho * facePrim.S_2 + ( timeCoefficientsScalar[1] * flux_2.rhoS_2 + timeCoefficientsScalar[2] * flux_3.rhoS_2 ) * parameters.dx * parameters.dx * facePrim.rho;
    }

#endif // USE_PASSIVE_SCALAR
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif