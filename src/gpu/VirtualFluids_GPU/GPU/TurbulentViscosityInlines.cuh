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
//! \file TurbulentViscosity.h
//! \ingroup GPU
//! \author Henry Korb, Henrik Asmuth
//======================================================================================

#ifndef TURBULENT_VISCOSITY_INLINES_CUH_
#define TURBULENT_VISCOSITY_INLINES_CUH_

#include <cuda.h>
#include <cuda_runtime.h>

#include "LBM/LB.h" 
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;

__inline__ __device__ real calcTurbulentViscositySmagorinsky(real Cs, real dxux, real dyuy, real dzuz, real Dxy, real Dxz , real Dyz)
{
    return Cs*Cs * sqrt( c2o1 * ( dxux*dxux + dyuy*dyuy + dzuz*dzuz ) + Dxy*Dxy + Dxz*Dxz + Dyz*Dyz );
}

__inline__ __device__ real calcTurbulentViscosityQR(real C, real dxux, real dyuy, real dzuz, real Dxy, real Dxz , real Dyz)
{
        // ! Verstappen's QR model
        //! Second invariant of the strain-rate tensor
        real Q = c1o2*( dxux*dxux + dyuy*dyuy + dzuz*dzuz ) + c1o4*( Dxy*Dxy + Dxz*Dxz + Dyz*Dyz);
        //! Third invariant of the strain-rate tensor (determinant)
        real R = - dxux*dyuy*dzuz - c1o4*( Dxy*Dxz*Dyz + dxux*Dyz*Dyz + dyuy*Dxz*Dxz + dzuz*Dxy*Dxy );
        
        return C * max(R, c0o1) / Q;
}

#endif //TURBULENT_VISCOSITY_H_e