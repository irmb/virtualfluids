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
//! \author Martin Schoenherr
//=======================================================================================
#include "Calculation/Calculation.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void InitNavierStokesIncompressible_Device(unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    unsigned int* geoD,
    real* rho,
    real* ux,
    real* uy,
    real* uz,
    unsigned int size_Mat,
    real* DD,
    bool EvenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<size_Mat)
   {
      ////////////////////////////////////////////////////////////////////////////////
      unsigned int BC;
      BC        =   geoD[k];

      if( BC != GEO_SOLID &&  BC != GEO_VOID)
      {
         Distributions27 D;
         if (EvenOrOdd==true)
         {
            D.f[dP00   ] = &DD[dP00   *size_Mat];
            D.f[dM00   ] = &DD[dM00   *size_Mat];
            D.f[d0P0   ] = &DD[d0P0   *size_Mat];
            D.f[d0M0   ] = &DD[d0M0   *size_Mat];
            D.f[d00P   ] = &DD[d00P   *size_Mat];
            D.f[d00M   ] = &DD[d00M   *size_Mat];
            D.f[dPP0  ] = &DD[dPP0  *size_Mat];
            D.f[dMM0  ] = &DD[dMM0  *size_Mat];
            D.f[dPM0  ] = &DD[dPM0  *size_Mat];
            D.f[dMP0  ] = &DD[dMP0  *size_Mat];
            D.f[dP0P  ] = &DD[dP0P  *size_Mat];
            D.f[dM0M  ] = &DD[dM0M  *size_Mat];
            D.f[dP0M  ] = &DD[dP0M  *size_Mat];
            D.f[dM0P  ] = &DD[dM0P  *size_Mat];
            D.f[d0PP  ] = &DD[d0PP  *size_Mat];
            D.f[d0MM  ] = &DD[d0MM  *size_Mat];
            D.f[d0PM  ] = &DD[d0PM  *size_Mat];
            D.f[d0MP  ] = &DD[d0MP  *size_Mat];
            D.f[d000] = &DD[d000*size_Mat];
            D.f[dPPP ] = &DD[dPPP *size_Mat];
            D.f[dMMP ] = &DD[dMMP *size_Mat];
            D.f[dPMP ] = &DD[dPMP *size_Mat];
            D.f[dMPP ] = &DD[dMPP *size_Mat];
            D.f[dPPM ] = &DD[dPPM *size_Mat];
            D.f[dMMM ] = &DD[dMMM *size_Mat];
            D.f[dPMM ] = &DD[dPMM *size_Mat];
            D.f[dMPM ] = &DD[dMPM *size_Mat];
         }
         else
         {
            D.f[dM00   ] = &DD[dP00   *size_Mat];
            D.f[dP00   ] = &DD[dM00   *size_Mat];
            D.f[d0M0   ] = &DD[d0P0   *size_Mat];
            D.f[d0P0   ] = &DD[d0M0   *size_Mat];
            D.f[d00M   ] = &DD[d00P   *size_Mat];
            D.f[d00P   ] = &DD[d00M   *size_Mat];
            D.f[dMM0  ] = &DD[dPP0  *size_Mat];
            D.f[dPP0  ] = &DD[dMM0  *size_Mat];
            D.f[dMP0  ] = &DD[dPM0  *size_Mat];
            D.f[dPM0  ] = &DD[dMP0  *size_Mat];
            D.f[dM0M  ] = &DD[dP0P  *size_Mat];
            D.f[dP0P  ] = &DD[dM0M  *size_Mat];
            D.f[dM0P  ] = &DD[dP0M  *size_Mat];
            D.f[dP0M  ] = &DD[dM0P  *size_Mat];
            D.f[d0MM  ] = &DD[d0PP  *size_Mat];
            D.f[d0PP  ] = &DD[d0MM  *size_Mat];
            D.f[d0MP  ] = &DD[d0PM  *size_Mat];
            D.f[d0PM  ] = &DD[d0MP  *size_Mat];
            D.f[d000] = &DD[d000*size_Mat];
            D.f[dMMM ] = &DD[dPPP *size_Mat];
            D.f[dPPM ] = &DD[dMMP *size_Mat];
            D.f[dMPM ] = &DD[dPMP *size_Mat];
            D.f[dPMM ] = &DD[dMPP *size_Mat];
            D.f[dMMP ] = &DD[dPPM *size_Mat];
            D.f[dPPP ] = &DD[dMMM *size_Mat];
            D.f[dMPP ] = &DD[dPMM *size_Mat];
            D.f[dPMP ] = &DD[dMPM *size_Mat];
         }
         //////////////////////////////////////////////////////////////////////////
         real drho = rho[k];//0.0f;//
         real  vx1 = ux[k]; //0.0f;//
         real  vx2 = uy[k]; //0.0f;//
         real  vx3 = uz[k]; //0.0f;//
         //////////////////////////////////////////////////////////////////////////
         //index
         //////////////////////////////////////////////////////////////////////////
         unsigned int kzero= k;
         unsigned int ke   = k;
         unsigned int kw   = neighborX[k];
         unsigned int kn   = k;
         unsigned int ks   = neighborY[k];
         unsigned int kt   = k;
         unsigned int kb   = neighborZ[k];
         unsigned int ksw  = neighborY[kw];
         unsigned int kne  = k;
         unsigned int kse  = ks;
         unsigned int knw  = kw;
         unsigned int kbw  = neighborZ[kw];
         unsigned int kte  = k;
         unsigned int kbe  = kb;
         unsigned int ktw  = kw;
         unsigned int kbs  = neighborZ[ks];
         unsigned int ktn  = k;
         unsigned int kbn  = kb;
         unsigned int kts  = ks;
         unsigned int ktse = ks;
         unsigned int kbnw = kbw;
         unsigned int ktnw = kw;
         unsigned int kbse = kbs;
         unsigned int ktsw = ksw;
         unsigned int kbne = kb;
         unsigned int ktne = k;
         unsigned int kbsw = neighborZ[ksw];
         //////////////////////////////////////////////////////////////////////////
         real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D.f[d000])[kzero] =   c8o27* (drho-cu_sq);
         (D.f[dP00   ])[ke   ] =   c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
         (D.f[dM00   ])[kw   ] =   c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
         (D.f[d0P0   ])[kn   ] =   c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
         (D.f[d0M0   ])[ks   ] =   c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D.f[d00P   ])[kt   ] =   c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
         (D.f[d00M   ])[kb   ] =   c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
         (D.f[dPP0  ])[kne  ] =   c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D.f[dMM0  ])[ksw  ] =   c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D.f[dPM0  ])[kse  ] =   c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D.f[dMP0  ])[knw  ] =   c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D.f[dP0P  ])[kte  ] =   c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D.f[dM0M  ])[kbw  ] =   c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D.f[dP0M  ])[kbe  ] =   c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D.f[dM0P  ])[ktw  ] =   c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D.f[d0PP  ])[ktn  ] =   c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D.f[d0MM  ])[kbs  ] =   c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D.f[d0PM  ])[kbn  ] =   c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D.f[d0MP  ])[kts  ] =   c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D.f[dPPP ])[ktne ] =   c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D.f[dMMM ])[kbsw ] =   c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D.f[dPPM ])[kbne ] =   c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D.f[dMMP ])[ktsw ] =   c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D.f[dPMP ])[ktse ] =   c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D.f[dMPM ])[kbnw ] =   c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D.f[dPMM ])[kbse ] =   c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D.f[dMPP ])[ktnw ] =   c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      }
      else
      {
          //////////////////////////////////////////////////////////////////////////
          Distributions27 D;
          D.f[d000] = &DD[d000*size_Mat];
          //////////////////////////////////////////////////////////////////////////
          (D.f[d000])[k] = c96o1;
          //////////////////////////////////////////////////////////////////////////
      }
   }
}