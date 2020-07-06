#include "LBM/LB.h" 
#include "LBM/D3Q27.h"
#include "Core/RealConstants.h"
#include "math.h"

extern "C" __global__ void LB_Init_Incomp_AD_7(unsigned int* neighborX,
	unsigned int* neighborY,
	unsigned int* neighborZ,
	unsigned int* geoD,
	real* Conc,
	real* ux,
	real* uy,
	real* uz,
	unsigned int size_Mat,
	real* DD7,
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

      if( BC != GEO_SOLID && BC != GEO_VOID)
      {
         Distributions7 D7;
         if (EvenOrOdd==true)
         {
            D7.f[0] = &DD7[0*size_Mat];
            D7.f[1] = &DD7[1*size_Mat];
            D7.f[2] = &DD7[2*size_Mat];
            D7.f[3] = &DD7[3*size_Mat];
            D7.f[4] = &DD7[4*size_Mat];
            D7.f[5] = &DD7[5*size_Mat];
            D7.f[6] = &DD7[6*size_Mat];
         }
         else
         {
            D7.f[0] = &DD7[0*size_Mat];
            D7.f[2] = &DD7[1*size_Mat];
            D7.f[1] = &DD7[2*size_Mat];
            D7.f[4] = &DD7[3*size_Mat];
            D7.f[3] = &DD7[4*size_Mat];
            D7.f[6] = &DD7[5*size_Mat];
            D7.f[5] = &DD7[6*size_Mat];
         }
         //////////////////////////////////////////////////////////////////////////
         real ConcD = Conc[k];
         real   vx1 = ux[k];
         real   vx2 = uy[k];
         real   vx3 = uz[k];
         real lambdaD     = -c3o1 + sqrt(c3o1);
         real Diffusivity = c1o20;
         real Lam         = -(c1o2+c1o1/lambdaD);
         real nue_d       = Lam/c3o1;
         real ae          = Diffusivity/nue_d - c1o1;
         real ux_sq       = vx1 * vx1;
         real uy_sq       = vx2 * vx2;
         real uz_sq       = vx3 * vx3;
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
         //////////////////////////////////////////////////////////////////////////

         (D7.f[0])[kzero] = ConcD*(c1o3*(ae*(-c3o1))-(ux_sq+uy_sq+uz_sq));
         (D7.f[1])[ke   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)+vx1*c1o2);
         (D7.f[2])[kw   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(ux_sq)-vx1*c1o2);
         (D7.f[3])[kn   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)+vx2*c1o2);
         (D7.f[4])[ks   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uy_sq)-vx2*c1o2);
         (D7.f[5])[kt   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)+vx3*c1o2);
         (D7.f[6])[kb   ] = ConcD*(c1o6*(ae+c1o1)+c1o2*(uz_sq)-vx3*c1o2);
      }
   }
}