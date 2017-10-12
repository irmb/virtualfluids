/* Device code */
#include "LBM/D3Q27.h"
//#include "LBM/LB.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////
__global__ void LBInitSP27( unsigned int* neighborX,
                                       unsigned int* neighborY,
                                       unsigned int* neighborZ,
                                       unsigned int* geoD,
                                       doubflo* rho,
                                       doubflo* ux,
                                       doubflo* uy,
                                       doubflo* uz,
                                       unsigned int size_Mat,
                                       doubflo* DD,
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
            D.f[dirE   ] = &DD[dirE   *size_Mat];
            D.f[dirW   ] = &DD[dirW   *size_Mat];
            D.f[dirN   ] = &DD[dirN   *size_Mat];
            D.f[dirS   ] = &DD[dirS   *size_Mat];
            D.f[dirT   ] = &DD[dirT   *size_Mat];
            D.f[dirB   ] = &DD[dirB   *size_Mat];
            D.f[dirNE  ] = &DD[dirNE  *size_Mat];
            D.f[dirSW  ] = &DD[dirSW  *size_Mat];
            D.f[dirSE  ] = &DD[dirSE  *size_Mat];
            D.f[dirNW  ] = &DD[dirNW  *size_Mat];
            D.f[dirTE  ] = &DD[dirTE  *size_Mat];
            D.f[dirBW  ] = &DD[dirBW  *size_Mat];
            D.f[dirBE  ] = &DD[dirBE  *size_Mat];
            D.f[dirTW  ] = &DD[dirTW  *size_Mat];
            D.f[dirTN  ] = &DD[dirTN  *size_Mat];
            D.f[dirBS  ] = &DD[dirBS  *size_Mat];
            D.f[dirBN  ] = &DD[dirBN  *size_Mat];
            D.f[dirTS  ] = &DD[dirTS  *size_Mat];
            D.f[dirZERO] = &DD[dirZERO*size_Mat];
            D.f[dirTNE ] = &DD[dirTNE *size_Mat];
            D.f[dirTSW ] = &DD[dirTSW *size_Mat];
            D.f[dirTSE ] = &DD[dirTSE *size_Mat];
            D.f[dirTNW ] = &DD[dirTNW *size_Mat];
            D.f[dirBNE ] = &DD[dirBNE *size_Mat];
            D.f[dirBSW ] = &DD[dirBSW *size_Mat];
            D.f[dirBSE ] = &DD[dirBSE *size_Mat];
            D.f[dirBNW ] = &DD[dirBNW *size_Mat];
         }
         else
         {
            D.f[dirW   ] = &DD[dirE   *size_Mat];
            D.f[dirE   ] = &DD[dirW   *size_Mat];
            D.f[dirS   ] = &DD[dirN   *size_Mat];
            D.f[dirN   ] = &DD[dirS   *size_Mat];
            D.f[dirB   ] = &DD[dirT   *size_Mat];
            D.f[dirT   ] = &DD[dirB   *size_Mat];
            D.f[dirSW  ] = &DD[dirNE  *size_Mat];
            D.f[dirNE  ] = &DD[dirSW  *size_Mat];
            D.f[dirNW  ] = &DD[dirSE  *size_Mat];
            D.f[dirSE  ] = &DD[dirNW  *size_Mat];
            D.f[dirBW  ] = &DD[dirTE  *size_Mat];
            D.f[dirTE  ] = &DD[dirBW  *size_Mat];
            D.f[dirTW  ] = &DD[dirBE  *size_Mat];
            D.f[dirBE  ] = &DD[dirTW  *size_Mat];
            D.f[dirBS  ] = &DD[dirTN  *size_Mat];
            D.f[dirTN  ] = &DD[dirBS  *size_Mat];
            D.f[dirTS  ] = &DD[dirBN  *size_Mat];
            D.f[dirBN  ] = &DD[dirTS  *size_Mat];
            D.f[dirZERO] = &DD[dirZERO*size_Mat];
            D.f[dirBSW ] = &DD[dirTNE *size_Mat];
            D.f[dirBNE ] = &DD[dirTSW *size_Mat];
            D.f[dirBNW ] = &DD[dirTSE *size_Mat];
            D.f[dirBSE ] = &DD[dirTNW *size_Mat];
            D.f[dirTSW ] = &DD[dirBNE *size_Mat];
            D.f[dirTNE ] = &DD[dirBSW *size_Mat];
            D.f[dirTNW ] = &DD[dirBSE *size_Mat];
            D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         }
         //////////////////////////////////////////////////////////////////////////
         doubflo drho = rho[k];//0.0f;//
         doubflo  vx1 = ux[k]; //0.0f;//
         doubflo  vx2 = uy[k]; //0.0f;//
         doubflo  vx3 = uz[k]; //0.0f;//
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
         doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

         (D.f[dirZERO])[kzero] =   c8over27* (drho-cu_sq);
         (D.f[dirE   ])[ke   ] =   c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
         (D.f[dirW   ])[kw   ] =   c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
         (D.f[dirN   ])[kn   ] =   c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
         (D.f[dirS   ])[ks   ] =   c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
         (D.f[dirT   ])[kt   ] =   c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
         (D.f[dirB   ])[kb   ] =   c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
         (D.f[dirNE  ])[kne  ] =   c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq);
         (D.f[dirSW  ])[ksw  ] =   c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq);
         (D.f[dirSE  ])[kse  ] =   c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq);
         (D.f[dirNW  ])[knw  ] =   c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq);
         (D.f[dirTE  ])[kte  ] =   c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq);
         (D.f[dirBW  ])[kbw  ] =   c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq);
         (D.f[dirBE  ])[kbe  ] =   c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq);
         (D.f[dirTW  ])[ktw  ] =   c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq);
         (D.f[dirTN  ])[ktn  ] =   c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq);
         (D.f[dirBS  ])[kbs  ] =   c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq);
         (D.f[dirBN  ])[kbn  ] =   c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq);
         (D.f[dirTS  ])[kts  ] =   c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq);
         (D.f[dirTNE ])[ktne ] =   c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq);
         (D.f[dirBSW ])[kbsw ] =   c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq);
         (D.f[dirBNE ])[kbne ] =   c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq);
         (D.f[dirTSW ])[ktsw ] =   c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq);
         (D.f[dirTSE ])[ktse ] =   c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq);
         (D.f[dirBNW ])[kbnw ] =   c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq);
         (D.f[dirBSE ])[kbse ] =   c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq);
         (D.f[dirTNW ])[ktnw ] =   c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq);
      }
   }
}





