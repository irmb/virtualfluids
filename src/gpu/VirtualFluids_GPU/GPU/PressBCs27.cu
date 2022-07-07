/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include "lbm/constants/NumericConstants.h"
#include "KernelUtilities.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QInflowScaleByPressDevice27(  real* rhoBC,
														 real* DD, 
														 int* k_Q, 
														 int* k_N, 
														 int numberOfBCnodes, 
														 real om1, 
														 unsigned int* neighborX,
														 unsigned int* neighborY,
														 unsigned int* neighborZ,
														 unsigned int size_Mat, 
														 bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int K1QK  = k_N[k];
      unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==true)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f1_E    = (D.f[E   ])[k1e   ];
      real f1_W    = (D.f[W   ])[k1w   ];
      real f1_N    = (D.f[N   ])[k1n   ];
      real f1_S    = (D.f[S   ])[k1s   ];
      real f1_T    = (D.f[T   ])[k1t   ];
      real f1_B    = (D.f[B   ])[k1b   ];
      real f1_NE   = (D.f[NE  ])[k1ne  ];
      real f1_SW   = (D.f[SW  ])[k1sw  ];
      real f1_SE   = (D.f[SE  ])[k1se  ];
      real f1_NW   = (D.f[NW  ])[k1nw  ];
      real f1_TE   = (D.f[TE  ])[k1te  ];
      real f1_BW   = (D.f[BW  ])[k1bw  ];
      real f1_BE   = (D.f[BE  ])[k1be  ];
      real f1_TW   = (D.f[TW  ])[k1tw  ];
      real f1_TN   = (D.f[TN  ])[k1tn  ];
      real f1_BS   = (D.f[BS  ])[k1bs  ];
      real f1_BN   = (D.f[BN  ])[k1bn  ];
      real f1_TS   = (D.f[TS  ])[k1ts  ];
      //real f1_ZERO = (D.f[REST])[k1zero];
      real f1_TNE  = (D.f[TNE ])[k1tne ];
      real f1_TSW  = (D.f[TSW ])[k1tsw ];
      real f1_TSE  = (D.f[TSE ])[k1tse ];
      real f1_TNW  = (D.f[TNW ])[k1tnw ];
      real f1_BNE  = (D.f[BNE ])[k1bne ];
      real f1_BSW  = (D.f[BSW ])[k1bsw ];
      real f1_BSE  = (D.f[BSE ])[k1bse ];
      real f1_BNW  = (D.f[BNW ])[k1bnw ];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E    = (D.f[E   ])[ke   ];
      real f_W    = (D.f[W   ])[kw   ];
      real f_N    = (D.f[N   ])[kn   ];
      real f_S    = (D.f[S   ])[ks   ];
      real f_T    = (D.f[T   ])[kt   ];
      real f_B    = (D.f[B   ])[kb   ];
      real f_NE   = (D.f[NE  ])[kne  ];
      real f_SW   = (D.f[SW  ])[ksw  ];
      real f_SE   = (D.f[SE  ])[kse  ];
      real f_NW   = (D.f[NW  ])[knw  ];
      real f_TE   = (D.f[TE  ])[kte  ];
      real f_BW   = (D.f[BW  ])[kbw  ];
      real f_BE   = (D.f[BE  ])[kbe  ];
      real f_TW   = (D.f[TW  ])[ktw  ];
      real f_TN   = (D.f[TN  ])[ktn  ];
      real f_BS   = (D.f[BS  ])[kbs  ];
      real f_BN   = (D.f[BN  ])[kbn  ];
      real f_TS   = (D.f[TS  ])[kts  ];
      //real f_ZERO = (D.f[REST])[kzero];
      real f_TNE  = (D.f[TNE ])[ktne ];
      real f_TSW  = (D.f[TSW ])[ktsw ];
      real f_TSE  = (D.f[TSE ])[ktse ];
      real f_TNW  = (D.f[TNW ])[ktnw ];
      real f_BNE  = (D.f[BNE ])[kbne ];
      real f_BSW  = (D.f[BSW ])[kbsw ];
      real f_BSE  = (D.f[BSE ])[kbse ];
      real f_BNW  = (D.f[BNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      // real vx1, vx2, vx3;
      real drho, drho1;
      //////////////////////////////////////////////////////////////////////////
	  //Dichte
      drho1  =  f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW + 
                f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((D.f[REST])[k1zero]); 
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 
      //////////////////////////////////////////////////////////////////////////
	  //Schallgeschwindigkeit
	  real cs = c1o1 / sqrtf(c3o1);
      //////////////////////////////////////////////////////////////////////////
	  real rhoInterpol = drho1 * cs + (c1o1 - cs) * drho; 
	  //real diffRho = (rhoBC[k] + one) / (rhoInterpol + one);
	  real diffRhoToAdd = rhoBC[k] - rhoInterpol;
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //no velocity
	  //////////////////////////////////////////
      f_E    = f1_E   * cs + (c1o1 - cs) * f_E   ;
      f_W    = f1_W   * cs + (c1o1 - cs) * f_W   ;
      f_N    = f1_N   * cs + (c1o1 - cs) * f_N   ;
      f_S    = f1_S   * cs + (c1o1 - cs) * f_S   ;
      f_T    = f1_T   * cs + (c1o1 - cs) * f_T   ;
      f_B    = f1_B   * cs + (c1o1 - cs) * f_B   ;
      f_NE   = f1_NE  * cs + (c1o1 - cs) * f_NE  ;
      f_SW   = f1_SW  * cs + (c1o1 - cs) * f_SW  ;
      f_SE   = f1_SE  * cs + (c1o1 - cs) * f_SE  ;
      f_NW   = f1_NW  * cs + (c1o1 - cs) * f_NW  ;
      f_TE   = f1_TE  * cs + (c1o1 - cs) * f_TE  ;
      f_BW   = f1_BW  * cs + (c1o1 - cs) * f_BW  ;
      f_BE   = f1_BE  * cs + (c1o1 - cs) * f_BE  ;
      f_TW   = f1_TW  * cs + (c1o1 - cs) * f_TW  ;
      f_TN   = f1_TN  * cs + (c1o1 - cs) * f_TN  ;
      f_BS   = f1_BS  * cs + (c1o1 - cs) * f_BS  ;
      f_BN   = f1_BN  * cs + (c1o1 - cs) * f_BN  ;
      f_TS   = f1_TS  * cs + (c1o1 - cs) * f_TS  ;
      f_TNE  = f1_TNE * cs + (c1o1 - cs) * f_TNE ;
      f_TSW  = f1_TSW * cs + (c1o1 - cs) * f_TSW ;
      f_TSE  = f1_TSE * cs + (c1o1 - cs) * f_TSE ;
      f_TNW  = f1_TNW * cs + (c1o1 - cs) * f_TNW ;
      f_BNE  = f1_BNE * cs + (c1o1 - cs) * f_BNE ;
      f_BSW  = f1_BSW * cs + (c1o1 - cs) * f_BSW ;
      f_BSE  = f1_BSE * cs + (c1o1 - cs) * f_BSE ;
      f_BNW  = f1_BNW * cs + (c1o1 - cs) * f_BNW ;
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //scale by press
	  //////////////////////////////////////////
	  //f_E    = (f_E   + c2over27 ) * diffRho - c2over27 ;
   //   f_W    = (f_W   + c2over27 ) * diffRho - c2over27 ;
   //   f_N    = (f_N   + c2over27 ) * diffRho - c2over27 ;
   //   f_S    = (f_S   + c2over27 ) * diffRho - c2over27 ;
   //   f_T    = (f_T   + c2over27 ) * diffRho - c2over27 ;
   //   f_B    = (f_B   + c2over27 ) * diffRho - c2over27 ;
	  //f_NE   = (f_NE  + c1over54 ) * diffRho - c1over54 ;
   //   f_SW   = (f_SW  + c1over54 ) * diffRho - c1over54 ;
   //   f_SE   = (f_SE  + c1over54 ) * diffRho - c1over54 ;
   //   f_NW   = (f_NW  + c1over54 ) * diffRho - c1over54 ;
   //   f_TE   = (f_TE  + c1over54 ) * diffRho - c1over54 ;
   //   f_BW   = (f_BW  + c1over54 ) * diffRho - c1over54 ;
   //   f_BE   = (f_BE  + c1over54 ) * diffRho - c1over54 ;
   //   f_TW   = (f_TW  + c1over54 ) * diffRho - c1over54 ;
   //   f_TN   = (f_TN  + c1over54 ) * diffRho - c1over54 ;
   //   f_BS   = (f_BS  + c1over54 ) * diffRho - c1over54 ;
   //   f_BN   = (f_BN  + c1over54 ) * diffRho - c1over54 ;
   //   f_TS   = (f_TS  + c1over54 ) * diffRho - c1over54 ;
   //   f_TNE  = (f_TNE + c1over216) * diffRho - c1over216;
   //   f_TSW  = (f_TSW + c1over216) * diffRho - c1over216;
   //   f_TSE  = (f_TSE + c1over216) * diffRho - c1over216;
   //   f_TNW  = (f_TNW + c1over216) * diffRho - c1over216;
   //   f_BNE  = (f_BNE + c1over216) * diffRho - c1over216;
   //   f_BSW  = (f_BSW + c1over216) * diffRho - c1over216;
   //   f_BSE  = (f_BSE + c1over216) * diffRho - c1over216;
   //   f_BNW  = (f_BNW + c1over216) * diffRho - c1over216;
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  // add press
	  //////////////////////////////////////////
	  f_E    = (f_E   + c2o27  * diffRhoToAdd);
      f_W    = (f_W   + c2o27  * diffRhoToAdd);
      f_N    = (f_N   + c2o27  * diffRhoToAdd);
      f_S    = (f_S   + c2o27  * diffRhoToAdd);
      f_T    = (f_T   + c2o27  * diffRhoToAdd);
      f_B    = (f_B   + c2o27  * diffRhoToAdd);
	  f_NE   = (f_NE  + c1o54  * diffRhoToAdd);
      f_SW   = (f_SW  + c1o54  * diffRhoToAdd);
      f_SE   = (f_SE  + c1o54  * diffRhoToAdd);
      f_NW   = (f_NW  + c1o54  * diffRhoToAdd);
      f_TE   = (f_TE  + c1o54  * diffRhoToAdd);
      f_BW   = (f_BW  + c1o54  * diffRhoToAdd);
      f_BE   = (f_BE  + c1o54  * diffRhoToAdd);
      f_TW   = (f_TW  + c1o54  * diffRhoToAdd);
      f_TN   = (f_TN  + c1o54  * diffRhoToAdd);
      f_BS   = (f_BS  + c1o54  * diffRhoToAdd);
      f_BN   = (f_BN  + c1o54  * diffRhoToAdd);
      f_TS   = (f_TS  + c1o54  * diffRhoToAdd);
      f_TNE  = (f_TNE + c1o216 * diffRhoToAdd);
      f_TSW  = (f_TSW + c1o216 * diffRhoToAdd);
      f_TSE  = (f_TSE + c1o216 * diffRhoToAdd);
      f_TNW  = (f_TNW + c1o216 * diffRhoToAdd);
      f_BNE  = (f_BNE + c1o216 * diffRhoToAdd);
      f_BSW  = (f_BSW + c1o216 * diffRhoToAdd);
      f_BSE  = (f_BSE + c1o216 * diffRhoToAdd);
      f_BNW  = (f_BNW + c1o216 * diffRhoToAdd);
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////
      //__syncthreads();
	  // -X
	  //(D.f[E   ])[ke   ] = f_E   ;
	  //(D.f[SE  ])[kse  ] = f_SE  ;
	  //(D.f[NE  ])[kne  ] = f_NE  ;
	  //(D.f[BE  ])[kbe  ] = f_BE  ;
	  //(D.f[TE  ])[kte  ] = f_TE  ;
	  //(D.f[TSE ])[ktse ] = f_TSE ;
	  //(D.f[TNE ])[ktne ] = f_TNE ;
	  //(D.f[BSE ])[kbse ] = f_BSE ;
	  //(D.f[BNE ])[kbne ] = f_BNE ;     
	  // X
	  (D.f[W   ])[kw   ] = f_W   ;
	  (D.f[SW  ])[ksw  ] = f_SW  ;
	  (D.f[NW  ])[knw  ] = f_NW  ;
	  (D.f[BW  ])[kbw  ] = f_BW  ;
	  (D.f[TW  ])[ktw  ] = f_TW  ;
	  (D.f[TSW ])[ktsw ] = f_TSW ;
	  (D.f[TNW ])[ktnw ] = f_TNW ;
	  (D.f[BSW ])[kbsw ] = f_BSW ;
	  (D.f[BNW ])[kbnw ] = f_BNW ;     
	  // Y
	  //(D.f[S   ])[ks   ] = f_S   ;
	  //(D.f[SE  ])[kse  ] = f_SE  ;
	  //(D.f[SW  ])[ksw  ] = f_SW  ;
	  //(D.f[TS  ])[kts  ] = f_TS  ;
	  //(D.f[BS  ])[kbs  ] = f_BS  ;
	  //(D.f[TSE ])[ktse ] = f_TSE ;
	  //(D.f[TSW ])[ktsw ] = f_TSW ;
	  //(D.f[BSE ])[kbse ] = f_BSE ;
	  //(D.f[BSW ])[kbsw ] = f_BSW ;     
	  // Z
	  //(D.f[B   ])[kb   ] = f_B   ;
	  //(D.f[BE  ])[kbe  ] = f_BE  ;
	  //(D.f[BW  ])[kbw  ] = f_BW  ;
	  //(D.f[BN  ])[kbn  ] = f_BN  ;
	  //(D.f[BS  ])[kbs  ] = f_BS  ;
	  //(D.f[BNE ])[kbne ] = f_BNE ;
	  //(D.f[BNW ])[kbnw ] = f_BNW ;
	  //(D.f[BSE ])[kbse ] = f_BSE ;
	  //(D.f[BSW ])[kbsw ] = f_BSW ;     
      //////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceIncompNEQ27( real* rhoBC,
													real* DD, 
													int* k_Q, 
													int* k_N, 
													int numberOfBCnodes, 
													real om1, 
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
													bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int K1QK  = k_N[k];
      unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==true) //// ACHTUNG PREColl !!!!!!!!!!!!!!
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[E   ])[k1e   ];
      f1_E    = (D.f[W   ])[k1w   ];
      f1_S    = (D.f[N   ])[k1n   ];
      f1_N    = (D.f[S   ])[k1s   ];
      f1_B    = (D.f[T   ])[k1t   ];
      f1_T    = (D.f[B   ])[k1b   ];
      f1_SW   = (D.f[NE  ])[k1ne  ];
      f1_NE   = (D.f[SW  ])[k1sw  ];
      f1_NW   = (D.f[SE  ])[k1se  ];
      f1_SE   = (D.f[NW  ])[k1nw  ];
      f1_BW   = (D.f[TE  ])[k1te  ];
      f1_TE   = (D.f[BW  ])[k1bw  ];
      f1_TW   = (D.f[BE  ])[k1be  ];
      f1_BE   = (D.f[TW  ])[k1tw  ];
      f1_BS   = (D.f[TN  ])[k1tn  ];
      f1_TN   = (D.f[BS  ])[k1bs  ];
      f1_TS   = (D.f[BN  ])[k1bn  ];
      f1_BN   = (D.f[TS  ])[k1ts  ];
      f1_ZERO = (D.f[REST])[k1zero];
      f1_BSW  = (D.f[TNE ])[k1tne ];
      f1_BNE  = (D.f[TSW ])[k1tsw ];
      f1_BNW  = (D.f[TSE ])[k1tse ];
      f1_BSE  = (D.f[TNW ])[k1tnw ];
      f1_TSW  = (D.f[BNE ])[k1bne ];
      f1_TNE  = (D.f[BSW ])[k1bsw ];
      f1_TNW  = (D.f[BSE ])[k1bse ];
      f1_TSE  = (D.f[BNW ])[k1bnw ];

      //////////////////////////////////////////////////////////////////////////
      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                          f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      real vx1      =  ((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
						  ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
						  (f1_E - f1_W); 


      real vx2    =   (-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
						 ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
						 (f1_N - f1_S); 

      real vx3    =   ((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
						 (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
						 (f1_T - f1_B); 

      real cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      f1_ZERO  -= c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
	   
	  drho1 = rhoBC[k];

	  //if(vx1 < zero){
		 // vx1 *= 0.9;
	  //}
	  //if(vx2 < zero){
		 // vx2 *= c1o10;//0.9;
	  //}

      f1_ZERO  += c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

	  //drho1 = (drho1 + rhoBC[k])/2.f;
	  //drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[E   ])[ke   ] = f1_W   ;  
      (D.f[W   ])[kw   ] = f1_E   ;	
      (D.f[N   ])[kn   ] = f1_S   ;	
      (D.f[S   ])[ks   ] = f1_N   ;	
      (D.f[T   ])[kt   ] = f1_B   ;	
      (D.f[B   ])[kb   ] = f1_T   ;	
      (D.f[NE  ])[kne  ] = f1_SW  ;	
      (D.f[SW  ])[ksw  ] = f1_NE  ;	
      (D.f[SE  ])[kse  ] = f1_NW  ;	
      (D.f[NW  ])[knw  ] = f1_SE  ;	
      (D.f[TE  ])[kte  ] = f1_BW  ;	
      (D.f[BW  ])[kbw  ] = f1_TE  ;	
      (D.f[BE  ])[kbe  ] = f1_TW  ;	
      (D.f[TW  ])[ktw  ] = f1_BE  ;	
      (D.f[TN  ])[ktn  ] = f1_BS  ;	
      (D.f[BS  ])[kbs  ] = f1_TN  ;	
      (D.f[BN  ])[kbn  ] = f1_TS  ;	
      (D.f[TS  ])[kts  ] = f1_BN  ;	
      (D.f[REST])[kzero] = f1_ZERO;	
      (D.f[TNE ])[ktne ] = f1_BSW ;	
      (D.f[TSW ])[ktsw ] = f1_BNE ;	
      (D.f[TSE ])[ktse ] = f1_BNW ;	
      (D.f[TNW ])[ktnw ] = f1_BSE ;	
      (D.f[BNE ])[kbne ] = f1_TSW ;	
      (D.f[BSW ])[kbsw ] = f1_TNE ;	
      (D.f[BSE ])[kbse ] = f1_TNW ;	
      (D.f[BNW ])[kbnw ] = f1_TSE ;       
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceNEQ27(real* rhoBC,
                                             real* distribution, 
                                             int* bcNodeIndices,
                                             int* bcNeighborIndices,
                                             int numberOfBCnodes,
                                             real omega1, 
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int numberOfLBnodes, 
                                             bool isEvenTimestep)
{
   //////////////////////////////////////////////////////////////////////////
	//! The pressure boundary condition is executed in the following steps
	//!
	////////////////////////////////////////////////////////////////////////////////
	//! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
	//!
   const unsigned x = threadIdx.x;    // global x-index 
   const unsigned y = blockIdx.x;     // global y-index 
   const unsigned z = blockIdx.y;     // global z-index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;

   //////////////////////////////////////////////////////////////////////////
   //! - Run for all indices in size of boundary condition (numberOfBCnodes)
   //!
   if(k < numberOfBCnodes)
   {
      //////////////////////////////////////////////////////////////////////////
      //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep is based on the esoteric twist algorithm \ref
      //! <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
      //!
      Distributions27 dist;
      getPointersToDistributions(dist, distribution, numberOfLBnodes, isEvenTimestep);

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local pressure
      //!
      real rhoBClocal = rhoBC[k];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing)
      //!
      unsigned int KQK  = bcNodeIndices[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //! - Set neighbor indices (necessary for indirect addressing) for neighboring node
      //!
      unsigned int K1QK  = bcNeighborIndices[k];
      unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set local distributions for neighboring node
      //!
      real f1_W    = (dist.f[E   ])[k1e   ];
      real f1_E    = (dist.f[W   ])[k1w   ];
      real f1_S    = (dist.f[N   ])[k1n   ];
      real f1_N    = (dist.f[S   ])[k1s   ];
      real f1_B    = (dist.f[T   ])[k1t   ];
      real f1_T    = (dist.f[B   ])[k1b   ];
      real f1_SW   = (dist.f[NE  ])[k1ne  ];
      real f1_NE   = (dist.f[SW  ])[k1sw  ];
      real f1_NW   = (dist.f[SE  ])[k1se  ];
      real f1_SE   = (dist.f[NW  ])[k1nw  ];
      real f1_BW   = (dist.f[TE  ])[k1te  ];
      real f1_TE   = (dist.f[BW  ])[k1bw  ];
      real f1_TW   = (dist.f[BE  ])[k1be  ];
      real f1_BE   = (dist.f[TW  ])[k1tw  ];
      real f1_BS   = (dist.f[TN  ])[k1tn  ];
      real f1_TN   = (dist.f[BS  ])[k1bs  ];
      real f1_TS   = (dist.f[BN  ])[k1bn  ];
      real f1_BN   = (dist.f[TS  ])[k1ts  ];
      real f1_ZERO = (dist.f[REST])[k1zero];
      real f1_BSW  = (dist.f[TNE ])[k1tne ];
      real f1_BNE  = (dist.f[TSW ])[k1tsw ];
      real f1_BNW  = (dist.f[TSE ])[k1tse ];
      real f1_BSE  = (dist.f[TNW ])[k1tnw ];
      real f1_TSW  = (dist.f[BNE ])[k1bne ];
      real f1_TNE  = (dist.f[BSW ])[k1bsw ];
      real f1_TNW  = (dist.f[BSE ])[k1bse ];
      real f1_TSE  = (dist.f[BNW ])[k1bnw ];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Calculate macroscopic quantities (for neighboring node)
      //!
      real drho1 = f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                   f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW + 
                   f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((dist.f[REST])[kzero]); 

      real vx1  = (((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                   ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
                   (f1_E - f1_W)) / (c1o1 + drho1);          

      real vx2  = ((-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                   ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
                   (f1_N - f1_S)) / (c1o1 + drho1); 

      real vx3  = (((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
                   (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
                   (f1_T - f1_B)) / (c1o1 + drho1); 

      real cusq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

      ////////////////////////////////////////////////////////////////////////////////
      //! subtract the equilibrium (eq) to obtain the non-equilibrium (neq) (for neighboring node)
      //!
      f1_ZERO  -= c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     -= c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    -= c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    -=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   -=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

      ////////////////////////////////////////////////////////////////////////////////
      //! redefine drho1 with rhoBClocal
      //!
      drho1 = rhoBClocal;

      ////////////////////////////////////////////////////////////////////////////////
      //! add the equilibrium (eq), which is calculated with rhoBClocal (for neighboring node)
      //!
      f1_ZERO  += c8o27*  (drho1-(drho1+c1o1)*cusq);
      f1_E     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f1_W     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f1_N     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f1_S     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f1_B     += c2o27*  (drho1+(drho1+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    += c1o54*  (drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    +=  c1o54* (drho1+(drho1+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   +=  c1o216*(drho1+(drho1+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      ////////////////////////////////////////////////////////////////////////////////
      //! write the new distributions to the bc nodes
      //!
      (dist.f[E   ])[ke   ] = f1_W   ;
      (dist.f[W   ])[kw   ] = f1_E   ;
      (dist.f[N   ])[kn   ] = f1_S   ;
      (dist.f[S   ])[ks   ] = f1_N   ;
      (dist.f[T   ])[kt   ] = f1_B   ;
      (dist.f[B   ])[kb   ] = f1_T   ;
      (dist.f[NE  ])[kne  ] = f1_SW  ;
      (dist.f[SW  ])[ksw  ] = f1_NE  ;
      (dist.f[SE  ])[kse  ] = f1_NW  ;
      (dist.f[NW  ])[knw  ] = f1_SE  ;
      (dist.f[TE  ])[kte  ] = f1_BW  ;
      (dist.f[BW  ])[kbw  ] = f1_TE  ;
      (dist.f[BE  ])[kbe  ] = f1_TW  ;
      (dist.f[TW  ])[ktw  ] = f1_BE  ;
      (dist.f[TN  ])[ktn  ] = f1_BS  ;
      (dist.f[BS  ])[kbs  ] = f1_TN  ;
      (dist.f[BN  ])[kbn  ] = f1_TS  ;
      (dist.f[TS  ])[kts  ] = f1_BN  ;
      (dist.f[REST])[kzero] = f1_ZERO;
      (dist.f[TNE ])[ktne ] = f1_BSW ;
      (dist.f[TSW ])[ktsw ] = f1_BNE ;
      (dist.f[TSE ])[ktse ] = f1_BNW ;
      (dist.f[TNW ])[ktnw ] = f1_BSE ;
      (dist.f[BNE ])[kbne ] = f1_TSW ;
      (dist.f[BSW ])[kbsw ] = f1_TNE ;
      (dist.f[BSE ])[kbse ] = f1_TNW ;
      (dist.f[BNW ])[kbnw ] = f1_TSE ;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void LB_BC_Press_East27( int nx, 
                                               int ny, 
                                               int tz, 
                                               unsigned int* bcMatD, 
                                               unsigned int* neighborX,
                                               unsigned int* neighborY,
                                               unsigned int* neighborZ,
                                               real* DD, 
                                               unsigned int size_Mat, 
                                               bool isEvenTimestep) 
{
   //thread-index
   int ty = blockIdx.x;
   int tx = threadIdx.x;

   int  k, k1, nxny;                   // Zugriff auf arrays im device

   int  x = tx + STARTOFFX;  // Globaler x-Index 
   int  y = ty + STARTOFFY;  // Globaler y-Index 
   int  z = tz + STARTOFFZ;  // Globaler z-Index 

   k = nx*(ny*z + y) + x;
   nxny = nx*ny;
   k1 = k-nxny;

   if( bcMatD[k] == GEO_PRESS && bcMatD[k1] == GEO_FLUID)
   {
      Distributions27 D;
      if (isEvenTimestep==true)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      //index
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
      //unsigned int kzero= k;
      //unsigned int ke   = k;
      //unsigned int kw   = k + 1;
      //unsigned int kn   = k;
      //unsigned int ks   = k + nx;
      //unsigned int kt   = k;
      //unsigned int kb   = k + nxny;
      //unsigned int ksw  = k + nx + 1;
      //unsigned int kne  = k;
      //unsigned int kse  = k + nx;
      //unsigned int knw  = k + 1;
      //unsigned int kbw  = k + nxny + 1;
      //unsigned int kte  = k;
      //unsigned int kbe  = k + nxny;
      //unsigned int ktw  = k + 1;
      //unsigned int kbs  = k + nxny + nx;
      //unsigned int ktn  = k;
      //unsigned int kbn  = k + nxny;
      //unsigned int kts  = k + nx;
      //unsigned int ktse = k + nx;
      //unsigned int kbnw = k + nxny + 1;
      //unsigned int ktnw = k + 1;
      //unsigned int kbse = k + nxny + nx;
      //unsigned int ktsw = k + nx + 1;
      //unsigned int kbne = k + nxny;
      //unsigned int ktne = k;
      //unsigned int kbsw = k + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int k1zero= k1;
      unsigned int k1e   = k1;
      unsigned int k1w   = neighborX[k1];
      unsigned int k1n   = k1;
      unsigned int k1s   = neighborY[k1];
      unsigned int k1t   = k1;
      unsigned int k1b   = neighborZ[k1];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = k1;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = k1;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = k1;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = k1;
      unsigned int k1bsw = neighborZ[k1sw];
      //unsigned int k1zero= k1;
      //unsigned int k1e   = k1;
      //unsigned int k1w   = k1 + 1;
      //unsigned int k1n   = k1;
      //unsigned int k1s   = k1 + nx;
      //unsigned int k1t   = k1;
      //unsigned int k1b   = k1 + nxny;
      //unsigned int k1sw  = k1 + nx + 1;
      //unsigned int k1ne  = k1;
      //unsigned int k1se  = k1 + nx;
      //unsigned int k1nw  = k1 + 1;
      //unsigned int k1bw  = k1 + nxny + 1;
      //unsigned int k1te  = k1;
      //unsigned int k1be  = k1 + nxny;
      //unsigned int k1tw  = k1 + 1;
      //unsigned int k1bs  = k1 + nxny + nx;
      //unsigned int k1tn  = k1;
      //unsigned int k1bn  = k1 + nxny;
      //unsigned int k1ts  = k1 + nx;
      //unsigned int k1tse = k1 + nx;
      //unsigned int k1bnw = k1 + nxny + 1;
      //unsigned int k1tnw = k1 + 1;
      //unsigned int k1bse = k1 + nxny + nx;
      //unsigned int k1tsw = k1 + nx + 1;
      //unsigned int k1bne = k1 + nxny;
      //unsigned int k1tne = k1;
      //unsigned int k1bsw = k1 + nxny + nx + 1;
      ////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                   f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[E   ])[k1e   ];
      f1_E    = (D.f[W   ])[k1w   ];
      f1_S    = (D.f[N   ])[k1n   ];
      f1_N    = (D.f[S   ])[k1s   ];
      f1_B    = (D.f[T   ])[k1t   ];
      f1_T    = (D.f[B   ])[k1b   ];
      f1_SW   = (D.f[NE  ])[k1ne  ];
      f1_NE   = (D.f[SW  ])[k1sw  ];
      f1_NW   = (D.f[SE  ])[k1se  ];
      f1_SE   = (D.f[NW  ])[k1nw  ];
      f1_BW   = (D.f[TE  ])[k1te  ];
      f1_TE   = (D.f[BW  ])[k1bw  ];
      f1_TW   = (D.f[BE  ])[k1be  ];
      f1_BE   = (D.f[TW  ])[k1tw  ];
      f1_BS   = (D.f[TN  ])[k1tn  ];
      f1_TN   = (D.f[BS  ])[k1bs  ];
      f1_TS   = (D.f[BN  ])[k1bn  ];
      f1_BN   = (D.f[TS  ])[k1ts  ];
      f1_ZERO = (D.f[REST])[k1zero];
      f1_BSW  = (D.f[TNE ])[k1tne ];
      f1_BNE  = (D.f[TSW ])[k1tsw ];
      f1_BNW  = (D.f[TSE ])[k1tse ];
      f1_BSE  = (D.f[TNW ])[k1tnw ];
      f1_TSW  = (D.f[BNE ])[k1bne ];
      f1_TNE  = (D.f[BSW ])[k1bsw ];
      f1_TNW  = (D.f[BSE ])[k1bse ];
      f1_TSE  = (D.f[BNW ])[k1bnw ];

      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                        f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      __syncthreads();

      (D.f[E   ])[ke   ] = f1_W   -c2o27*drho1;
      (D.f[W   ])[kw   ] = f1_E   -c2o27*drho1;
      (D.f[N   ])[kn   ] = f1_S   -c2o27*drho1;
      (D.f[S   ])[ks   ] = f1_N   -c2o27*drho1;
      (D.f[T   ])[kt   ] = f1_B   -c2o27*drho1;
      (D.f[B   ])[kb   ] = f1_T   -c2o27*drho1;
      (D.f[NE  ])[kne  ] = f1_SW  -c1o54*drho1;
      (D.f[SW  ])[ksw  ] = f1_NE  -c1o54*drho1;
      (D.f[SE  ])[kse  ] = f1_NW  -c1o54*drho1;
      (D.f[NW  ])[knw  ] = f1_SE  -c1o54*drho1;
      (D.f[TE  ])[kte  ] = f1_BW  -c1o54*drho1;
      (D.f[BW  ])[kbw  ] = f1_TE  -c1o54*drho1;
      (D.f[BE  ])[kbe  ] = f1_TW  -c1o54*drho1;
      (D.f[TW  ])[ktw  ] = f1_BE  -c1o54*drho1;
      (D.f[TN  ])[ktn  ] = f1_BS  -c1o54*drho1;
      (D.f[BS  ])[kbs  ] = f1_TN  -c1o54*drho1;
      (D.f[BN  ])[kbn  ] = f1_TS  -c1o54*drho1;
      (D.f[TS  ])[kts  ] = f1_BN  -c1o54*drho1;
      (D.f[REST])[kzero] = f1_ZERO-c8o27*drho1;
      (D.f[TNE ])[ktne ] = f1_BSW -c1o216*drho1;
      (D.f[TSW ])[ktsw ] = f1_BNE -c1o216*drho1;
      (D.f[TSE ])[ktse ] = f1_BNW -c1o216*drho1;
      (D.f[TNW ])[ktnw ] = f1_BSE -c1o216*drho1;
      (D.f[BNE ])[kbne ] = f1_TSW -c1o216*drho1;
      (D.f[BSW ])[kbsw ] = f1_TNE -c1o216*drho1;
      (D.f[BSE ])[kbse ] = f1_TNW -c1o216*drho1;
      (D.f[BNW ])[kbnw ] = f1_TSE -c1o216*drho1;       
   }
   __syncthreads();
}          
//////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDevice27(real* rhoBC,
                                           real* DD, 
                                           int* k_Q, 
                                           real* QQ,
                                           unsigned int numberOfBCnodes, 
                                           real om1, 
                                           unsigned int* neighborX,
                                           unsigned int* neighborY,
                                           unsigned int* neighborZ,
                                           unsigned int size_Mat, 
                                           bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[E   ] = &DD[E   *size_Mat];
      D.f[W   ] = &DD[W   *size_Mat];
      D.f[N   ] = &DD[N   *size_Mat];
      D.f[S   ] = &DD[S   *size_Mat];
      D.f[T   ] = &DD[T   *size_Mat];
      D.f[B   ] = &DD[B   *size_Mat];
      D.f[NE  ] = &DD[NE  *size_Mat];
      D.f[SW  ] = &DD[SW  *size_Mat];
      D.f[SE  ] = &DD[SE  *size_Mat];
      D.f[NW  ] = &DD[NW  *size_Mat];
      D.f[TE  ] = &DD[TE  *size_Mat];
      D.f[BW  ] = &DD[BW  *size_Mat];
      D.f[BE  ] = &DD[BE  *size_Mat];
      D.f[TW  ] = &DD[TW  *size_Mat];
      D.f[TN  ] = &DD[TN  *size_Mat];
      D.f[BS  ] = &DD[BS  *size_Mat];
      D.f[BN  ] = &DD[BN  *size_Mat];
      D.f[TS  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[TNE *size_Mat];
      D.f[TSW ] = &DD[TSW *size_Mat];
      D.f[TSE ] = &DD[TSE *size_Mat];
      D.f[TNW ] = &DD[TNW *size_Mat];
      D.f[BNE ] = &DD[BNE *size_Mat];
      D.f[BSW ] = &DD[BSW *size_Mat];
      D.f[BSE ] = &DD[BSE *size_Mat];
      D.f[BNW ] = &DD[BNW *size_Mat];
   } 
   else
   {
      D.f[W   ] = &DD[E   *size_Mat];
      D.f[E   ] = &DD[W   *size_Mat];
      D.f[S   ] = &DD[N   *size_Mat];
      D.f[N   ] = &DD[S   *size_Mat];
      D.f[B   ] = &DD[T   *size_Mat];
      D.f[T   ] = &DD[B   *size_Mat];
      D.f[SW  ] = &DD[NE  *size_Mat];
      D.f[NE  ] = &DD[SW  *size_Mat];
      D.f[NW  ] = &DD[SE  *size_Mat];
      D.f[SE  ] = &DD[NW  *size_Mat];
      D.f[BW  ] = &DD[TE  *size_Mat];
      D.f[TE  ] = &DD[BW  *size_Mat];
      D.f[TW  ] = &DD[BE  *size_Mat];
      D.f[BE  ] = &DD[TW  *size_Mat];
      D.f[BS  ] = &DD[TN  *size_Mat];
      D.f[TN  ] = &DD[BS  *size_Mat];
      D.f[TS  ] = &DD[BN  *size_Mat];
      D.f[BN  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   * numberOfBCnodes];
      q_dirW   = &QQ[W   * numberOfBCnodes];
      q_dirN   = &QQ[N   * numberOfBCnodes];
      q_dirS   = &QQ[S   * numberOfBCnodes];
      q_dirT   = &QQ[T   * numberOfBCnodes];
      q_dirB   = &QQ[B   * numberOfBCnodes];
      q_dirNE  = &QQ[NE  * numberOfBCnodes];
      q_dirSW  = &QQ[SW  * numberOfBCnodes];
      q_dirSE  = &QQ[SE  * numberOfBCnodes];
      q_dirNW  = &QQ[NW  * numberOfBCnodes];
      q_dirTE  = &QQ[TE  * numberOfBCnodes];
      q_dirBW  = &QQ[BW  * numberOfBCnodes];
      q_dirBE  = &QQ[BE  * numberOfBCnodes];
      q_dirTW  = &QQ[TW  * numberOfBCnodes];
      q_dirTN  = &QQ[TN  * numberOfBCnodes];
      q_dirBS  = &QQ[BS  * numberOfBCnodes];
      q_dirBN  = &QQ[BN  * numberOfBCnodes];
      q_dirTS  = &QQ[TS  * numberOfBCnodes];
      q_dirTNE = &QQ[TNE * numberOfBCnodes];
      q_dirTSW = &QQ[TSW * numberOfBCnodes];
      q_dirTSE = &QQ[TSE * numberOfBCnodes];
      q_dirTNW = &QQ[TNW * numberOfBCnodes];
      q_dirBNE = &QQ[BNE * numberOfBCnodes];
      q_dirBSW = &QQ[BSW * numberOfBCnodes];
      q_dirBSE = &QQ[BSE * numberOfBCnodes];
      q_dirBNW = &QQ[BNW * numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[E   ])[ke   ];
      f_E    = (D.f[W   ])[kw   ];
      f_S    = (D.f[N   ])[kn   ];
      f_N    = (D.f[S   ])[ks   ];
      f_B    = (D.f[T   ])[kt   ];
      f_T    = (D.f[B   ])[kb   ];
      f_SW   = (D.f[NE  ])[kne  ];
      f_NE   = (D.f[SW  ])[ksw  ];
      f_NW   = (D.f[SE  ])[kse  ];
      f_SE   = (D.f[NW  ])[knw  ];
      f_BW   = (D.f[TE  ])[kte  ];
      f_TE   = (D.f[BW  ])[kbw  ];
      f_TW   = (D.f[BE  ])[kbe  ];
      f_BE   = (D.f[TW  ])[ktw  ];
      f_BS   = (D.f[TN  ])[ktn  ];
      f_TN   = (D.f[BS  ])[kbs  ];
      f_TS   = (D.f[BN  ])[kbn  ];
      f_BN   = (D.f[TS  ])[kts  ];
      f_BSW  = (D.f[TNE ])[ktne ];
      f_BNE  = (D.f[TSW ])[ktsw ];
      f_BNW  = (D.f[TSE ])[ktse ];
      f_BSE  = (D.f[TNW ])[ktnw ];
      f_TSW  = (D.f[BNE ])[kbne ];
      f_TNE  = (D.f[BSW ])[kbsw ];
      f_TNW  = (D.f[BSE ])[kbse ];
      f_TSE  = (D.f[BNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      real q, vx1, vx2, vx3, drho;
      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                  (f_E - f_W); 


      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S); 

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      drho = rhoBC[k];
      //deltaRho = (rhoBC[k] + one) / (deltaRho + one);
      ////////////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[W])[kw]=c2o27* (drho+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq); 
         //(D.f[E])[ke]=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[E])[ke]=c2o27* (drho+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq); 
         //(D.f[W])[kw]=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[S])[ks]=c2o27* (drho+c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); 
         //(D.f[N])[kn]=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[N])[kn]=c2o27* (drho+c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cu_sq); 
         //(D.f[S])[ks]=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[B])[kb]=c2o27* (drho+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
         //(D.f[T])[kt]=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[T])[kt]=c2o27* (drho+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq); 
         //(D.f[B])[kb]=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[SW])[ksw]=c1o54* (drho+c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         //(D.f[NE])[kne]=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[NE])[kne]=c1o54* (drho+c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         //(D.f[SW])[ksw]=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[NW])[knw]=c1o54* (drho+c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         //(D.f[SE])[kse]=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[SE])[kse]=c1o54* (drho+c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         //(D.f[NW])[knw]=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BW])[kbw]=c1o54* (drho+c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         //(D.f[TE])[kte]=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TE])[kte]=c1o54* (drho+c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         //(D.f[BW])[kbw]=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TW])[ktw]=c1o54* (drho+c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         //(D.f[BE])[kbe]=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BE])[kbe]=c1o54* (drho+c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         //(D.f[TW])[ktw]=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BS])[kbs]=c1o54* (drho+c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         //(D.f[TN])[ktn]=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TN])[ktn]=c1o54* (drho+c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         //(D.f[BS])[kbs]=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TS])[kts]=c1o54* (drho+c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         //(D.f[BN])[kbn]=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BN])[kbn]=c1o54* (drho+c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         //(D.f[TS])[kts]=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BSW])[kbsw]=c1o216*(drho+c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         //(D.f[TNE])[ktne]=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TNE])[ktne]=c1o216*(drho+c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         //(D.f[BSW])[kbsw]=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TSW])[ktsw]=c1o216*(drho+c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         //(D.f[BNE])[kbne]=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BNE])[kbne]=c1o216*(drho+c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         //(D.f[TSW])[ktsw]=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BNW])[kbnw]=c1o216*(drho+c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         //(D.f[TSE])[ktse]=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TSE])[ktse]=c1o216*(drho+c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         //(D.f[BNW])[kbnw]=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TNW])[ktnw]=c1o216*(drho+c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         //(D.f[BSE])[kbse]=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BSE])[kbse]=c1o216*(drho+c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         //(D.f[TNW])[ktnw]=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceAntiBB27(   real* rhoBC,
												   real* vx,
												   real* vy,
												   real* vz,
												   real* DD, 
												   int* k_Q, 
												   real* QQ,
												   int numberOfBCnodes, 
												   real om1, 
												   unsigned int* neighborX,
												   unsigned int* neighborY,
												   unsigned int* neighborZ,
												   unsigned int size_Mat, 
												   bool isEvenTimestep)
{
   Distributions27 D;
   if (isEvenTimestep==true)
   {
      D.f[E   ] = &DD[E   *size_Mat];
      D.f[W   ] = &DD[W   *size_Mat];
      D.f[N   ] = &DD[N   *size_Mat];
      D.f[S   ] = &DD[S   *size_Mat];
      D.f[T   ] = &DD[T   *size_Mat];
      D.f[B   ] = &DD[B   *size_Mat];
      D.f[NE  ] = &DD[NE  *size_Mat];
      D.f[SW  ] = &DD[SW  *size_Mat];
      D.f[SE  ] = &DD[SE  *size_Mat];
      D.f[NW  ] = &DD[NW  *size_Mat];
      D.f[TE  ] = &DD[TE  *size_Mat];
      D.f[BW  ] = &DD[BW  *size_Mat];
      D.f[BE  ] = &DD[BE  *size_Mat];
      D.f[TW  ] = &DD[TW  *size_Mat];
      D.f[TN  ] = &DD[TN  *size_Mat];
      D.f[BS  ] = &DD[BS  *size_Mat];
      D.f[BN  ] = &DD[BN  *size_Mat];
      D.f[TS  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[TNE *size_Mat];
      D.f[TSW ] = &DD[TSW *size_Mat];
      D.f[TSE ] = &DD[TSE *size_Mat];
      D.f[TNW ] = &DD[TNW *size_Mat];
      D.f[BNE ] = &DD[BNE *size_Mat];
      D.f[BSW ] = &DD[BSW *size_Mat];
      D.f[BSE ] = &DD[BSE *size_Mat];
      D.f[BNW ] = &DD[BNW *size_Mat];
   } 
   else
   {
      D.f[W   ] = &DD[E   *size_Mat];
      D.f[E   ] = &DD[W   *size_Mat];
      D.f[S   ] = &DD[N   *size_Mat];
      D.f[N   ] = &DD[S   *size_Mat];
      D.f[B   ] = &DD[T   *size_Mat];
      D.f[T   ] = &DD[B   *size_Mat];
      D.f[SW  ] = &DD[NE  *size_Mat];
      D.f[NE  ] = &DD[SW  *size_Mat];
      D.f[NW  ] = &DD[SE  *size_Mat];
      D.f[SE  ] = &DD[NW  *size_Mat];
      D.f[BW  ] = &DD[TE  *size_Mat];
      D.f[TE  ] = &DD[BW  *size_Mat];
      D.f[TW  ] = &DD[BE  *size_Mat];
      D.f[BE  ] = &DD[TW  *size_Mat];
      D.f[BS  ] = &DD[TN  *size_Mat];
      D.f[TN  ] = &DD[BS  *size_Mat];
      D.f[TS  ] = &DD[BN  *size_Mat];
      D.f[BN  ] = &DD[TS  *size_Mat];
      D.f[REST] = &DD[REST*size_Mat];
      D.f[TNE ] = &DD[BSW *size_Mat];
      D.f[TSW ] = &DD[BNE *size_Mat];
      D.f[TSE ] = &DD[BNW *size_Mat];
      D.f[TNW ] = &DD[BSE *size_Mat];
      D.f[BNE ] = &DD[TSW *size_Mat];
      D.f[BSW ] = &DD[TNE *size_Mat];
      D.f[BSE ] = &DD[TNW *size_Mat];
      D.f[BNW ] = &DD[TSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[E   *numberOfBCnodes];
      q_dirW   = &QQ[W   *numberOfBCnodes];
      q_dirN   = &QQ[N   *numberOfBCnodes];
      q_dirS   = &QQ[S   *numberOfBCnodes];
      q_dirT   = &QQ[T   *numberOfBCnodes];
      q_dirB   = &QQ[B   *numberOfBCnodes];
      q_dirNE  = &QQ[NE  *numberOfBCnodes];
      q_dirSW  = &QQ[SW  *numberOfBCnodes];
      q_dirSE  = &QQ[SE  *numberOfBCnodes];
      q_dirNW  = &QQ[NW  *numberOfBCnodes];
      q_dirTE  = &QQ[TE  *numberOfBCnodes];
      q_dirBW  = &QQ[BW  *numberOfBCnodes];
      q_dirBE  = &QQ[BE  *numberOfBCnodes];
      q_dirTW  = &QQ[TW  *numberOfBCnodes];
      q_dirTN  = &QQ[TN  *numberOfBCnodes];
      q_dirBS  = &QQ[BS  *numberOfBCnodes];
      q_dirBN  = &QQ[BN  *numberOfBCnodes];
      q_dirTS  = &QQ[TS  *numberOfBCnodes];
      q_dirTNE = &QQ[TNE *numberOfBCnodes];
      q_dirTSW = &QQ[TSW *numberOfBCnodes];
      q_dirTSE = &QQ[TSE *numberOfBCnodes];
      q_dirTNW = &QQ[TNW *numberOfBCnodes];
      q_dirBNE = &QQ[BNE *numberOfBCnodes];
      q_dirBSW = &QQ[BSW *numberOfBCnodes];
      q_dirBSE = &QQ[BSE *numberOfBCnodes];
      q_dirBNW = &QQ[BNW *numberOfBCnodes];
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW, f_ZERO;

      f_W    = (D.f[E   ])[ke   ];
      f_E    = (D.f[W   ])[kw   ];
      f_S    = (D.f[N   ])[kn   ];
      f_N    = (D.f[S   ])[ks   ];
      f_B    = (D.f[T   ])[kt   ];
      f_T    = (D.f[B   ])[kb   ];
      f_SW   = (D.f[NE  ])[kne  ];
      f_NE   = (D.f[SW  ])[ksw  ];
      f_NW   = (D.f[SE  ])[kse  ];
      f_SE   = (D.f[NW  ])[knw  ];
      f_BW   = (D.f[TE  ])[kte  ];
      f_TE   = (D.f[BW  ])[kbw  ];
      f_TW   = (D.f[BE  ])[kbe  ];
      f_BE   = (D.f[TW  ])[ktw  ];
      f_BS   = (D.f[TN  ])[ktn  ];
      f_TN   = (D.f[BS  ])[kbs  ];
      f_TS   = (D.f[BN  ])[kbn  ];
      f_BN   = (D.f[TS  ])[kts  ];
      f_BSW  = (D.f[TNE ])[ktne ];
      f_BNE  = (D.f[TSW ])[ktsw ];
      f_BNW  = (D.f[TSE ])[ktse ];
      f_BSE  = (D.f[TNW ])[ktnw ];
      f_TSW  = (D.f[BNE ])[kbne ];
      f_TNE  = (D.f[BSW ])[kbsw ];
      f_TNW  = (D.f[BSE ])[kbse ];
      f_TSE  = (D.f[BNW ])[kbnw ];
      f_ZERO = (D.f[REST])[kzero];
      ////////////////////////////////////////////////////////////////////////////////
      //real vx1, vx2, vx3, drho;
      //vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //            ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //            (f_E - f_W); 


      //vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //            ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //            (f_N - f_S); 

      //vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //            (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //            (f_T - f_B); 

      //real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      real drho    = f_ZERO+f_E+f_W+f_N+f_S+f_T+f_B+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+
						f_TNE+f_TSW+f_TSE+f_TNW+f_BNE+f_BSW+f_BSE+f_BNW;
      drho = drho - rhoBC[k];
	  drho *= 0.01f;
      ////////////////////////////////////////////////////////////////////////////////
	  real q;
      //deltaRho = (rhoBC[k] + one) / (deltaRho + one);
      ////////////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[W])[kw]=f_W-c2o27*drho; 
      }

      q = q_dirW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[E])[ke]=f_E-c2o27*drho;
      }

      q = q_dirN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[S])[ks]=f_S-c2o27*drho; 
      }

      q = q_dirS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[N])[kn]=f_N-c2o27*drho; 
      }

      q = q_dirT[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[B])[kb]=f_B-c2o27*drho; 
      }

      q = q_dirB[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[T])[kt]=f_T-c2o27*drho;
      }

      q = q_dirNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[SW])[ksw]=f_SW-c1o54*drho;
      }

      q = q_dirSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[NE])[kne]=f_NE-c1o54*drho;
      }

      q = q_dirSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[NW])[knw]=f_NW-c1o54*drho;
      }

      q = q_dirNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[SE])[kse]=f_SE-c1o54*drho;
      }

      q = q_dirTE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BW])[kbw]=f_BW-c1o54*drho; 
      }

      q = q_dirBW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TE])[kte]=f_TE-c1o54*drho; 
      }

      q = q_dirBE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TW])[ktw]=f_TW-c1o54*drho;
      }

      q = q_dirTW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BE])[kbe]=f_BE-c1o54*drho;
      }

      q = q_dirTN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BS])[kbs]=f_BS-c1o54*drho;
      }

      q = q_dirBS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TN])[ktn]=f_TN-c1o54*drho;
      }

      q = q_dirBN[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TS])[kts]=f_TS-c1o54*drho;
      }

      q = q_dirTS[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BN])[kbn]=f_BN-c1o54*drho;
      }

      q = q_dirTNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BSW])[kbsw]=f_BSW-c1o216*drho;
      }

      q = q_dirBSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TNE])[ktne]=f_TNE-c1o216*drho;
      }

      q = q_dirBNE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TSW])[ktsw]=f_TSW-c1o216*drho;
      }

      q = q_dirTSW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BNE])[kbne]=f_BNE-c1o216*drho;
      }

      q = q_dirTSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BNW])[kbnw]=f_BNW-c1o216*drho;
      }

      q = q_dirBNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TSE])[ktse]=f_TSE-c1o216*drho;
      }

      q = q_dirBSE[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[TNW])[ktnw]=f_TNW-c1o216*drho;
      }

      q = q_dirTNW[k];
      if (q>=c0o1 && q<=c1o1)
      {
         (D.f[BSE])[kbse]=f_BSE-c1o216*drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceFixBackflow27( real* rhoBC,
                                                      real* DD, 
                                                      int* k_Q, 
                                                      int numberOfBCnodes, 
                                                      real om1, 
                                                      unsigned int* neighborX,
                                                      unsigned int* neighborY,
                                                      unsigned int* neighborZ,
                                                      unsigned int size_Mat, 
                                                      bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real deltaRho;
      ////////////////////////////////////////////////////////////////////////////////
      deltaRho = rhoBC[k];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         (D.f[W])[kw]       = c2o27  * deltaRho;
         (D.f[E])[ke]       = c2o27  * deltaRho;
         (D.f[S])[ks]       = c2o27  * deltaRho;
         (D.f[N])[kn]       = c2o27  * deltaRho;
         (D.f[B])[kb]       = c2o27  * deltaRho;
         (D.f[T])[kt]       = c2o27  * deltaRho;
         (D.f[SW])[ksw]     = c1o54  * deltaRho;
         (D.f[NE])[kne]     = c1o54  * deltaRho;
         (D.f[NW])[knw]     = c1o54  * deltaRho;
         (D.f[SE])[kse]     = c1o54  * deltaRho;
         (D.f[BW])[kbw]     = c1o54  * deltaRho;
         (D.f[TE])[kte]     = c1o54  * deltaRho;
         (D.f[TW])[ktw]     = c1o54  * deltaRho;
         (D.f[BE])[kbe]     = c1o54  * deltaRho;
         (D.f[BS])[kbs]     = c1o54  * deltaRho;
         (D.f[TN])[ktn]     = c1o54  * deltaRho;
         (D.f[TS])[kts]     = c1o54  * deltaRho;
         (D.f[BN])[kbn]     = c1o54  * deltaRho;
         (D.f[BSW])[kbsw]   = c1o216 * deltaRho;
         (D.f[TNE])[ktne]   = c1o216 * deltaRho;
         (D.f[TSW])[ktsw]   = c1o216 * deltaRho;
         (D.f[BNE])[kbne]   = c1o216 * deltaRho;
         (D.f[BNW])[kbnw]   = c1o216 * deltaRho;
         (D.f[TSE])[ktse]   = c1o216 * deltaRho;
         (D.f[TNW])[ktnw]   = c1o216 * deltaRho;
         (D.f[BSE])[kbse]   = c1o216 * deltaRho;
         (D.f[REST])[kzero] = c8o27  * deltaRho;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceDirDepBot27(  real* rhoBC,
                                                     real* DD, 
                                                     int* k_Q, 
                                                     int numberOfBCnodes, 
                                                     real om1, 
                                                     unsigned int* neighborX,
                                                     unsigned int* neighborY,
                                                     unsigned int* neighborZ,
                                                     unsigned int size_Mat, 
                                                     bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      real rho;
      ////////////////////////////////////////////////////////////////////////////////
      rho = rhoBC[k];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E,f_W,f_N,f_S,f_T,f_NE,f_SW,f_SE,f_NW,f_TE,f_TW,f_TN,f_TS,f_ZERO,f_TNE,f_TSW,f_TSE,f_TNW;//,
            //f_B,f_BW,f_BE,f_BS,f_BN,f_BSW,f_BNE,f_BNW,f_BSE;

      f_E    = (D.f[E   ])[ke   ];
      f_W    = (D.f[W   ])[kw   ];
      f_N    = (D.f[N   ])[kn   ];
      f_S    = (D.f[S   ])[ks   ];
      f_T    = (D.f[T   ])[kt   ];
      f_NE   = (D.f[NE  ])[kne  ];
      f_SW   = (D.f[SW  ])[ksw  ];
      f_SE   = (D.f[SE  ])[kse  ];
      f_NW   = (D.f[NW  ])[knw  ];
      f_TE   = (D.f[TE  ])[kte  ];
      f_TW   = (D.f[TW  ])[ktw  ];
      f_TN   = (D.f[TN  ])[ktn  ];
      f_TS   = (D.f[TS  ])[kts  ];
      f_ZERO = (D.f[REST])[kzero];
      f_TNE  = (D.f[TNE ])[ktne ];
      f_TSW  = (D.f[TSW ])[ktsw ];
      f_TSE  = (D.f[TSE ])[ktse ];
      f_TNW  = (D.f[TNW ])[ktnw ];
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //f_B   = (four*rho- four*f_SW-     eight*f_TSW-four*f_W-   eight*f_TW- four*f_NW-     eight*f_TNW-four*f_S-   eight*f_TS-four*f_ZERO+     f_T-four*f_N-   eight*f_TN- four*f_SE-     eight*f_TSE-four*f_E-   eight*f_TE- four*f_NE-     eight*f_TNE)/nine;
      //f_BW  = ( two*rho+      f_SW-      four*f_TSW+     f_W-    four*f_TW+      f_NW-      four*f_TNW- two*f_S-    four*f_TS- two*f_ZERO-four*f_T- two*f_N-    four*f_TN- five*f_SE-      four*f_TSE-five*f_E+fourteen*f_TE- five*f_NE-      four*f_TNE)/eighteen;
      //f_BE  = ( two*rho- five*f_SW-      four*f_TSW-five*f_W+fourteen*f_TW- five*f_NW-      four*f_TNW- two*f_S-    four*f_TS- two*f_ZERO-four*f_T- two*f_N-    four*f_TN+      f_SE-      four*f_TSE+     f_E-    four*f_TE+      f_NE-      four*f_TNE)/eighteen;
      //f_BS  = ( two*rho+      f_SW-      four*f_TSW- two*f_W-    four*f_TW- five*f_NW-      four*f_TNW+     f_S-    four*f_TS- two*f_ZERO-four*f_T-five*f_N+fourteen*f_TN+      f_SE-      four*f_TSE- two*f_E-    four*f_TE- five*f_NE-      four*f_TNE)/eighteen;
      //f_BN  = ( two*rho- five*f_SW-      four*f_TSW- two*f_W-    four*f_TW+      f_NW-      four*f_TNW-five*f_S+fourteen*f_TS- two*f_ZERO-four*f_T+     f_N-    four*f_TN- five*f_SE-      four*f_TSE- two*f_E-    four*f_TE+      f_NE-      four*f_TNE)/eighteen;
      //f_BSW = ( two*rho+ four*f_SW-      four*f_TSW+     f_W-    four*f_TW-  two*f_NW-      four*f_TNW+     f_S-    four*f_TS- two*f_ZERO-four*f_T-five*f_N-    four*f_TN-  two*f_SE-      four*f_TSE-five*f_E-    four*f_TE-eight*f_NE+sixtyeight*f_TNE)/seventytwo;
      //f_BNE = ( two*rho-eight*f_SW+sixtyeight*f_TSW-five*f_W-    four*f_TW-  two*f_NW-      four*f_TNW-five*f_S-    four*f_TS- two*f_ZERO-four*f_T+     f_N-    four*f_TN-  two*f_SE-      four*f_TSE+     f_E-    four*f_TE+ four*f_NE-      four*f_TNE)/seventytwo;
      //f_BNW = ( two*rho-  two*f_SW-      four*f_TSW+     f_W-    four*f_TW+ four*f_NW-      four*f_TNW-five*f_S-    four*f_TS- two*f_ZERO-four*f_T+     f_N-    four*f_TN-eight*f_SE+sixtyeight*f_TSE-five*f_E-    four*f_TE-  two*f_NE-      four*f_TNE)/seventytwo;
      //f_BSE = ( two*rho-  two*f_SW-      four*f_TSW-five*f_W-    four*f_TW-eight*f_NW+sixtyeight*f_TNW+     f_S-    four*f_TS- two*f_ZERO-four*f_T-five*f_N-    four*f_TN+ four*f_SE-      four*f_TSE+     f_E-    four*f_TE-  two*f_NE-      four*f_TNE)/seventytwo;

      //real drho   =    f_ZERO+f_E+f_W+f_N+f_S+f_T+f_B+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW+f_BNE+f_BSW+f_BSE+f_BNW;
      //real vx1     =  (f_E -f_W +f_NE-f_SW+f_SE-f_NW+f_TE-f_BW+f_BE-f_TW+ f_TNE-f_TSW+f_TSE-f_TNW+ f_BNE-f_BSW+f_BSE-f_BNW);
      //real vx2     =  (f_N -f_S +f_NE-f_SW-f_SE+f_NW+f_TN-f_BS+f_BN-f_TS+ f_TNE-f_TSW-f_TSE+f_TNW+ f_BNE-f_BSW-f_BSE+f_BNW);
      //real vx3     =  (f_T -f_B +f_TE-f_BW-f_BE+f_TW+f_TN-f_BS-f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW- f_BNE-f_BSW-f_BSE-f_BNW);

      //real cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //(D.f[REST])[kzero] = c8over27*  (drho-cusq);
      //(D.f[E])[ke]    = c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
      //(D.f[W])[kw]    = c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
      //(D.f[N])[kn]     = c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
      //(D.f[S])[ks]    = c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
      //(D.f[T])[kt]    = c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
      //(D.f[B])[kb]    = c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
      //(D.f[NE])[kne]   = c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
      //(D.f[SW])[ksw]   = c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
      //(D.f[SE])[kse]   =  c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
      //(D.f[NW])[knw]   =  c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
      //(D.f[TE])[kte]   =  c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
      //(D.f[BW])[kbw]   =  c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
      //(D.f[BE])[kbe]   =  c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
      //(D.f[TW])[ktw]   =  c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
      //(D.f[TN])[ktn]   =  c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
      //(D.f[BS])[kbs]   =  c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
      //(D.f[BN])[kbn]   =  c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
      //(D.f[TS])[kts]   =  c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
      //(D.f[TNE])[ktne]  =  c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
      //(D.f[BSW])[kbsw]  =  c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
      //(D.f[BNE])[kbne]  =  c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
      //(D.f[TSW])[ktsw]  =  c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
      //(D.f[TSE])[ktse]  =  c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
      //(D.f[BNW])[kbnw]  =  c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
      //(D.f[BSE])[kbse]  =  c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
      //(D.f[TNW])[ktnw]  =  c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);
      real drho   =    f_ZERO+f_E+f_W+f_N+f_S+f_T+f_NE+f_SW+f_SE+f_NW+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      real dTop   =    f_T+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      (D.f[B])[kb]     = (f_T+c2o27)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c2o27;
      (D.f[BW])[kbw]   = (f_TW+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[BE])[kbe]   = (f_TE+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[BS])[kbs]   = (f_TS+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[BN])[kbn]   = (f_TN+c1o54)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o54;
      (D.f[BSW])[kbsw] = (f_TSW+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[BNE])[kbne] = (f_TNE+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[BNW])[kbnw] = (f_TNW+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
      (D.f[BSE])[kbse] = (f_TSE+c1o216)*(rho-drho+c1o1/c6o1)/(dTop+c1o1/c6o1)-c1o216;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressNoRhoDevice27(  real* rhoBC,
												 real* DD, 
												 int* k_Q, 
												 int* k_N, 
												 int numberOfBCnodes, 
												 real om1, 
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
												 bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      //unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int K1QK  = k_N[k];
      //unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==true)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f1_E    = (D.f[E   ])[k1e   ];
      real f1_W    = (D.f[W   ])[k1w   ];
      real f1_N    = (D.f[N   ])[k1n   ];
      real f1_S    = (D.f[S   ])[k1s   ];
      real f1_T    = (D.f[T   ])[k1t   ];
      real f1_B    = (D.f[B   ])[k1b   ];
      real f1_NE   = (D.f[NE  ])[k1ne  ];
      real f1_SW   = (D.f[SW  ])[k1sw  ];
      real f1_SE   = (D.f[SE  ])[k1se  ];
      real f1_NW   = (D.f[NW  ])[k1nw  ];
      real f1_TE   = (D.f[TE  ])[k1te  ];
      real f1_BW   = (D.f[BW  ])[k1bw  ];
      real f1_BE   = (D.f[BE  ])[k1be  ];
      real f1_TW   = (D.f[TW  ])[k1tw  ];
      real f1_TN   = (D.f[TN  ])[k1tn  ];
      real f1_BS   = (D.f[BS  ])[k1bs  ];
      real f1_BN   = (D.f[BN  ])[k1bn  ];
      real f1_TS   = (D.f[TS  ])[k1ts  ];
      //real f1_ZERO = (D.f[REST])[k1zero];
      real f1_TNE  = (D.f[TNE ])[k1tne ];
      real f1_TSW  = (D.f[TSW ])[k1tsw ];
      real f1_TSE  = (D.f[TSE ])[k1tse ];
      real f1_TNW  = (D.f[TNW ])[k1tnw ];
      real f1_BNE  = (D.f[BNE ])[k1bne ];
      real f1_BSW  = (D.f[BSW ])[k1bsw ];
      real f1_BSE  = (D.f[BSE ])[k1bse ];
      real f1_BNW  = (D.f[BNW ])[k1bnw ];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real f_E    = (D.f[E   ])[ke   ];
      real f_W    = (D.f[W   ])[kw   ];
      real f_N    = (D.f[N   ])[kn   ];
      real f_S    = (D.f[S   ])[ks   ];
      real f_T    = (D.f[T   ])[kt   ];
      real f_B    = (D.f[B   ])[kb   ];
      real f_NE   = (D.f[NE  ])[kne  ];
      real f_SW   = (D.f[SW  ])[ksw  ];
      real f_SE   = (D.f[SE  ])[kse  ];
      real f_NW   = (D.f[NW  ])[knw  ];
      real f_TE   = (D.f[TE  ])[kte  ];
      real f_BW   = (D.f[BW  ])[kbw  ];
      real f_BE   = (D.f[BE  ])[kbe  ];
      real f_TW   = (D.f[TW  ])[ktw  ];
      real f_TN   = (D.f[TN  ])[ktn  ];
      real f_BS   = (D.f[BS  ])[kbs  ];
      real f_BN   = (D.f[BN  ])[kbn  ];
      real f_TS   = (D.f[TS  ])[kts  ];
      //real f_ZERO = (D.f[REST])[kzero];
      real f_TNE  = (D.f[TNE ])[ktne ];
      real f_TSW  = (D.f[TSW ])[ktsw ];
      real f_TSE  = (D.f[TSE ])[ktse ];
      real f_TNW  = (D.f[TNW ])[ktnw ];
      real f_BNE  = (D.f[BNE ])[kbne ];
      real f_BSW  = (D.f[BSW ])[kbsw ];
      real f_BSE  = (D.f[BSE ])[kbse ];
      real f_BNW  = (D.f[BNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////

      //real vx1, vx2, vx3, drho;
      //real vx1, vx2, vx3, drho, drho1;
      //////////////////////////////////////////////////////////////////////////
	  //Dichte
    //   drho1  =  f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
    //             f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW + 
    //             f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((D.f[REST])[k1zero]); 
    //   drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
    //             f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
    //             f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 
      
      //////////////////////////////////////////////////////////////////////////
	  //Ux

	  //vx1    =  (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
   //               ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
   //               (f_E - f_W)) /(one + drho); 


   //   vx2    =   ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
   //               ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
   //               (f_N - f_S)) /(one + drho); 

   //   vx3    =   (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
   //               (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
   //               (f_T - f_B)) /(one + drho); 


      //real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

   //   //////////////////////////////////////////////////////////////////////////
	  ////real omega = om1;
   //   real cusq  = c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
   //   //////////////////////////////////////////////////////////////////////////
	  ////T�st MK
	  ////if(vx1 < zero) vx1 = zero;
   //   //////////////////////////////////////////////////////////////////////////
   //   real fZERO = c8over27*  (drho1-(one + drho1)*(cusq))                                                           ;
   //   real fE    = c2over27*  (drho1+(one + drho1)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
   //   real fW    = c2over27*  (drho1+(one + drho1)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
   //   real fN    = c2over27*  (drho1+(one + drho1)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
   //   real fS    = c2over27*  (drho1+(one + drho1)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
   //   real fT    = c2over27*  (drho1+(one + drho1)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
   //   real fB    = c2over27*  (drho1+(one + drho1)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
   //   real fNE   = c1over54*  (drho1+(one + drho1)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
   //   real fSW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
   //   real fSE   = c1over54*  (drho1+(one + drho1)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
   //   real fNW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
   //   real fTE   = c1over54*  (drho1+(one + drho1)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
   //   real fBW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
   //   real fBE   = c1over54*  (drho1+(one + drho1)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
   //   real fTW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
   //   real fTN   = c1over54*  (drho1+(one + drho1)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
   //   real fBS   = c1over54*  (drho1+(one + drho1)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
   //   real fBN   = c1over54*  (drho1+(one + drho1)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
   //   real fTS   = c1over54*  (drho1+(one + drho1)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
   //   real fTNE  = c1over216* (drho1+(one + drho1)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
   //   real fBSW  = c1over216* (drho1+(one + drho1)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
   //   real fBNE  = c1over216* (drho1+(one + drho1)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
   //   real fTSW  = c1over216* (drho1+(one + drho1)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
   //   real fTSE  = c1over216* (drho1+(one + drho1)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
   //   real fBNW  = c1over216* (drho1+(one + drho1)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
   //   real fBSE  = c1over216* (drho1+(one + drho1)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
   //   real fTNW  = c1over216* (drho1+(one + drho1)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

	  real cs = c1o1 / sqrtf(c3o1);
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //no velocity
	  //////////////////////////////////////////
      f_E    = f1_E   * cs + (c1o1 - cs) * f_E   ;
      f_W    = f1_W   * cs + (c1o1 - cs) * f_W   ;
      f_N    = f1_N   * cs + (c1o1 - cs) * f_N   ;
      f_S    = f1_S   * cs + (c1o1 - cs) * f_S   ;
      f_T    = f1_T   * cs + (c1o1 - cs) * f_T   ;
      f_B    = f1_B   * cs + (c1o1 - cs) * f_B   ;
      f_NE   = f1_NE  * cs + (c1o1 - cs) * f_NE  ;
      f_SW   = f1_SW  * cs + (c1o1 - cs) * f_SW  ;
      f_SE   = f1_SE  * cs + (c1o1 - cs) * f_SE  ;
      f_NW   = f1_NW  * cs + (c1o1 - cs) * f_NW  ;
      f_TE   = f1_TE  * cs + (c1o1 - cs) * f_TE  ;
      f_BW   = f1_BW  * cs + (c1o1 - cs) * f_BW  ;
      f_BE   = f1_BE  * cs + (c1o1 - cs) * f_BE  ;
      f_TW   = f1_TW  * cs + (c1o1 - cs) * f_TW  ;
      f_TN   = f1_TN  * cs + (c1o1 - cs) * f_TN  ;
      f_BS   = f1_BS  * cs + (c1o1 - cs) * f_BS  ;
      f_BN   = f1_BN  * cs + (c1o1 - cs) * f_BN  ;
      f_TS   = f1_TS  * cs + (c1o1 - cs) * f_TS  ;
      f_TNE  = f1_TNE * cs + (c1o1 - cs) * f_TNE ;
      f_TSW  = f1_TSW * cs + (c1o1 - cs) * f_TSW ;
      f_TSE  = f1_TSE * cs + (c1o1 - cs) * f_TSE ;
      f_TNW  = f1_TNW * cs + (c1o1 - cs) * f_TNW ;
      f_BNE  = f1_BNE * cs + (c1o1 - cs) * f_BNE ;
      f_BSW  = f1_BSW * cs + (c1o1 - cs) * f_BSW ;
      f_BSE  = f1_BSE * cs + (c1o1 - cs) * f_BSE ;
      f_BNW  = f1_BNW * cs + (c1o1 - cs) * f_BNW ;
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //with velocity
	  //if(true){//vx1 >= zero){
		 // real csMvx = one / sqrtf(three) - vx1;
		 // //real csMvy = one / sqrtf(three) - vx2;
		 // ///////////////////////////////////////////
		 // // X
		 // f_W   = f1_W   * csMvx + (one - csMvx) * f_W   ;//- c2over27  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_NW  = f1_NW  * csMvx + (one - csMvx) * f_NW  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_SW  = f1_SW  * csMvx + (one - csMvx) * f_SW  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_TW  = f1_TW  * csMvx + (one - csMvx) * f_TW  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_BW  = f1_BW  * csMvx + (one - csMvx) * f_BW  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_TNW = f1_TNW * csMvx + (one - csMvx) * f_TNW ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_TSW = f1_TSW * csMvx + (one - csMvx) * f_TSW ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_BNW = f1_BNW * csMvx + (one - csMvx) * f_BNW ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // f_BSW = f1_BSW * csMvx + (one - csMvx) * f_BSW ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx1);
		 // ///////////////////////////////////////////
		 // // Y
		 // //f_S   = f1_S   * csMvy + (one - csMvy) * f_S   ;//- c2over27  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_SE  = f1_SE  * csMvy + (one - csMvy) * f_SE  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_SW  = f1_SW  * csMvy + (one - csMvy) * f_SW  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_TS  = f1_TS  * csMvy + (one - csMvy) * f_TS  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_BS  = f1_BS  * csMvy + (one - csMvy) * f_BS  ;//- c1over54  * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_TSE = f1_TSE * csMvy + (one - csMvy) * f_TSE ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_TSW = f1_TSW * csMvy + (one - csMvy) * f_TSW ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_BSE = f1_BSE * csMvy + (one - csMvy) * f_BSE ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_BSW = f1_BSW * csMvy + (one - csMvy) * f_BSW ;//- c1over216 * ((drho + drho1)*c1o2-((drho + drho1)*c1o2 )*three*vx2);
		 // //f_S   = f1_S   * csMvy + (one - csMvy) * f_S;
		 // //f_SE  = f1_SE  * csMvy + (one - csMvy) * f_SE;
		 // //f_SW  = f1_SW  * csMvy + (one - csMvy) * f_SW;
		 // //f_TS  = f1_TS  * csMvy + (one - csMvy) * f_TS;
		 // //f_BS  = f1_BS  * csMvy + (one - csMvy) * f_BS;
		 // //f_TSE = f1_TSE * csMvy + (one - csMvy) * f_TSE;
		 // //f_TSW = f1_TSW * csMvy + (one - csMvy) * f_TSW;
		 // //f_BSE = f1_BSE * csMvy + (one - csMvy) * f_BSE;
		 // //f_BSW = f1_BSW * csMvy + (one - csMvy) * f_BSW;
		 // //////////////////////////////////////////////////////////////////////////
	  //}
	  //else
	  //{
		 // ///////////////////////////////////////////
		 // // X
		 // vx1   = vx1 * 0.9;
		 // f_W   = f_E   - six * c2over27  * ( vx1        );
		 // f_NW  = f_SE  - six * c1over54  * ( vx1-vx2    );
		 // f_SW  = f_NE  - six * c1over54  * ( vx1+vx2    );
		 // f_TW  = f_BE  - six * c1over54  * ( vx1    -vx3);
		 // f_BW  = f_TE  - six * c1over54  * ( vx1    +vx3);
		 // f_TNW = f_BSE - six * c1over216 * ( vx1-vx2-vx3);
		 // f_TSW = f_BNE - six * c1over216 * ( vx1+vx2-vx3);
		 // f_BNW = f_TSE - six * c1over216 * ( vx1-vx2+vx3);
		 // f_BSW = f_TNE - six * c1over216 * ( vx1+vx2+vx3);
		 // ///////////////////////////////////////////
		 // // Y
		 // //vx2   = vx2 * 0.9;
		 // //f_S   = f_N   - six * c2over27  * (     vx2    );
		 // //f_SE  = f_NW  - six * c1over54  * (-vx1+vx2    );
		 // //f_SW  = f_NE  - six * c1over54  * ( vx1+vx2    );
		 // //f_TS  = f_BN  - six * c1over54  * (     vx2-vx3);
		 // //f_BS  = f_TN  - six * c1over54  * (     vx2+vx3);
		 // //f_TSE = f_BNW - six * c1over216 * (-vx1+vx2-vx3);
		 // //f_TSW = f_BNE - six * c1over216 * ( vx1+vx2-vx3);
		 // //f_BSE = f_TNW - six * c1over216 * (-vx1+vx2+vx3);
		 // //f_BSW = f_TNE - six * c1over216 * ( vx1+vx2+vx3);
		 // ///////////////////////////////////////////
	  //}
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //////////////////////////////////////////////////////////////////////////
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////
      //__syncthreads();
	  // -X
	  //(D.f[E   ])[ke   ] = f_E   ;
	  //(D.f[SE  ])[kse  ] = f_SE  ;
	  //(D.f[NE  ])[kne  ] = f_NE  ;
	  //(D.f[BE  ])[kbe  ] = f_BE  ;
	  //(D.f[TE  ])[kte  ] = f_TE  ;
	  //(D.f[TSE ])[ktse ] = f_TSE ;
	  //(D.f[TNE ])[ktne ] = f_TNE ;
	  //(D.f[BSE ])[kbse ] = f_BSE ;
	  //(D.f[BNE ])[kbne ] = f_BNE ;     
	  // X
	  (D.f[W   ])[kw   ] = f_W   ;
	  (D.f[SW  ])[ksw  ] = f_SW  ;
	  (D.f[NW  ])[knw  ] = f_NW  ;
	  (D.f[BW  ])[kbw  ] = f_BW  ;
	  (D.f[TW  ])[ktw  ] = f_TW  ;
	  (D.f[TSW ])[ktsw ] = f_TSW ;
	  (D.f[TNW ])[ktnw ] = f_TNW ;
	  (D.f[BSW ])[kbsw ] = f_BSW ;
	  (D.f[BNW ])[kbnw ] = f_BNW ;     
	  // Y
	  //(D.f[S   ])[ks   ] = f_S   ;
	  //(D.f[SE  ])[kse  ] = f_SE  ;
	  //(D.f[SW  ])[ksw  ] = f_SW  ;
	  //(D.f[TS  ])[kts  ] = f_TS  ;
	  //(D.f[BS  ])[kbs  ] = f_BS  ;
	  //(D.f[TSE ])[ktse ] = f_TSE ;
	  //(D.f[TSW ])[ktsw ] = f_TSW ;
	  //(D.f[BSE ])[kbse ] = f_BSE ;
	  //(D.f[BSW ])[kbsw ] = f_BSW ;     
	  // Z
	  //(D.f[B   ])[kb   ] = f_B   ;
	  //(D.f[BE  ])[kbe  ] = f_BE  ;
	  //(D.f[BW  ])[kbw  ] = f_BW  ;
	  //(D.f[BN  ])[kbn  ] = f_BN  ;
	  //(D.f[BS  ])[kbs  ] = f_BS  ;
	  //(D.f[BNE ])[kbne ] = f_BNE ;
	  //(D.f[BNW ])[kbnw ] = f_BNW ;
	  //(D.f[BSE ])[kbse ] = f_BSE ;
	  //(D.f[BSW ])[kbsw ] = f_BSW ;     
      //////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceOld27(real* rhoBC,
                                             real* DD, 
                                             int* k_Q, 
                                             int* k_N, 
                                             int numberOfBCnodes, 
                                             real om1, 
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int size_Mat, 
                                             bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int K1QK  = k_N[k];
      unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[E   ])[k1e   ];
      f1_E    = (D.f[W   ])[k1w   ];
      f1_S    = (D.f[N   ])[k1n   ];
      f1_N    = (D.f[S   ])[k1s   ];
      f1_B    = (D.f[T   ])[k1t   ];
      f1_T    = (D.f[B   ])[k1b   ];
      f1_SW   = (D.f[NE  ])[k1ne  ];
      f1_NE   = (D.f[SW  ])[k1sw  ];
      f1_NW   = (D.f[SE  ])[k1se  ];
      f1_SE   = (D.f[NW  ])[k1nw  ];
      f1_BW   = (D.f[TE  ])[k1te  ];
      f1_TE   = (D.f[BW  ])[k1bw  ];
      f1_TW   = (D.f[BE  ])[k1be  ];
      f1_BE   = (D.f[TW  ])[k1tw  ];
      f1_BS   = (D.f[TN  ])[k1tn  ];
      f1_TN   = (D.f[BS  ])[k1bs  ];
      f1_TS   = (D.f[BN  ])[k1bn  ];
      f1_BN   = (D.f[TS  ])[k1ts  ];
      f1_ZERO = (D.f[REST])[k1zero];
      f1_BSW  = (D.f[TNE ])[k1tne ];
      f1_BNE  = (D.f[TSW ])[k1tsw ];
      f1_BNW  = (D.f[TSE ])[k1tse ];
      f1_BSE  = (D.f[TNW ])[k1tnw ];
      f1_TSW  = (D.f[BNE ])[k1bne ];
      f1_TNE  = (D.f[BSW ])[k1bsw ];
      f1_TNW  = (D.f[BSE ])[k1bse ];
      f1_TSE  = (D.f[BNW ])[k1bnw ];

      //////////////////////////////////////////////////////////////////////////
      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                          f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

	  //drho1 = (drho1 + rhoBC[k])/2.f;
	  drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[E   ])[ke   ] = f1_W   -c2o27*drho1;   //  c1o100;  // zero;  //
      (D.f[W   ])[kw   ] = f1_E   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[N   ])[kn   ] = f1_S   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[S   ])[ks   ] = f1_N   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[T   ])[kt   ] = f1_B   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[B   ])[kb   ] = f1_T   -c2o27*drho1;	//  c1o100;  // zero;  //
      (D.f[NE  ])[kne  ] = f1_SW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[SW  ])[ksw  ] = f1_NE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[SE  ])[kse  ] = f1_NW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[NW  ])[knw  ] = f1_SE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TE  ])[kte  ] = f1_BW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BW  ])[kbw  ] = f1_TE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BE  ])[kbe  ] = f1_TW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TW  ])[ktw  ] = f1_BE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TN  ])[ktn  ] = f1_BS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BS  ])[kbs  ] = f1_TN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BN  ])[kbn  ] = f1_TS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TS  ])[kts  ] = f1_BN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[REST])[kzero] = f1_ZERO-c8o27*drho1;	//  c1o100;  // zero;  //
      (D.f[TNE ])[ktne ] = f1_BSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[TSW ])[ktsw ] = f1_BNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[TSE ])[ktse ] = f1_BNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[TNW ])[ktnw ] = f1_BSE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BNE ])[kbne ] = f1_TSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BSW ])[kbsw ] = f1_TNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BSE ])[kbse ] = f1_TNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BNW ])[kbnw ] = f1_TSE -c1o216*drho1;  //  c1o100;  // zero;  //      
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceEQZ27(real* rhoBC,
                                             real* DD, 
                                             int* k_Q, 
                                             int* k_N,
											 real* kTestRE,
                                             int numberOfBCnodes, 
                                             real om1, 
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int size_Mat, 
                                             bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int K1QK  = k_N[k];
      unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==true)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////
    //   Distributions27 kDistTest;
    //      kDistTest.f[E   ] = &kTestRE[E   *numberOfBCnodes];
    //      kDistTest.f[W   ] = &kTestRE[W   *numberOfBCnodes];
    //      kDistTest.f[N   ] = &kTestRE[N   *numberOfBCnodes];
    //      kDistTest.f[S   ] = &kTestRE[S   *numberOfBCnodes];
    //      kDistTest.f[T   ] = &kTestRE[T   *numberOfBCnodes];
    //      kDistTest.f[B   ] = &kTestRE[B   *numberOfBCnodes];
    //      kDistTest.f[NE  ] = &kTestRE[NE  *numberOfBCnodes];
    //      kDistTest.f[SW  ] = &kTestRE[SW  *numberOfBCnodes];
    //      kDistTest.f[SE  ] = &kTestRE[SE  *numberOfBCnodes];
    //      kDistTest.f[NW  ] = &kTestRE[NW  *numberOfBCnodes];
    //      kDistTest.f[TE  ] = &kTestRE[TE  *numberOfBCnodes];
    //      kDistTest.f[BW  ] = &kTestRE[BW  *numberOfBCnodes];
    //      kDistTest.f[BE  ] = &kTestRE[BE  *numberOfBCnodes];
    //      kDistTest.f[TW  ] = &kTestRE[TW  *numberOfBCnodes];
    //      kDistTest.f[TN  ] = &kTestRE[TN  *numberOfBCnodes];
    //      kDistTest.f[BS  ] = &kTestRE[BS  *numberOfBCnodes];
    //      kDistTest.f[BN  ] = &kTestRE[BN  *numberOfBCnodes];
    //      kDistTest.f[TS  ] = &kTestRE[TS  *numberOfBCnodes];
    //      kDistTest.f[REST] = &kTestRE[REST*numberOfBCnodes];
    //      kDistTest.f[TNE ] = &kTestRE[TNE *numberOfBCnodes];
    //      kDistTest.f[TSW ] = &kTestRE[TSW *numberOfBCnodes];
    //      kDistTest.f[TSE ] = &kTestRE[TSE *numberOfBCnodes];
    //      kDistTest.f[TNW ] = &kTestRE[TNW *numberOfBCnodes];
    //      kDistTest.f[BNE ] = &kTestRE[BNE *numberOfBCnodes];
    //      kDistTest.f[BSW ] = &kTestRE[BSW *numberOfBCnodes];
    //      kDistTest.f[BSE ] = &kTestRE[BSE *numberOfBCnodes];
    //      kDistTest.f[BNW ] = &kTestRE[BNW *numberOfBCnodes];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //real f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   //f1_W    = (D.f[E   ])[k1e   ];
   //   //f1_E    = (D.f[W   ])[k1w   ];
   //   //f1_S    = (D.f[N   ])[k1n   ];
   //   //f1_N    = (D.f[S   ])[k1s   ];
   //   //f1_B    = (D.f[T   ])[k1t   ];
   //   //f1_T    = (D.f[B   ])[k1b   ];
   //   //f1_SW   = (D.f[NE  ])[k1ne  ];
   //   //f1_NE   = (D.f[SW  ])[k1sw  ];
   //   //f1_NW   = (D.f[SE  ])[k1se  ];
   //   //f1_SE   = (D.f[NW  ])[k1nw  ];
   //   //f1_BW   = (D.f[TE  ])[k1te  ];
   //   //f1_TE   = (D.f[BW  ])[k1bw  ];
   //   //f1_TW   = (D.f[BE  ])[k1be  ];
   //   //f1_BE   = (D.f[TW  ])[k1tw  ];
   //   //f1_BS   = (D.f[TN  ])[k1tn  ];
   //   //f1_TN   = (D.f[BS  ])[k1bs  ];
   //   //f1_TS   = (D.f[BN  ])[k1bn  ];
   //   //f1_BN   = (D.f[TS  ])[k1ts  ];
   //   //f1_ZERO = (D.f[REST])[k1zero];
   //   //f1_BSW  = (D.f[TNE ])[k1tne ];
   //   //f1_BNE  = (D.f[TSW ])[k1tsw ];
   //   //f1_BNW  = (D.f[TSE ])[k1tse ];
   //   //f1_BSE  = (D.f[TNW ])[k1tnw ];
   //   //f1_TSW  = (D.f[BNE ])[k1bne ];
   //   //f1_TNE  = (D.f[BSW ])[k1bsw ];
   //   //f1_TNW  = (D.f[BSE ])[k1bse ];
   //   //f1_TSE  = (D.f[BNW ])[k1bnw ];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   real f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   f1_E    = (D.f[E   ])[k1e   ];
   //   f1_W    = (D.f[W   ])[k1w   ];
   //   f1_N    = (D.f[N   ])[k1n   ];
   //   f1_S    = (D.f[S   ])[k1s   ];
   //   f1_T    = (D.f[T   ])[k1t   ];
   //   f1_B    = (D.f[B   ])[k1b   ];
   //   f1_NE   = (D.f[NE  ])[k1ne  ];
   //   f1_SW   = (D.f[SW  ])[k1sw  ];
   //   f1_SE   = (D.f[SE  ])[k1se  ];
   //   f1_NW   = (D.f[NW  ])[k1nw  ];
   //   f1_TE   = (D.f[TE  ])[k1te  ];
   //   f1_BW   = (D.f[BW  ])[k1bw  ];
   //   f1_BE   = (D.f[BE  ])[k1be  ];
   //   f1_TW   = (D.f[TW  ])[k1tw  ];
   //   f1_TN   = (D.f[TN  ])[k1tn  ];
   //   f1_BS   = (D.f[BS  ])[k1bs  ];
   //   f1_BN   = (D.f[BN  ])[k1bn  ];
   //   f1_TS   = (D.f[TS  ])[k1ts  ];
   //   f1_ZERO = (D.f[REST])[k1zero];
   //   f1_TNE  = (D.f[TNE ])[k1tne ];
   //   f1_TSW  = (D.f[TSW ])[k1tsw ];
   //   f1_TSE  = (D.f[TSE ])[k1tse ];
   //   f1_TNW  = (D.f[TNW ])[k1tnw ];
   //   f1_BNE  = (D.f[BNE ])[k1bne ];
   //   f1_BSW  = (D.f[BSW ])[k1bsw ];
   //   f1_BSE  = (D.f[BSE ])[k1bse ];
   //   f1_BNW  = (D.f[BNW ])[k1bnw ];
   //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   //////////////////////////////////////////////////////////////////////////
   //   real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+ f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;
	  //real vx1      = (((f1_TNE-f1_BSW)+(f1_BSE-f1_TNW)+(f1_BNE-f1_TSW)+(f1_TSE-f1_BNW)) + (((f1_NE-f1_SW)+(f1_TE-f1_BW))+((f1_SE-f1_NW)+(f1_BE-f1_TW))) + (f1_E-f1_W)) / (one + drho1);
	  //real vx2      = (((f1_TNE-f1_BSW)+(f1_TNW-f1_BSE)+(f1_BNE-f1_TSW)+(f1_BNW-f1_TSE)) + (((f1_NE-f1_SW)+(f1_TN-f1_BS))+((f1_BN-f1_TS)+(f1_NW-f1_SE))) + (f1_N-f1_S)) / (one + drho1);
	  //real vx3      = (((f1_TNE-f1_BSW)+(f1_TNW-f1_BSE)+(f1_TSW-f1_BNE)+(f1_TSE-f1_BNW)) + (((f1_TE-f1_BW)+(f1_TN-f1_BS))+((f1_TW-f1_BE)+(f1_TS-f1_BN))) + (f1_T-f1_B)) / (one + drho1);
   //   //////////////////////////////////////////////////////////////////////////
	  ////real omega = om1;
   //   real cusq  = c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
   //   //////////////////////////////////////////////////////////////////////////
	  ////T�st MK
	  ////if(vx1 < zero) vx1 = zero;
   //   //////////////////////////////////////////////////////////////////////////
	  ////becomes higher with neighbor source and lower with local source
   //   //real fZERO = c8over27*  (rhoBC[k]-(one + rhoBC[k])*(cusq))                                                           ;
   //   //real fE    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
   //   //real fW    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
   //   //real fN    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
   //   //real fS    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
   //   //real fT    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
   //   //real fB    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
   //   //real fNE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
   //   //real fSW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
   //   //real fSE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
   //   //real fNW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
   //   //real fTE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
   //   //real fBW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
   //   //real fBE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
   //   //real fTW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
   //   //real fTN   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
   //   //real fBS   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
   //   //real fBN   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
   //   //real fTS   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
   //   //real fTNE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
   //   //real fBSW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
   //   //real fBNE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
   //   //real fTSW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
   //   //real fTSE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
   //   //real fBNW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
   //   //real fBSE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
   //   //real fTNW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
   //   //////////////////////////////////////////////////////////////////////////
	  //// based on VirtualFluids (kucher + fard)
   //   real fZERO = c8over27  * rhoBC[k] * (one                                                                      - cusq);
   //   real fE    = c2over27  * rhoBC[k] * (one + three * ( vx1        ) + c9over2 * ( vx1        ) * ( vx1        ) - cusq);
   //   real fW    = c2over27  * rhoBC[k] * (one + three * (-vx1        ) + c9over2 * (-vx1        ) * (-vx1        ) - cusq);
   //   real fN    = c2over27  * rhoBC[k] * (one + three * (     vx2    ) + c9over2 * (     vx2    ) * (     vx2    ) - cusq);
   //   real fS    = c2over27  * rhoBC[k] * (one + three * (    -vx2    ) + c9over2 * (    -vx2    ) * (    -vx2    ) - cusq);
   //   real fT    = c2over27  * rhoBC[k] * (one + three * (         vx3) + c9over2 * (         vx3) * (         vx3) - cusq);
   //   real fB    = c2over27  * rhoBC[k] * (one + three * (        -vx3) + c9over2 * (        -vx3) * (        -vx3) - cusq);
   //   real fNE   = c1over54  * rhoBC[k] * (one + three * ( vx1+vx2    ) + c9over2 * ( vx1+vx2    ) * ( vx1+vx2    ) - cusq);
   //   real fSW   = c1over54  * rhoBC[k] * (one + three * (-vx1-vx2    ) + c9over2 * (-vx1-vx2    ) * (-vx1-vx2    ) - cusq);
   //   real fSE   = c1over54  * rhoBC[k] * (one + three * ( vx1-vx2    ) + c9over2 * ( vx1-vx2    ) * ( vx1-vx2    ) - cusq);
   //   real fNW   = c1over54  * rhoBC[k] * (one + three * (-vx1+vx2    ) + c9over2 * (-vx1+vx2    ) * (-vx1+vx2    ) - cusq);
   //   real fTE   = c1over54  * rhoBC[k] * (one + three * ( vx1    +vx3) + c9over2 * ( vx1    +vx3) * ( vx1    +vx3) - cusq);
   //   real fBW   = c1over54  * rhoBC[k] * (one + three * (-vx1    -vx3) + c9over2 * (-vx1    -vx3) * (-vx1    -vx3) - cusq);
   //   real fBE   = c1over54  * rhoBC[k] * (one + three * ( vx1    -vx3) + c9over2 * ( vx1    -vx3) * ( vx1    -vx3) - cusq);
   //   real fTW   = c1over54  * rhoBC[k] * (one + three * (-vx1    +vx3) + c9over2 * (-vx1    +vx3) * (-vx1    +vx3) - cusq);
   //   real fTN   = c1over54  * rhoBC[k] * (one + three * (     vx2+vx3) + c9over2 * (     vx2+vx3) * (     vx2+vx3) - cusq);
   //   real fBS   = c1over54  * rhoBC[k] * (one + three * (    -vx2-vx3) + c9over2 * (    -vx2-vx3) * (    -vx2-vx3) - cusq);
   //   real fBN   = c1over54  * rhoBC[k] * (one + three * (     vx2-vx3) + c9over2 * (     vx2-vx3) * (     vx2-vx3) - cusq);
   //   real fTS   = c1over54  * rhoBC[k] * (one + three * (    -vx2+vx3) + c9over2 * (    -vx2+vx3) * (    -vx2+vx3) - cusq);
   //   real fTNE  = c1over216 * rhoBC[k] * (one + three * ( vx1+vx2+vx3) + c9over2 * ( vx1+vx2+vx3) * ( vx1+vx2+vx3) - cusq);
   //   real fBSW  = c1over216 * rhoBC[k] * (one + three * (-vx1-vx2-vx3) + c9over2 * (-vx1-vx2-vx3) * (-vx1-vx2-vx3) - cusq);
   //   real fBNE  = c1over216 * rhoBC[k] * (one + three * ( vx1+vx2-vx3) + c9over2 * ( vx1+vx2-vx3) * ( vx1+vx2-vx3) - cusq);
   //   real fTSW  = c1over216 * rhoBC[k] * (one + three * (-vx1-vx2+vx3) + c9over2 * (-vx1-vx2+vx3) * (-vx1-vx2+vx3) - cusq);
   //   real fTSE  = c1over216 * rhoBC[k] * (one + three * ( vx1-vx2+vx3) + c9over2 * ( vx1-vx2+vx3) * ( vx1-vx2+vx3) - cusq);
   //   real fBNW  = c1over216 * rhoBC[k] * (one + three * (-vx1+vx2-vx3) + c9over2 * (-vx1+vx2-vx3) * (-vx1+vx2-vx3) - cusq);
   //   real fBSE  = c1over216 * rhoBC[k] * (one + three * ( vx1-vx2-vx3) + c9over2 * ( vx1-vx2-vx3) * ( vx1-vx2-vx3) - cusq);
   //   real fTNW  = c1over216 * rhoBC[k] * (one + three * (-vx1+vx2+vx3) + c9over2 * (-vx1+vx2+vx3) * (-vx1+vx2+vx3) - cusq);
   ////   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////test
   ////   real fZERO = c8over27  * ((drho1 + rhoBC[k]) / two) * (one                                                                      - cusq);
   ////   real fE    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1        ) + c9over2 * ( vx1        ) * ( vx1        ) - cusq);
   ////   real fW    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1        ) + c9over2 * (-vx1        ) * (-vx1        ) - cusq);
   ////   real fN    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2    ) + c9over2 * (     vx2    ) * (     vx2    ) - cusq);
   ////   real fS    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2    ) + c9over2 * (    -vx2    ) * (    -vx2    ) - cusq);
   ////   real fT    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (         vx3) + c9over2 * (         vx3) * (         vx3) - cusq);
   ////   real fB    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (        -vx3) + c9over2 * (        -vx3) * (        -vx3) - cusq);
   ////   real fNE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2    ) + c9over2 * ( vx1+vx2    ) * ( vx1+vx2    ) - cusq);
   ////   real fSW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2    ) + c9over2 * (-vx1-vx2    ) * (-vx1-vx2    ) - cusq);
   ////   real fSE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2    ) + c9over2 * ( vx1-vx2    ) * ( vx1-vx2    ) - cusq);
   ////   real fNW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2    ) + c9over2 * (-vx1+vx2    ) * (-vx1+vx2    ) - cusq);
   ////   real fTE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1    +vx3) + c9over2 * ( vx1    +vx3) * ( vx1    +vx3) - cusq);
   ////   real fBW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1    -vx3) + c9over2 * (-vx1    -vx3) * (-vx1    -vx3) - cusq);
   ////   real fBE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1    -vx3) + c9over2 * ( vx1    -vx3) * ( vx1    -vx3) - cusq);
   ////   real fTW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1    +vx3) + c9over2 * (-vx1    +vx3) * (-vx1    +vx3) - cusq);
   ////   real fTN   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2+vx3) + c9over2 * (     vx2+vx3) * (     vx2+vx3) - cusq);
   ////   real fBS   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2-vx3) + c9over2 * (    -vx2-vx3) * (    -vx2-vx3) - cusq);
   ////   real fBN   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2-vx3) + c9over2 * (     vx2-vx3) * (     vx2-vx3) - cusq);
   ////   real fTS   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2+vx3) + c9over2 * (    -vx2+vx3) * (    -vx2+vx3) - cusq);
   ////   real fTNE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2+vx3) + c9over2 * ( vx1+vx2+vx3) * ( vx1+vx2+vx3) - cusq);
   ////   real fBSW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2-vx3) + c9over2 * (-vx1-vx2-vx3) * (-vx1-vx2-vx3) - cusq);
   ////   real fBNE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2-vx3) + c9over2 * ( vx1+vx2-vx3) * ( vx1+vx2-vx3) - cusq);
   ////   real fTSW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2+vx3) + c9over2 * (-vx1-vx2+vx3) * (-vx1-vx2+vx3) - cusq);
   ////   real fTSE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2+vx3) + c9over2 * ( vx1-vx2+vx3) * ( vx1-vx2+vx3) - cusq);
   ////   real fBNW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2-vx3) + c9over2 * (-vx1+vx2-vx3) * (-vx1+vx2-vx3) - cusq);
   ////   real fBSE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2-vx3) + c9over2 * ( vx1-vx2-vx3) * ( vx1-vx2-vx3) - cusq);
   ////   real fTNW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2+vx3) + c9over2 * (-vx1+vx2+vx3) * (-vx1+vx2+vx3) - cusq);

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // based on BGK Plus Comp
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//double mfabb = (D.f[E   ])[k1e   ];
			//double mfcbb = (D.f[W   ])[k1w   ];
			//double mfbab = (D.f[N   ])[k1n   ];
			//double mfbcb = (D.f[S   ])[k1s   ];
			//double mfbba = (D.f[T   ])[k1t   ];
			//double mfbbc = (D.f[B   ])[k1b   ];
			//double mfaab = (D.f[NE  ])[k1ne  ];
			//double mfccb = (D.f[SW  ])[k1sw  ];
			//double mfacb = (D.f[SE  ])[k1se  ];
			//double mfcab = (D.f[NW  ])[k1nw  ];
			//double mfaba = (D.f[TE  ])[k1te  ];
			//double mfcbc = (D.f[BW  ])[k1bw  ];
			//double mfabc = (D.f[BE  ])[k1be  ];
			//double mfcba = (D.f[TW  ])[k1tw  ];
			//double mfbaa = (D.f[TN  ])[k1tn  ];
			//double mfbcc = (D.f[BS  ])[k1bs  ];
			//double mfbac = (D.f[BN  ])[k1bn  ];
			//double mfbca = (D.f[TS  ])[k1ts  ];
			//double mfbbb = (D.f[REST])[k1zero];
			//double mfaaa = (D.f[TNE ])[k1tne ];
			//double mfcca = (D.f[TSW ])[k1tsw ];
			//double mfaca = (D.f[TSE ])[k1tse ];
			//double mfcaa = (D.f[TNW ])[k1tnw ];
			//double mfaac = (D.f[BNE ])[k1bne ];
			//double mfccc = (D.f[BSW ])[k1bsw ];
			//double mfacc = (D.f[BSE ])[k1bse ];
			//double mfcac = (D.f[BNW ])[k1bnw ];
			real mfabb = (D.f[E   ])[k1e   ];
			real mfcbb = (D.f[W   ])[k1w   ];
			real mfbab = (D.f[N   ])[k1n   ];
			real mfbcb = (D.f[S   ])[k1s   ];
			real mfbba = (D.f[T   ])[k1t   ];
			real mfbbc = (D.f[B   ])[k1b   ];
			real mfaab = (D.f[NE  ])[k1ne  ];
			real mfccb = (D.f[SW  ])[k1sw  ];
			real mfacb = (D.f[SE  ])[k1se  ];
			real mfcab = (D.f[NW  ])[k1nw  ];
			real mfaba = (D.f[TE  ])[k1te  ];
			real mfcbc = (D.f[BW  ])[k1bw  ];
			real mfabc = (D.f[BE  ])[k1be  ];
			real mfcba = (D.f[TW  ])[k1tw  ];
			real mfbaa = (D.f[TN  ])[k1tn  ];
			real mfbcc = (D.f[BS  ])[k1bs  ];
			real mfbac = (D.f[BN  ])[k1bn  ];
			real mfbca = (D.f[TS  ])[k1ts  ];
			real mfbbb = (D.f[REST])[k1zero];
			real mfaaa = (D.f[TNE ])[k1tne ];
			real mfcca = (D.f[TSW ])[k1tsw ];
			real mfaca = (D.f[TSE ])[k1tse ];
			real mfcaa = (D.f[TNW ])[k1tnw ];
			real mfaac = (D.f[BNE ])[k1bne ];
			real mfccc = (D.f[BSW ])[k1bsw ];
			real mfacc = (D.f[BSE ])[k1bse ];
			real mfcac = (D.f[BNW ])[k1bnw ];

			//real mfcbb = (D.f[E   ])[ke   ];
			//real mfabb = (D.f[W   ])[kw   ];
			//real mfbcb = (D.f[N   ])[kn   ];
			//real mfbab = (D.f[S   ])[ks   ];
			//real mfbbc = (D.f[T   ])[kt   ];
			//real mfbba = (D.f[B   ])[kb   ];
			//real mfccb = (D.f[NE  ])[kne  ];
			//real mfaab = (D.f[SW  ])[ksw  ];
			//real mfcab = (D.f[SE  ])[kse  ];
			//real mfacb = (D.f[NW  ])[knw  ];
			//real mfcbc = (D.f[TE  ])[kte  ];
			//real mfaba = (D.f[BW  ])[kbw  ];
			//real mfcba = (D.f[BE  ])[kbe  ];
			//real mfabc = (D.f[TW  ])[ktw  ];
			//real mfbcc = (D.f[TN  ])[ktn  ];
			//real mfbaa = (D.f[BS  ])[kbs  ];
			//real mfbca = (D.f[BN  ])[kbn  ];
			//real mfbac = (D.f[TS  ])[kts  ];
			//real mfbbb = (D.f[REST])[kzero];
			//real mfccc = (D.f[TNE ])[ktne ];
			//real mfaac = (D.f[TSW ])[ktsw ];
			//real mfcac = (D.f[TSE ])[ktse ];
			//real mfacc = (D.f[TNW ])[ktnw ];
			//real mfcca = (D.f[BNE ])[kbne ];
			//real mfaaa = (D.f[BSW ])[kbsw ];
			//real mfcaa = (D.f[BSE ])[kbse ];
			//real mfaca = (D.f[BNW ])[kbnw ];
			////////////////////////////////////////////////////////////////////////////////////
			//real rho   = (((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//				(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//				((mfabb+mfcbb) + (mfbab+mfbcb)) + (mfbba+mfbbc)) + mfbbb) + one;//!!!!Achtung + one
			////////////////////////////////////////////////////////////////////////////////////
			real rho = rhoBC[k];
			////////////////////////////////////////////////////////////////////////////////////
			real OoRho = c1o1 / (rho * 1.5f);
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) * OoRho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) * OoRho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) * OoRho;
			/////////////////////////
			//Test Values
			//double vvx    = 0.016;
			//double vvy    = zero;
			//double vvz    = zero;
			////////////////////////////////////////////////////////////////////////////////////////
			////round off error test
			//if(vvx!=zero){
			//	(kDistTest.f[E   ])[k] = mfabb;
			//	(kDistTest.f[W   ])[k] = mfcbb;
			//	(kDistTest.f[N   ])[k] = mfbab;
			//	(kDistTest.f[S   ])[k] = mfbcb;
			//	(kDistTest.f[T   ])[k] = mfbba;
			//	(kDistTest.f[B   ])[k] = mfbbc;
			//	(kDistTest.f[NE  ])[k] = mfaab;
			//	(kDistTest.f[SW  ])[k] = mfccb;
			//	(kDistTest.f[SE  ])[k] = mfacb;
			//	(kDistTest.f[NW  ])[k] = mfcab;
			//	(kDistTest.f[TE  ])[k] = mfaba;
			//	(kDistTest.f[BW  ])[k] = mfcbc;
			//	(kDistTest.f[BE  ])[k] = mfabc;
			//	(kDistTest.f[TW  ])[k] = mfcba;
			//	(kDistTest.f[TN  ])[k] = mfbaa;
			//	(kDistTest.f[BS  ])[k] = mfbcc;
			//	(kDistTest.f[BN  ])[k] = mfbac;
			//	(kDistTest.f[TS  ])[k] = mfbca;
			//	(kDistTest.f[REST])[k] = KQK;
			//	(kDistTest.f[TNE ])[k] = mfaaa;
			//	(kDistTest.f[TSW ])[k] = mfcca;
			//	(kDistTest.f[TSE ])[k] = mfaca;
			//	(kDistTest.f[TNW ])[k] = mfcaa;
			//	(kDistTest.f[BNE ])[k] = mfaac;
			//	(kDistTest.f[BSW ])[k] = mfccc;
			//	(kDistTest.f[BSE ])[k] = mfacc;
			//	(kDistTest.f[BNW ])[k] = mfcac;
			//}else{
			//	(kDistTest.f[E   ])[k] = zero;
			//	(kDistTest.f[W   ])[k] = zero;
			//	(kDistTest.f[N   ])[k] = zero;
			//	(kDistTest.f[S   ])[k] = zero;
			//	(kDistTest.f[T   ])[k] = zero;
			//	(kDistTest.f[B   ])[k] = zero;
			//	(kDistTest.f[NE  ])[k] = zero;
			//	(kDistTest.f[SW  ])[k] = zero;
			//	(kDistTest.f[SE  ])[k] = zero;
			//	(kDistTest.f[NW  ])[k] = zero;
			//	(kDistTest.f[TE  ])[k] = zero;
			//	(kDistTest.f[BW  ])[k] = zero;
			//	(kDistTest.f[BE  ])[k] = zero;
			//	(kDistTest.f[TW  ])[k] = zero;
			//	(kDistTest.f[TN  ])[k] = zero;
			//	(kDistTest.f[BS  ])[k] = zero;
			//	(kDistTest.f[BN  ])[k] = zero;
			//	(kDistTest.f[TS  ])[k] = zero;
			//	(kDistTest.f[REST])[k] = zero;
			//	(kDistTest.f[TNE ])[k] = zero;
			//	(kDistTest.f[TSW ])[k] = zero;
			//	(kDistTest.f[TSE ])[k] = zero;
			//	(kDistTest.f[TNW ])[k] = zero;
			//	(kDistTest.f[BNE ])[k] = zero;
			//	(kDistTest.f[BSW ])[k] = zero;
			//	(kDistTest.f[BSE ])[k] = zero;
			//	(kDistTest.f[BNW ])[k] = zero;
			//}

			//////////////////////////////////////////////////////////////////////////////////////
			//// first bad fix for negative x velocity
			////if(vvx > zero) vvx = zero;
			//////////////////////////////////////////////////////////////////////////////////////
			////// second bad fix for negative x velocity
			////if(vvx > zero){
			////	vvx = -vvx;
			////	vvy = -vvy;
			////	vvz = -vvz;
			////}
			////////////////////////////////////////////////////////////////////////////////////
			double vx2    = vvx * vvx;
			double vy2    = vvy * vvy;
			double vz2    = vvz * vvz;
			//////////////////////////////////////////////////////////////////////////////////
			//original
            real XXb    = -c2o3 + vx2;
            real XXc    = -c1o2 * (XXb + c1o1 + vvx);
            real XXa    = XXc + vvx;
            real YYb    = -c2o3 + vy2;
            real YYc    = -c1o2 * (YYb + c1o1 + vvy);
            real YYa    = YYc + vvy;
            real ZZb    = -c2o3 + vz2;
            real ZZc    = -c1o2 * (ZZb + c1o1 + vvz);
            real ZZa    = ZZc + vvz;
			//////////////////////////////////////////////////////////////////////////////////
			//unkonditioniert
            mfcbb = -(rhoBC[k] + c1o1) * XXc * YYb * ZZb - c2o27; 
			mfabb = -(rhoBC[k] + c1o1) * XXa * YYb * ZZb - c2o27;
			mfbcb = -(rhoBC[k] + c1o1) * XXb * YYc * ZZb - c2o27;
			mfbab = -(rhoBC[k] + c1o1) * XXb * YYa * ZZb - c2o27;
			mfbbc = -(rhoBC[k] + c1o1) * XXb * YYb * ZZc - c2o27;
			mfbba = -(rhoBC[k] + c1o1) * XXb * YYb * ZZa - c2o27;
			mfccb = -(rhoBC[k] + c1o1) * XXc * YYc * ZZb - c1o54;
			mfaab = -(rhoBC[k] + c1o1) * XXa * YYa * ZZb - c1o54;
			mfcab = -(rhoBC[k] + c1o1) * XXc * YYa * ZZb - c1o54;
			mfacb = -(rhoBC[k] + c1o1) * XXa * YYc * ZZb - c1o54;
			mfcbc = -(rhoBC[k] + c1o1) * XXc * YYb * ZZc - c1o54;
			mfaba = -(rhoBC[k] + c1o1) * XXa * YYb * ZZa - c1o54;
			mfcba = -(rhoBC[k] + c1o1) * XXc * YYb * ZZa - c1o54;
			mfabc = -(rhoBC[k] + c1o1) * XXa * YYb * ZZc - c1o54;
			mfbcc = -(rhoBC[k] + c1o1) * XXb * YYc * ZZc - c1o54;
			mfbaa = -(rhoBC[k] + c1o1) * XXb * YYa * ZZa - c1o54;
			mfbca = -(rhoBC[k] + c1o1) * XXb * YYc * ZZa - c1o54;
			mfbac = -(rhoBC[k] + c1o1) * XXb * YYa * ZZc - c1o54;
			mfbbb = -(rhoBC[k] + c1o1) * XXb * YYb * ZZb - c8o27;
			mfccc = -(rhoBC[k] + c1o1) * XXc * YYc * ZZc - c1o216;
			mfaac = -(rhoBC[k] + c1o1) * XXa * YYa * ZZc - c1o216;
			mfcac = -(rhoBC[k] + c1o1) * XXc * YYa * ZZc - c1o216;
			mfacc = -(rhoBC[k] + c1o1) * XXa * YYc * ZZc - c1o216;
			mfcca = -(rhoBC[k] + c1o1) * XXc * YYc * ZZa - c1o216;
			mfaaa = -(rhoBC[k] + c1o1) * XXa * YYa * ZZa - c1o216;
			mfcaa = -(rhoBC[k] + c1o1) * XXc * YYa * ZZa - c1o216;
			mfaca = -(rhoBC[k] + c1o1) * XXa * YYc * ZZa - c1o216;
			//////////////////////////////////////////////////////////
			////konditioniert
			//double OneOver216RhoPlusOne = c1over216*(rhoBC[k]+one);
			//double OnoOver216Rho        = c1over216*rhoBC[k];
			//mfcbb = OnoOver216Rho*sixteen + OneOver216RhoPlusOne*twelve*(-(two*vy2) - two*vz2 + three*vy2*vz2 + vvx*(-two + three*vy2)*(-two + three*vz2) + vx2*(-two + three*vy2)*(-two + three*vz2));
			//mfabb = OnoOver216Rho*sixteen - OneOver216RhoPlusOne*twelve*(two*vy2 + two*vz2 - three*vy2*vz2 + vvx*(-two + three*vy2)*(-two + three*vz2) + vx2*(-four + six*vy2 + six*vz2 - nine*vy2*vz2));
			//mfbcb = four*(-(four*OneOver216RhoPlusOne) + four*OnoOver216Rho + OneOver216RhoPlusOne*(-two + three*vx2)*(one + three*vvy + three*vy2)*(-two + three*vz2));
			//mfbab = four*(four*OnoOver216Rho - OneOver216RhoPlusOne*three*(vvy*(-two + three*vx2)*(-two + three*vz2) - one*vx2*(one + three*vy2)*(-two + three*vz2) + two*(-(two*vy2) + vz2 + three*vy2*vz2)));
			//mfbbc = four*(-(four*OneOver216RhoPlusOne) + four*OnoOver216Rho + OneOver216RhoPlusOne*(-two + three*vx2)*(-two + three*vy2)*(one + three*vvz + three*vz2));
			//mfbba = four*(four*OnoOver216Rho - OneOver216RhoPlusOne*three*(vvz*(-two + three*vx2)*(-two + three*vy2) - one*vx2*(-two + three*vy2)*(one + three*vz2) + two*(vy2 - two*vz2 + three*vy2*vz2)));
			//mfccb = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) - two*vy2 - six*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(one + three*vvy + three*vy2)*(-two + three*vz2))));
			//mfaab = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) - two*vy2 - six*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-two + three*vz2))));
			//mfcab = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 + two*vy2 + six*vx2*vy2 - one*vz2 - three*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-two + three*vz2)));
			//mfacb = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 + two*vy2 + six*vx2*vy2 - one*vz2 - three*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-two + three*vz2) + vvx*(one + three*vvy + three*vy2)*(-two + three*vz2)));
			//mfcbc = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) + vy2 + three*vx2*vy2 + vvz*(one + three*vx2)*(-two + three*vy2) - two*vz2 - six*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(one + three*vvz + three*vz2))));
			//mfaba = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(-(two*vx2) + vy2 + three*vx2*vy2 - one*vvz*(one + three*vx2)*(-two + three*vy2) - two*vz2 - six*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(-one + three*vvz - three*vz2))));
			//mfcba = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 - one*vy2 - three*vx2*vy2 + vvz*(one + three*vx2)*(-two + three*vy2) + two*vz2 + six*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(-one + three*vvz - three*vz2)));
			//mfabc = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(two*vx2 - one*vy2 - three*vx2*vy2 - one*vvz*(one + three*vx2)*(-two + three*vy2) + two*vz2 + six*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 + vvx*(-two + three*vy2)*(one + three*vvz + three*vz2)));
			//mfbcc = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(vx2 - two*vy2 + three*vx2*vy2 + vvz*(-two + three*vx2)*(one + three*vy2) - two*vz2 + three*vx2*vz2 - six*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(one + three*vvz + three*vz2))));
			//mfbaa = -(two*(-(OnoOver216Rho*two) + OneOver216RhoPlusOne*three*(vx2 - two*vy2 + three*vx2*vy2 - one*vvz*(-two + three*vx2)*(one + three*vy2) - two*vz2 + three*vx2*vz2 - six*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(-one + three*vvz - three*vz2))));
			//mfbca = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(-(one*vx2) + two*vy2 - three*vx2*vy2 + vvz*(-two + three*vx2)*(one + three*vy2) + two*vz2 - three*vx2*vz2 + six*vy2*vz2 - nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(-one + three*vvz - three*vz2)));
			//mfbac = two*(OnoOver216Rho*two + OneOver216RhoPlusOne*three*(-(one*vx2) + two*vy2 - three*vx2*vy2 - one*vvz*(-two + three*vx2)*(one + three*vy2) + two*vz2 - three*vx2*vz2 + six*vy2*vz2 - nine*vx2*vy2*vz2 + vvy*(-two + three*vx2)*(one + three*vvz + three*vz2)));
			//mfbbb = eight*(eight*OnoOver216Rho + OneOver216RhoPlusOne*three*(four*vy2 + four*vz2 - six*vy2*vz2 + vx2*(-two + three*vy2)*(-two + three*vz2)));
			//mfccc = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(one + three*vvz + three*vz2) + vvx*(one + three*vvy + three*vy2)*(one + three*vvz + three*vz2));
			//mfaac = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(one + three*vvz + three*vz2) + vvx*(-one + three*vvy - three*vy2)*(one + three*vvz + three*vz2));
			//mfcac = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(one + three*vvz + three*vz2) - one*vvx*(-one + three*vvy - three*vy2)*(one + three*vvz + three*vz2));
			//mfacc = OnoOver216Rho + OneOver216RhoPlusOne*three*(vvz + vx2 + three*vvz*vx2 + vy2 + three*vvz*vy2 + three*vx2*vy2 + nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(one + three*vvz + three*vz2) - one*vvx*(one + three*vvy + three*vy2)*(one + three*vvz + three*vz2));
			//mfcca = OnoOver216Rho + OneOver216RhoPlusOne*three*(-(one*vvz) + vx2 - three*vvz*vx2 + vy2 - three*vvz*vy2 + three*vx2*vy2 - nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) - one*vvx*(one + three*vvy + three*vy2)*(-one + three*vvz - three*vz2));
			//mfaaa = OnoOver216Rho - OneOver216RhoPlusOne*three*(vvz - one*vx2 + three*vvz*vx2 - one*vy2 + three*vvz*vy2 - three*vx2*vy2 + nine*vvz*vx2*vy2 - one*vz2 - three*vx2*vz2 - three*vy2*vz2 - nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-one + three*vvz - three*vz2));
			//mfcaa = OnoOver216Rho + OneOver216RhoPlusOne*three*(-(one*vvz) + vx2 - three*vvz*vx2 + vy2 - three*vvz*vy2 + three*vx2*vy2 - nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 + vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) + vvx*(-one + three*vvy - three*vy2)*(-one + three*vvz - three*vz2));
			//mfaca = OnoOver216Rho + OneOver216RhoPlusOne*three*(-(one*vvz) + vx2 - three*vvz*vx2 + vy2 - three*vvz*vy2 + three*vx2*vy2 - nine*vvz*vx2*vy2 + vz2 + three*vx2*vz2 + three*vy2*vz2 + nine*vx2*vy2*vz2 - one*vvy*(one + three*vx2)*(-one + three*vvz - three*vz2) + vvx*(one + three*vvy + three*vy2)*(-one + three*vvz - three*vz2));

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //if (isEvenTimestep==true)
      //{
      //   D.f[E   ] = &DD[E   *size_Mat];
      //   D.f[W   ] = &DD[W   *size_Mat];
      //   D.f[N   ] = &DD[N   *size_Mat];
      //   D.f[S   ] = &DD[S   *size_Mat];
      //   D.f[T   ] = &DD[T   *size_Mat];
      //   D.f[B   ] = &DD[B   *size_Mat];
      //   D.f[NE  ] = &DD[NE  *size_Mat];
      //   D.f[SW  ] = &DD[SW  *size_Mat];
      //   D.f[SE  ] = &DD[SE  *size_Mat];
      //   D.f[NW  ] = &DD[NW  *size_Mat];
      //   D.f[TE  ] = &DD[TE  *size_Mat];
      //   D.f[BW  ] = &DD[BW  *size_Mat];
      //   D.f[BE  ] = &DD[BE  *size_Mat];
      //   D.f[TW  ] = &DD[TW  *size_Mat];
      //   D.f[TN  ] = &DD[TN  *size_Mat];
      //   D.f[BS  ] = &DD[BS  *size_Mat];
      //   D.f[BN  ] = &DD[BN  *size_Mat];
      //   D.f[TS  ] = &DD[TS  *size_Mat];
      //   D.f[REST] = &DD[REST*size_Mat];
      //   D.f[TNE ] = &DD[TNE *size_Mat];
      //   D.f[TSW ] = &DD[TSW *size_Mat];
      //   D.f[TSE ] = &DD[TSE *size_Mat];
      //   D.f[TNW ] = &DD[TNW *size_Mat];
      //   D.f[BNE ] = &DD[BNE *size_Mat];
      //   D.f[BSW ] = &DD[BSW *size_Mat];
      //   D.f[BSE ] = &DD[BSE *size_Mat];
      //   D.f[BNW ] = &DD[BNW *size_Mat];
      //} 
      //else
      //{
      //   D.f[W   ] = &DD[E   *size_Mat];
      //   D.f[E   ] = &DD[W   *size_Mat];
      //   D.f[S   ] = &DD[N   *size_Mat];
      //   D.f[N   ] = &DD[S   *size_Mat];
      //   D.f[B   ] = &DD[T   *size_Mat];
      //   D.f[T   ] = &DD[B   *size_Mat];
      //   D.f[SW  ] = &DD[NE  *size_Mat];
      //   D.f[NE  ] = &DD[SW  *size_Mat];
      //   D.f[NW  ] = &DD[SE  *size_Mat];
      //   D.f[SE  ] = &DD[NW  *size_Mat];
      //   D.f[BW  ] = &DD[TE  *size_Mat];
      //   D.f[TE  ] = &DD[BW  *size_Mat];
      //   D.f[TW  ] = &DD[BE  *size_Mat];
      //   D.f[BE  ] = &DD[TW  *size_Mat];
      //   D.f[BS  ] = &DD[TN  *size_Mat];
      //   D.f[TN  ] = &DD[BS  *size_Mat];
      //   D.f[TS  ] = &DD[BN  *size_Mat];
      //   D.f[BN  ] = &DD[TS  *size_Mat];
      //   D.f[REST] = &DD[REST*size_Mat];
      //   D.f[TNE ] = &DD[BSW *size_Mat];
      //   D.f[TSW ] = &DD[BNE *size_Mat];
      //   D.f[TSE ] = &DD[BNW *size_Mat];
      //   D.f[TNW ] = &DD[BSE *size_Mat];
      //   D.f[BNE ] = &DD[TSW *size_Mat];
      //   D.f[BSW ] = &DD[TNE *size_Mat];
      //   D.f[BSE ] = &DD[TNW *size_Mat];
      //   D.f[BNW ] = &DD[TSE *size_Mat];
      //}
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //__syncthreads();

			(D.f[E   ])[ke   ] = mfabb;//mfcbb;
			(D.f[W   ])[kw   ] = mfcbb;//mfabb;
			(D.f[N   ])[kn   ] = mfbab;//mfbcb;
			(D.f[S   ])[ks   ] = mfbcb;//mfbab;
			(D.f[T   ])[kt   ] = mfbba;//mfbbc;
			(D.f[B   ])[kb   ] = mfbbc;//mfbba;
			(D.f[NE  ])[kne  ] = mfaab;//mfccb;
			(D.f[SW  ])[ksw  ] = mfccb;//mfaab;
			(D.f[SE  ])[kse  ] = mfacb;//mfcab;
			(D.f[NW  ])[knw  ] = mfcab;//mfacb;
			(D.f[TE  ])[kte  ] = mfaba;//mfcbc;
			(D.f[BW  ])[kbw  ] = mfcbc;//mfaba;
			(D.f[BE  ])[kbe  ] = mfabc;//mfcba;
			(D.f[TW  ])[ktw  ] = mfcba;//mfabc;
			(D.f[TN  ])[ktn  ] = mfbaa;//mfbcc;
			(D.f[BS  ])[kbs  ] = mfbcc;//mfbaa;
			(D.f[BN  ])[kbn  ] = mfbac;//mfbca;
			(D.f[TS  ])[kts  ] = mfbca;//mfbac;
			(D.f[REST])[kzero] = mfbbb;//mfbbb;
			(D.f[TNE ])[ktne ] = mfaaa;//mfccc;
			(D.f[TSW ])[ktsw ] = mfcca;//mfaac;
			(D.f[TSE ])[ktse ] = mfaca;//mfcac;
			(D.f[TNW ])[ktnw ] = mfcaa;//mfacc;
			(D.f[BNE ])[kbne ] = mfaac;//mfcca;
			(D.f[BSW ])[kbsw ] = mfccc;//mfaaa;
			(D.f[BSE ])[kbse ] = mfacc;//mfcaa;
			(D.f[BNW ])[kbnw ] = mfcac;//mfaca;
			//(D.f[E   ])[ke   ] = mfcbb;
			//(D.f[W   ])[kw   ] = mfabb;
			//(D.f[N   ])[kn   ] = mfbcb;
			//(D.f[S   ])[ks   ] = mfbab;
			//(D.f[T   ])[kt   ] = mfbbc;
			//(D.f[B   ])[kb   ] = mfbba;
			//(D.f[NE  ])[kne  ] = mfccb;
			//(D.f[SW  ])[ksw  ] = mfaab;
			//(D.f[SE  ])[kse  ] = mfcab;
			//(D.f[NW  ])[knw  ] = mfacb;
			//(D.f[TE  ])[kte  ] = mfcbc;
			//(D.f[BW  ])[kbw  ] = mfaba;
			//(D.f[BE  ])[kbe  ] = mfcba;
			//(D.f[TW  ])[ktw  ] = mfabc;
			//(D.f[TN  ])[ktn  ] = mfbcc;
			//(D.f[BS  ])[kbs  ] = mfbaa;
			//(D.f[BN  ])[kbn  ] = mfbca;
			//(D.f[TS  ])[kts  ] = mfbac;
			//(D.f[REST])[kzero] = mfbbb;
			//(D.f[TNE ])[ktne ] = mfccc;
			//(D.f[TSW ])[ktsw ] = mfaac;
			//(D.f[TSE ])[ktse ] = mfcac;
			//(D.f[TNW ])[ktnw ] = mfacc;
			//(D.f[BNE ])[kbne ] = mfcca;
			//(D.f[BSW ])[kbsw ] = mfaaa;
			//(D.f[BSE ])[kbse ] = mfcaa;
			//(D.f[BNW ])[kbnw ] = mfaca;

      //(D.f[E   ])[ke   ] = fE ;  //f1_E ;   //fW;    //fE ;  
      //(D.f[W   ])[kw   ] = fW ;  //f1_W ;   //fE;    //fW ;  
      //(D.f[N   ])[kn   ] = fN ;  //f1_N ;   //fS;    //fN ;  
      //(D.f[S   ])[ks   ] = fS ;  //f1_S ;   //fN;    //fS ;  
      //(D.f[T   ])[kt   ] = fT ;  //f1_T ;   //fB;    //fT ;  
      //(D.f[B   ])[kb   ] = fB ;  //f1_B ;   //fT;    //fB ;  
      //(D.f[NE  ])[kne  ] = fNE;  //f1_NE;   //fSW;   //fNE;  
      //(D.f[SW  ])[ksw  ] = fSW;  //f1_SW;   //fNE;   //fSW;  
      //(D.f[SE  ])[kse  ] = fSE;  //f1_SE;   //fNW;   //fSE;  
      //(D.f[NW  ])[knw  ] = fNW;  //f1_NW;   //fSE;   //fNW;  
      //(D.f[TE  ])[kte  ] = fTE;  //f1_TE;   //fBW;   //fTE;  
      //(D.f[BW  ])[kbw  ] = fBW;  //f1_BW;   //fTE;   //fBW;  
      //(D.f[BE  ])[kbe  ] = fBE;  //f1_BE;   //fTW;   //fBE;  
      //(D.f[TW  ])[ktw  ] = fTW;  //f1_TW;   //fBE;   //fTW;  
      //(D.f[TN  ])[ktn  ] = fTN;  //f1_TN;   //fBS;   //fTN;  
      //(D.f[BS  ])[kbs  ] = fBS;  //f1_BS;   //fTN;   //fBS;  
      //(D.f[BN  ])[kbn  ] = fBN;  //f1_BN;   //fTS;   //fBN;  
      //(D.f[TS  ])[kts  ] = fTS;  //f1_TS;   //fBN;   //fTS;  
      //(D.f[REST])[kzero] = fZERO;//f1_ZERO; //fZERO; //fZERO;
      //(D.f[TNE ])[ktne ] = fTNE; //f1_TNE;  //fBSW;  //fTNE; 
      //(D.f[BSW ])[kbsw ] = fBSW; //f1_BSW;  //fTNE;  //fBSW; 
      //(D.f[BNE ])[kbne ] = fBNE; //f1_BNE;  //fTSW;  //fBNE; 
      //(D.f[TSW ])[ktsw ] = fTSW; //f1_TSW;  //fBNE;  //fTSW; 
      //(D.f[TSE ])[ktse ] = fTSE; //f1_TSE;  //fBNW;  //fTSE; 
      //(D.f[BNW ])[kbnw ] = fBNW; //f1_BNW;  //fTSE;  //fBNW; 
      //(D.f[BSE ])[kbse ] = fBSE; //f1_BSE;  //fTNW;  //fBSE; 
      //(D.f[TNW ])[ktnw ] = fTNW; //f1_TNW;  //fBSE;  //fTNW; 
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceZero27(	 real* DD, 
												 int* k_Q, 
												 int numberOfBCnodes, 
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
												 bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //__syncthreads();
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      (D.f[E   ])[ke   ] =c0o1;
      (D.f[W   ])[kw   ] =c0o1;
      (D.f[N   ])[kn   ] =c0o1;
      (D.f[S   ])[ks   ] =c0o1;
      (D.f[T   ])[kt   ] =c0o1;
      (D.f[B   ])[kb   ] =c0o1;
      (D.f[NE  ])[kne  ] =c0o1;
      (D.f[SW  ])[ksw  ] =c0o1;
      (D.f[SE  ])[kse  ] =c0o1;
      (D.f[NW  ])[knw  ] =c0o1;
      (D.f[TE  ])[kte  ] =c0o1;
      (D.f[BW  ])[kbw  ] =c0o1;
      (D.f[BE  ])[kbe  ] =c0o1;
      (D.f[TW  ])[ktw  ] =c0o1;
      (D.f[TN  ])[ktn  ] =c0o1;
      (D.f[BS  ])[kbs  ] =c0o1;
      (D.f[BN  ])[kbn  ] =c0o1;
      (D.f[TS  ])[kts  ] =c0o1;
      (D.f[REST])[kzero] =c0o1;
      (D.f[TNE ])[ktne ] =c0o1;
      (D.f[TSW ])[ktsw ] =c0o1;
      (D.f[TSE ])[ktse ] =c0o1;
      (D.f[TNW ])[ktnw ] =c0o1;
      (D.f[BNE ])[kbne ] =c0o1;
      (D.f[BSW ])[kbsw ] =c0o1;
      (D.f[BSE ])[kbse ] =c0o1;
      (D.f[BNW ])[kbnw ] =c0o1;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceFake27(	 real* rhoBC,
												 real* DD, 
												 int* k_Q, 
												 int* k_N, 
												 int numberOfBCnodes, 
												 real om1, 
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
												 bool isEvenTimestep)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<numberOfBCnodes)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //index
      unsigned int KQK  = k_Q[k];
      unsigned int kzero= KQK;
      unsigned int ke   = KQK;
      unsigned int kw   = neighborX[KQK];
      unsigned int kn   = KQK;
      unsigned int ks   = neighborY[KQK];
      unsigned int kt   = KQK;
      unsigned int kb   = neighborZ[KQK];
      unsigned int ksw  = neighborY[kw];
      unsigned int kne  = KQK;
      unsigned int kse  = ks;
      unsigned int knw  = kw;
      unsigned int kbw  = neighborZ[kw];
      unsigned int kte  = KQK;
      unsigned int kbe  = kb;
      unsigned int ktw  = kw;
      unsigned int kbs  = neighborZ[ks];
      unsigned int ktn  = KQK;
      unsigned int kbn  = kb;
      unsigned int kts  = ks;
      unsigned int ktse = ks;
      unsigned int kbnw = kbw;
      unsigned int ktnw = kw;
      unsigned int kbse = kbs;
      unsigned int ktsw = ksw;
      unsigned int kbne = kb;
      unsigned int ktne = KQK;
      unsigned int kbsw = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
      //index1
      unsigned int K1QK  = k_N[k];
      unsigned int k1zero= K1QK;
      unsigned int k1e   = K1QK;
      unsigned int k1w   = neighborX[K1QK];
      unsigned int k1n   = K1QK;
      unsigned int k1s   = neighborY[K1QK];
      unsigned int k1t   = K1QK;
      unsigned int k1b   = neighborZ[K1QK];
      unsigned int k1sw  = neighborY[k1w];
      unsigned int k1ne  = K1QK;
      unsigned int k1se  = k1s;
      unsigned int k1nw  = k1w;
      unsigned int k1bw  = neighborZ[k1w];
      unsigned int k1te  = K1QK;
      unsigned int k1be  = k1b;
      unsigned int k1tw  = k1w;
      unsigned int k1bs  = neighborZ[k1s];
      unsigned int k1tn  = K1QK;
      unsigned int k1bn  = k1b;
      unsigned int k1ts  = k1s;
      unsigned int k1tse = k1s;
      unsigned int k1bnw = k1bw;
      unsigned int k1tnw = k1w;
      unsigned int k1bse = k1bs;
      unsigned int k1tsw = k1sw;
      unsigned int k1bne = k1b;
      unsigned int k1tne = K1QK;
      unsigned int k1bsw = neighborZ[k1sw];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (isEvenTimestep==false)
      {
         D.f[E   ] = &DD[E   *size_Mat];
         D.f[W   ] = &DD[W   *size_Mat];
         D.f[N   ] = &DD[N   *size_Mat];
         D.f[S   ] = &DD[S   *size_Mat];
         D.f[T   ] = &DD[T   *size_Mat];
         D.f[B   ] = &DD[B   *size_Mat];
         D.f[NE  ] = &DD[NE  *size_Mat];
         D.f[SW  ] = &DD[SW  *size_Mat];
         D.f[SE  ] = &DD[SE  *size_Mat];
         D.f[NW  ] = &DD[NW  *size_Mat];
         D.f[TE  ] = &DD[TE  *size_Mat];
         D.f[BW  ] = &DD[BW  *size_Mat];
         D.f[BE  ] = &DD[BE  *size_Mat];
         D.f[TW  ] = &DD[TW  *size_Mat];
         D.f[TN  ] = &DD[TN  *size_Mat];
         D.f[BS  ] = &DD[BS  *size_Mat];
         D.f[BN  ] = &DD[BN  *size_Mat];
         D.f[TS  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[TNE *size_Mat];
         D.f[TSW ] = &DD[TSW *size_Mat];
         D.f[TSE ] = &DD[TSE *size_Mat];
         D.f[TNW ] = &DD[TNW *size_Mat];
         D.f[BNE ] = &DD[BNE *size_Mat];
         D.f[BSW ] = &DD[BSW *size_Mat];
         D.f[BSE ] = &DD[BSE *size_Mat];
         D.f[BNW ] = &DD[BNW *size_Mat];
      } 
      else
      {
         D.f[W   ] = &DD[E   *size_Mat];
         D.f[E   ] = &DD[W   *size_Mat];
         D.f[S   ] = &DD[N   *size_Mat];
         D.f[N   ] = &DD[S   *size_Mat];
         D.f[B   ] = &DD[T   *size_Mat];
         D.f[T   ] = &DD[B   *size_Mat];
         D.f[SW  ] = &DD[NE  *size_Mat];
         D.f[NE  ] = &DD[SW  *size_Mat];
         D.f[NW  ] = &DD[SE  *size_Mat];
         D.f[SE  ] = &DD[NW  *size_Mat];
         D.f[BW  ] = &DD[TE  *size_Mat];
         D.f[TE  ] = &DD[BW  *size_Mat];
         D.f[TW  ] = &DD[BE  *size_Mat];
         D.f[BE  ] = &DD[TW  *size_Mat];
         D.f[BS  ] = &DD[TN  *size_Mat];
         D.f[TN  ] = &DD[BS  *size_Mat];
         D.f[TS  ] = &DD[BN  *size_Mat];
         D.f[BN  ] = &DD[TS  *size_Mat];
         D.f[REST] = &DD[REST*size_Mat];
         D.f[TNE ] = &DD[BSW *size_Mat];
         D.f[TSW ] = &DD[BNE *size_Mat];
         D.f[TSE ] = &DD[BNW *size_Mat];
         D.f[TNW ] = &DD[BSE *size_Mat];
         D.f[BNE ] = &DD[TSW *size_Mat];
         D.f[BSW ] = &DD[TNE *size_Mat];
         D.f[BSE ] = &DD[TNW *size_Mat];
         D.f[BNW ] = &DD[TSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      real        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
         f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[E   ])[k1e   ];
      f1_E    = (D.f[W   ])[k1w   ];
      f1_S    = (D.f[N   ])[k1n   ];
      f1_N    = (D.f[S   ])[k1s   ];
      f1_B    = (D.f[T   ])[k1t   ];
      f1_T    = (D.f[B   ])[k1b   ];
      f1_SW   = (D.f[NE  ])[k1ne  ];
      f1_NE   = (D.f[SW  ])[k1sw  ];
      f1_NW   = (D.f[SE  ])[k1se  ];
      f1_SE   = (D.f[NW  ])[k1nw  ];
      f1_BW   = (D.f[TE  ])[k1te  ];
      f1_TE   = (D.f[BW  ])[k1bw  ];
      f1_TW   = (D.f[BE  ])[k1be  ];
      f1_BE   = (D.f[TW  ])[k1tw  ];
      f1_BS   = (D.f[TN  ])[k1tn  ];
      f1_TN   = (D.f[BS  ])[k1bs  ];
      f1_TS   = (D.f[BN  ])[k1bn  ];
      f1_BN   = (D.f[TS  ])[k1ts  ];
      f1_ZERO = (D.f[REST])[k1zero];
      f1_BSW  = (D.f[TNE ])[k1tne ];
      f1_BNE  = (D.f[TSW ])[k1tsw ];
      f1_BNW  = (D.f[TSE ])[k1tse ];
      f1_BSE  = (D.f[TNW ])[k1tnw ];
      f1_TSW  = (D.f[BNE ])[k1bne ];
      f1_TNE  = (D.f[BSW ])[k1bsw ];
      f1_TNW  = (D.f[BSE ])[k1bse ];
      f1_TSE  = (D.f[BNW ])[k1bnw ];

      ////////////////////////////////////////////////////////////////////////////////
      real vx1, vx2, vx3;
      vx1    =  ((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                  ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
                  (f1_E - f1_W); 


      vx2    =   (-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                  ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
                  (f1_N - f1_S); 

      vx3    =   ((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
                  (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
                  (f1_T - f1_B); 

      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      real drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
         f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

	  //drho1 = (drho1 + rhoBC[k])/2.f;
	  drho1 = drho1 - rhoBC[k];

      __syncthreads();

      (D.f[E   ])[ke   ] = c2o27* (rhoBC[k]+c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cu_sq);
      (D.f[W   ])[kw   ] = c2o27* (rhoBC[k]+c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cu_sq);
      (D.f[N   ])[kn   ] = c2o27* (rhoBC[k]+c3o1*(    -vx2    )+c9o2*(    -vx2    )*(    -vx2    )-cu_sq);
      (D.f[S   ])[ks   ] = c2o27* (rhoBC[k]+c3o1*(     vx2    )+c9o2*(     vx2    )*(     vx2    )-cu_sq);
      (D.f[T   ])[kt   ] = c2o27* (rhoBC[k]+c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cu_sq);
      (D.f[B   ])[kb   ] = c2o27* (rhoBC[k]+c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cu_sq);
      (D.f[NE  ])[kne  ] = f1_SW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[SW  ])[ksw  ] = f1_NE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[SE  ])[kse  ] = f1_NW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[NW  ])[knw  ] = f1_SE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TE  ])[kte  ] = f1_BW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BW  ])[kbw  ] = f1_TE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BE  ])[kbe  ] = f1_TW  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TW  ])[ktw  ] = f1_BE  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TN  ])[ktn  ] = f1_BS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BS  ])[kbs  ] = f1_TN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[BN  ])[kbn  ] = f1_TS  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[TS  ])[kts  ] = f1_BN  -c1o54*drho1;	//  c1o100;  // zero;  //
      (D.f[REST])[kzero] = f1_ZERO-c8o27*drho1;	//  c1o100;  // zero;  //
      (D.f[TNE ])[ktne ] = f1_BSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[TSW ])[ktsw ] = f1_BNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[TSE ])[ktse ] = f1_BNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[TNW ])[ktnw ] = f1_BSE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BNE ])[kbne ] = f1_TSW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BSW ])[kbsw ] = f1_TNE -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BSE ])[kbse ] = f1_TNW -c1o216*drho1;	//  c1o100;  // zero;  //
      (D.f[BNW ])[kbnw ] = f1_TSE -c1o216*drho1;  //  c1o100;  // zero;  //      
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDevice27_IntBB(real* rho,
												real* DD, 
												int* k_Q, 
												real* QQ,
												unsigned int numberOfBCnodes, 
												real om1, 
												unsigned int* neighborX,
												unsigned int* neighborY,
												unsigned int* neighborZ,
												unsigned int size_Mat, 
												bool isEvenTimestep)
{
	Distributions27 D;
	if (isEvenTimestep==true)
	{
		D.f[E   ] = &DD[E   *size_Mat];
		D.f[W   ] = &DD[W   *size_Mat];
		D.f[N   ] = &DD[N   *size_Mat];
		D.f[S   ] = &DD[S   *size_Mat];
		D.f[T   ] = &DD[T   *size_Mat];
		D.f[B   ] = &DD[B   *size_Mat];
		D.f[NE  ] = &DD[NE  *size_Mat];
		D.f[SW  ] = &DD[SW  *size_Mat];
		D.f[SE  ] = &DD[SE  *size_Mat];
		D.f[NW  ] = &DD[NW  *size_Mat];
		D.f[TE  ] = &DD[TE  *size_Mat];
		D.f[BW  ] = &DD[BW  *size_Mat];
		D.f[BE  ] = &DD[BE  *size_Mat];
		D.f[TW  ] = &DD[TW  *size_Mat];
		D.f[TN  ] = &DD[TN  *size_Mat];
		D.f[BS  ] = &DD[BS  *size_Mat];
		D.f[BN  ] = &DD[BN  *size_Mat];
		D.f[TS  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[TNE *size_Mat];
		D.f[TSW ] = &DD[TSW *size_Mat];
		D.f[TSE ] = &DD[TSE *size_Mat];
		D.f[TNW ] = &DD[TNW *size_Mat];
		D.f[BNE ] = &DD[BNE *size_Mat];
		D.f[BSW ] = &DD[BSW *size_Mat];
		D.f[BSE ] = &DD[BSE *size_Mat];
		D.f[BNW ] = &DD[BNW *size_Mat];
	} 
	else
	{
		D.f[W   ] = &DD[E   *size_Mat];
		D.f[E   ] = &DD[W   *size_Mat];
		D.f[S   ] = &DD[N   *size_Mat];
		D.f[N   ] = &DD[S   *size_Mat];
		D.f[B   ] = &DD[T   *size_Mat];
		D.f[T   ] = &DD[B   *size_Mat];
		D.f[SW  ] = &DD[NE  *size_Mat];
		D.f[NE  ] = &DD[SW  *size_Mat];
		D.f[NW  ] = &DD[SE  *size_Mat];
		D.f[SE  ] = &DD[NW  *size_Mat];
		D.f[BW  ] = &DD[TE  *size_Mat];
		D.f[TE  ] = &DD[BW  *size_Mat];
		D.f[TW  ] = &DD[BE  *size_Mat];
		D.f[BE  ] = &DD[TW  *size_Mat];
		D.f[BS  ] = &DD[TN  *size_Mat];
		D.f[TN  ] = &DD[BS  *size_Mat];
		D.f[TS  ] = &DD[BN  *size_Mat];
		D.f[BN  ] = &DD[TS  *size_Mat];
		D.f[REST] = &DD[REST*size_Mat];
		D.f[TNE ] = &DD[BSW *size_Mat];
		D.f[TSW ] = &DD[BNE *size_Mat];
		D.f[TSE ] = &DD[BNW *size_Mat];
		D.f[TNW ] = &DD[BSE *size_Mat];
		D.f[BNE ] = &DD[TSW *size_Mat];
		D.f[BSW ] = &DD[TNE *size_Mat];
		D.f[BSE ] = &DD[TNW *size_Mat];
		D.f[BNW ] = &DD[TSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if(k < numberOfBCnodes)
	{
		////////////////////////////////////////////////////////////////////////////////
		//real VeloX = vx[k];
		//real VeloY = vy[k];
		//real VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////
		real *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[E   * numberOfBCnodes];
		q_dirW   = &QQ[W   * numberOfBCnodes];
		q_dirN   = &QQ[N   * numberOfBCnodes];
		q_dirS   = &QQ[S   * numberOfBCnodes];
		q_dirT   = &QQ[T   * numberOfBCnodes];
		q_dirB   = &QQ[B   * numberOfBCnodes];
		q_dirNE  = &QQ[NE  * numberOfBCnodes];
		q_dirSW  = &QQ[SW  * numberOfBCnodes];
		q_dirSE  = &QQ[SE  * numberOfBCnodes];
		q_dirNW  = &QQ[NW  * numberOfBCnodes];
		q_dirTE  = &QQ[TE  * numberOfBCnodes];
		q_dirBW  = &QQ[BW  * numberOfBCnodes];
		q_dirBE  = &QQ[BE  * numberOfBCnodes];
		q_dirTW  = &QQ[TW  * numberOfBCnodes];
		q_dirTN  = &QQ[TN  * numberOfBCnodes];
		q_dirBS  = &QQ[BS  * numberOfBCnodes];
		q_dirBN  = &QQ[BN  * numberOfBCnodes];
		q_dirTS  = &QQ[TS  * numberOfBCnodes];
		q_dirTNE = &QQ[TNE * numberOfBCnodes];
		q_dirTSW = &QQ[TSW * numberOfBCnodes];
		q_dirTSE = &QQ[TSE * numberOfBCnodes];
		q_dirTNW = &QQ[TNW * numberOfBCnodes];
		q_dirBNE = &QQ[BNE * numberOfBCnodes];
		q_dirBSW = &QQ[BSW * numberOfBCnodes];
		q_dirBSE = &QQ[BSE * numberOfBCnodes];
		q_dirBNW = &QQ[BNW * numberOfBCnodes];
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK  = k_Q[k];
		unsigned int kzero= KQK;
		unsigned int ke   = KQK;
		unsigned int kw   = neighborX[KQK];
		unsigned int kn   = KQK;
		unsigned int ks   = neighborY[KQK];
		unsigned int kt   = KQK;
		unsigned int kb   = neighborZ[KQK];
		unsigned int ksw  = neighborY[kw];
		unsigned int kne  = KQK;
		unsigned int kse  = ks;
		unsigned int knw  = kw;
		unsigned int kbw  = neighborZ[kw];
		unsigned int kte  = KQK;
		unsigned int kbe  = kb;
		unsigned int ktw  = kw;
		unsigned int kbs  = neighborZ[ks];
		unsigned int ktn  = KQK;
		unsigned int kbn  = kb;
		unsigned int kts  = ks;
		unsigned int ktse = ks;
		unsigned int kbnw = kbw;
		unsigned int ktnw = kw;
		unsigned int kbse = kbs;
		unsigned int ktsw = ksw;
		unsigned int kbne = kb;
		unsigned int ktne = KQK;
		unsigned int kbsw = neighborZ[ksw];
		////////////////////////////////////////////////////////////////////////////////
		real f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
			f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		f_W    = (D.f[E   ])[ke   ];
		f_E    = (D.f[W   ])[kw   ];
		f_S    = (D.f[N   ])[kn   ];
		f_N    = (D.f[S   ])[ks   ];
		f_B    = (D.f[T   ])[kt   ];
		f_T    = (D.f[B   ])[kb   ];
		f_SW   = (D.f[NE  ])[kne  ];
		f_NE   = (D.f[SW  ])[ksw  ];
		f_NW   = (D.f[SE  ])[kse  ];
		f_SE   = (D.f[NW  ])[knw  ];
		f_BW   = (D.f[TE  ])[kte  ];
		f_TE   = (D.f[BW  ])[kbw  ];
		f_TW   = (D.f[BE  ])[kbe  ];
		f_BE   = (D.f[TW  ])[ktw  ];
		f_BS   = (D.f[TN  ])[ktn  ];
		f_TN   = (D.f[BS  ])[kbs  ];
		f_TS   = (D.f[BN  ])[kbn  ];
		f_BN   = (D.f[TS  ])[kts  ];
		f_BSW  = (D.f[TNE ])[ktne ];
		f_BNE  = (D.f[TSW ])[ktsw ];
		f_BNW  = (D.f[TSE ])[ktse ];
		f_BSE  = (D.f[TNW ])[ktnw ];
		f_TSW  = (D.f[BNE ])[kbne ];
		f_TNE  = (D.f[BSW ])[kbsw ];
		f_TNW  = (D.f[BSE ])[kbse ];
		f_TSE  = (D.f[BNW ])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////
		real vx1, vx2, vx3, drho, feq, q;
		drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
			f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
			f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[REST])[kzero]); 

		vx1    = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
			((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
			(f_E - f_W))/(c1o1+drho); 


		vx2    =  ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
			((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
			(f_N - f_S))/(c1o1+drho); 

		vx3    =  (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
			(-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
			(f_T - f_B))/(c1o1+drho); 

		real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

		//////////////////////////////////////////////////////////////////////////
		if (isEvenTimestep==false)
		{
			D.f[E   ] = &DD[E   *size_Mat];
			D.f[W   ] = &DD[W   *size_Mat];
			D.f[N   ] = &DD[N   *size_Mat];
			D.f[S   ] = &DD[S   *size_Mat];
			D.f[T   ] = &DD[T   *size_Mat];
			D.f[B   ] = &DD[B   *size_Mat];
			D.f[NE  ] = &DD[NE  *size_Mat];
			D.f[SW  ] = &DD[SW  *size_Mat];
			D.f[SE  ] = &DD[SE  *size_Mat];
			D.f[NW  ] = &DD[NW  *size_Mat];
			D.f[TE  ] = &DD[TE  *size_Mat];
			D.f[BW  ] = &DD[BW  *size_Mat];
			D.f[BE  ] = &DD[BE  *size_Mat];
			D.f[TW  ] = &DD[TW  *size_Mat];
			D.f[TN  ] = &DD[TN  *size_Mat];
			D.f[BS  ] = &DD[BS  *size_Mat];
			D.f[BN  ] = &DD[BN  *size_Mat];
			D.f[TS  ] = &DD[TS  *size_Mat];
			D.f[REST] = &DD[REST*size_Mat];
			D.f[TNE ] = &DD[TNE *size_Mat];
			D.f[TSW ] = &DD[TSW *size_Mat];
			D.f[TSE ] = &DD[TSE *size_Mat];
			D.f[TNW ] = &DD[TNW *size_Mat];
			D.f[BNE ] = &DD[BNE *size_Mat];
			D.f[BSW ] = &DD[BSW *size_Mat];
			D.f[BSE ] = &DD[BSE *size_Mat];
			D.f[BNW ] = &DD[BNW *size_Mat];
		} 
		else
		{
			D.f[W   ] = &DD[E   *size_Mat];
			D.f[E   ] = &DD[W   *size_Mat];
			D.f[S   ] = &DD[N   *size_Mat];
			D.f[N   ] = &DD[S   *size_Mat];
			D.f[B   ] = &DD[T   *size_Mat];
			D.f[T   ] = &DD[B   *size_Mat];
			D.f[SW  ] = &DD[NE  *size_Mat];
			D.f[NE  ] = &DD[SW  *size_Mat];
			D.f[NW  ] = &DD[SE  *size_Mat];
			D.f[SE  ] = &DD[NW  *size_Mat];
			D.f[BW  ] = &DD[TE  *size_Mat];
			D.f[TE  ] = &DD[BW  *size_Mat];
			D.f[TW  ] = &DD[BE  *size_Mat];
			D.f[BE  ] = &DD[TW  *size_Mat];
			D.f[BS  ] = &DD[TN  *size_Mat];
			D.f[TN  ] = &DD[BS  *size_Mat];
			D.f[TS  ] = &DD[BN  *size_Mat];
			D.f[BN  ] = &DD[TS  *size_Mat];
			D.f[REST] = &DD[REST*size_Mat];
			D.f[TNE ] = &DD[BSW *size_Mat];
			D.f[TSW ] = &DD[BNE *size_Mat];
			D.f[TSE ] = &DD[BNW *size_Mat];
			D.f[TNW ] = &DD[BSE *size_Mat];
			D.f[BNE ] = &DD[TSW *size_Mat];
			D.f[BSW ] = &DD[TNE *size_Mat];
			D.f[BSE ] = &DD[TNW *size_Mat];
			D.f[BNW ] = &DD[TSE *size_Mat];
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Test
		//(D.f[REST])[k]=c1o10;
		real rhoDiff = drho - rho[k];
		real VeloX = vx1;
		real VeloY = vx2;
		real VeloZ = vx3;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		q = q_dirE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c9o2*( vx1        )*( vx1        )-cu_sq); 
			(D.f[W])[kw]=(c1o1-q)/(c1o1+q)*(f_E-f_W+(f_E+f_W-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_E+f_W)-c2o27*(rhoDiff + c6o1*( VeloX     )))/(c1o1+q);
		}

		q = q_dirW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c9o2*(-vx1        )*(-vx1        )-cu_sq); 
			(D.f[E])[ke]=(c1o1-q)/(c1o1+q)*(f_W-f_E+(f_W+f_E-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_W+f_E)-c2o27*(rhoDiff + c6o1*(-VeloX     )))/(c1o1+q);
		}

		q = q_dirN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c9o2*(     vx2    )*(     vx2    )-cu_sq); 
			(D.f[S])[ks]=(c1o1-q)/(c1o1+q)*(f_N-f_S+(f_N+f_S-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_N+f_S)-c2o27*(rhoDiff + c6o1*( VeloY     )))/(c1o1+q);
		}

		q = q_dirS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c9o2*(    -vx2    )*(    -vx2    )-cu_sq); 
			(D.f[N])[kn]=(c1o1-q)/(c1o1+q)*(f_S-f_N+(f_S+f_N-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_S+f_N)-c2o27*(rhoDiff + c6o1*(-VeloY     )))/(c1o1+q);
		}

		q = q_dirT[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c9o2*(         vx3)*(         vx3)-cu_sq); 
			(D.f[B])[kb]=(c1o1-q)/(c1o1+q)*(f_T-f_B+(f_T+f_B-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_T+f_B)-c2o27*(rhoDiff + c6o1*( VeloZ     )))/(c1o1+q);
		}

		q = q_dirB[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c2o27* (drho+c9o2*(        -vx3)*(        -vx3)-cu_sq); 
			(D.f[T])[kt]=(c1o1-q)/(c1o1+q)*(f_B-f_T+(f_B+f_T-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_B+f_T)-c2o27*(rhoDiff + c6o1*(-VeloZ     )))/(c1o1+q);
		}

		q = q_dirNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
			(D.f[SW])[ksw]=(c1o1-q)/(c1o1+q)*(f_NE-f_SW+(f_NE+f_SW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NE+f_SW)-c1o54*(rhoDiff + c6o1*(VeloX+VeloY)))/(c1o1+q);
		}

		q = q_dirSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
			(D.f[NE])[kne]=(c1o1-q)/(c1o1+q)*(f_SW-f_NE+(f_SW+f_NE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SW+f_NE)-c1o54*(rhoDiff + c6o1*(-VeloX-VeloY)))/(c1o1+q);
		}

		q = q_dirSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
			(D.f[NW])[knw]=(c1o1-q)/(c1o1+q)*(f_SE-f_NW+(f_SE+f_NW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_SE+f_NW)-c1o54*(rhoDiff + c6o1*( VeloX-VeloY)))/(c1o1+q);
		}

		q = q_dirNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
			(D.f[SE])[kse]=(c1o1-q)/(c1o1+q)*(f_NW-f_SE+(f_NW+f_SE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_NW+f_SE)-c1o54*(rhoDiff + c6o1*(-VeloX+VeloY)))/(c1o1+q);
		}

		q = q_dirTE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
			(D.f[BW])[kbw]=(c1o1-q)/(c1o1+q)*(f_TE-f_BW+(f_TE+f_BW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TE+f_BW)-c1o54*(rhoDiff + c6o1*( VeloX+VeloZ)))/(c1o1+q);
		}

		q = q_dirBW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
			(D.f[TE])[kte]=(c1o1-q)/(c1o1+q)*(f_BW-f_TE+(f_BW+f_TE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BW+f_TE)-c1o54*(rhoDiff + c6o1*(-VeloX-VeloZ)))/(c1o1+q);
		}

		q = q_dirBE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
			(D.f[TW])[ktw]=(c1o1-q)/(c1o1+q)*(f_BE-f_TW+(f_BE+f_TW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BE+f_TW)-c1o54*(rhoDiff + c6o1*( VeloX-VeloZ)))/(c1o1+q);
		}

		q = q_dirTW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
			(D.f[BE])[kbe]=(c1o1-q)/(c1o1+q)*(f_TW-f_BE+(f_TW+f_BE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TW+f_BE)-c1o54*(rhoDiff + c6o1*(-VeloX+VeloZ)))/(c1o1+q);
		}

		q = q_dirTN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
			(D.f[BS])[kbs]=(c1o1-q)/(c1o1+q)*(f_TN-f_BS+(f_TN+f_BS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TN+f_BS)-c1o54*(rhoDiff + c6o1*( VeloY+VeloZ)))/(c1o1+q);
		}

		q = q_dirBS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
			(D.f[TN])[ktn]=(c1o1-q)/(c1o1+q)*(f_BS-f_TN+(f_BS+f_TN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BS+f_TN)-c1o54*(rhoDiff + c6o1*( -VeloY-VeloZ)))/(c1o1+q);
		}

		q = q_dirBN[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
			(D.f[TS])[kts]=(c1o1-q)/(c1o1+q)*(f_BN-f_TS+(f_BN+f_TS-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BN+f_TS)-c1o54*(rhoDiff + c6o1*( VeloY-VeloZ)))/(c1o1+q);
		}

		q = q_dirTS[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o54* (drho+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
			(D.f[BN])[kbn]=(c1o1-q)/(c1o1+q)*(f_TS-f_BN+(f_TS+f_BN-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TS+f_BN)-c1o54*(rhoDiff + c6o1*( -VeloY+VeloZ)))/(c1o1+q);
		}

		q = q_dirTNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
			(D.f[BSW])[kbsw]=(c1o1-q)/(c1o1+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNE+f_BSW)-c1o216*(rhoDiff + c6o1*( VeloX+VeloY+VeloZ)))/(c1o1+q);
		}

		q = q_dirBSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
			(D.f[TNE])[ktne]=(c1o1-q)/(c1o1+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSW+f_TNE)-c1o216*(rhoDiff + c6o1*(-VeloX-VeloY-VeloZ)))/(c1o1+q);
		}

		q = q_dirBNE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
			(D.f[TSW])[ktsw]=(c1o1-q)/(c1o1+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNE+f_TSW)-c1o216*(rhoDiff + c6o1*( VeloX+VeloY-VeloZ)))/(c1o1+q);
		}

		q = q_dirTSW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
			(D.f[BNE])[kbne]=(c1o1-q)/(c1o1+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSW+f_BNE)-c1o216*(rhoDiff + c6o1*(-VeloX-VeloY+VeloZ)))/(c1o1+q);
		}

		q = q_dirTSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
			(D.f[BNW])[kbnw]=(c1o1-q)/(c1o1+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TSE+f_BNW)-c1o216*(rhoDiff + c6o1*( VeloX-VeloY+VeloZ)))/(c1o1+q);
		}

		q = q_dirBNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
			(D.f[TSE])[ktse]=(c1o1-q)/(c1o1+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BNW+f_TSE)-c1o216*(rhoDiff + c6o1*(-VeloX+VeloY-VeloZ)))/(c1o1+q);
		}

		q = q_dirBSE[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
			(D.f[TNW])[ktnw]=(c1o1-q)/(c1o1+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_BSE+f_TNW)-c1o216*(rhoDiff + c6o1*( VeloX-VeloY-VeloZ)))/(c1o1+q);
		}

		q = q_dirTNW[k];
		if (q>=c0o1 && q<=c1o1)
		{
			feq=c1o216*(drho+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
			(D.f[BSE])[kbse]=(c1o1-q)/(c1o1+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-c2o1*feq*om1)/(c1o1-om1))*c1o2+(q*(f_TNW+f_BSE)-c1o216*(rhoDiff + c6o1*(-VeloX+VeloY+VeloZ)))/(c1o1+q);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


