/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QInflowScaleByPressDevice27(  doubflo* rhoBC,
														 doubflo* DD, 
														 int* k_Q, 
														 int* k_N, 
														 int kQ, 
														 doubflo om1, 
														 unsigned int* neighborX,
														 unsigned int* neighborY,
														 unsigned int* neighborZ,
														 unsigned int size_Mat, 
														 bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==true)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo f1_E    = (D.f[dirE   ])[k1e   ];
      doubflo f1_W    = (D.f[dirW   ])[k1w   ];
      doubflo f1_N    = (D.f[dirN   ])[k1n   ];
      doubflo f1_S    = (D.f[dirS   ])[k1s   ];
      doubflo f1_T    = (D.f[dirT   ])[k1t   ];
      doubflo f1_B    = (D.f[dirB   ])[k1b   ];
      doubflo f1_NE   = (D.f[dirNE  ])[k1ne  ];
      doubflo f1_SW   = (D.f[dirSW  ])[k1sw  ];
      doubflo f1_SE   = (D.f[dirSE  ])[k1se  ];
      doubflo f1_NW   = (D.f[dirNW  ])[k1nw  ];
      doubflo f1_TE   = (D.f[dirTE  ])[k1te  ];
      doubflo f1_BW   = (D.f[dirBW  ])[k1bw  ];
      doubflo f1_BE   = (D.f[dirBE  ])[k1be  ];
      doubflo f1_TW   = (D.f[dirTW  ])[k1tw  ];
      doubflo f1_TN   = (D.f[dirTN  ])[k1tn  ];
      doubflo f1_BS   = (D.f[dirBS  ])[k1bs  ];
      doubflo f1_BN   = (D.f[dirBN  ])[k1bn  ];
      doubflo f1_TS   = (D.f[dirTS  ])[k1ts  ];
      doubflo f1_ZERO = (D.f[dirZERO])[k1zero];
      doubflo f1_TNE  = (D.f[dirTNE ])[k1tne ];
      doubflo f1_TSW  = (D.f[dirTSW ])[k1tsw ];
      doubflo f1_TSE  = (D.f[dirTSE ])[k1tse ];
      doubflo f1_TNW  = (D.f[dirTNW ])[k1tnw ];
      doubflo f1_BNE  = (D.f[dirBNE ])[k1bne ];
      doubflo f1_BSW  = (D.f[dirBSW ])[k1bsw ];
      doubflo f1_BSE  = (D.f[dirBSE ])[k1bse ];
      doubflo f1_BNW  = (D.f[dirBNW ])[k1bnw ];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo f_E    = (D.f[dirE   ])[ke   ];
      doubflo f_W    = (D.f[dirW   ])[kw   ];
      doubflo f_N    = (D.f[dirN   ])[kn   ];
      doubflo f_S    = (D.f[dirS   ])[ks   ];
      doubflo f_T    = (D.f[dirT   ])[kt   ];
      doubflo f_B    = (D.f[dirB   ])[kb   ];
      doubflo f_NE   = (D.f[dirNE  ])[kne  ];
      doubflo f_SW   = (D.f[dirSW  ])[ksw  ];
      doubflo f_SE   = (D.f[dirSE  ])[kse  ];
      doubflo f_NW   = (D.f[dirNW  ])[knw  ];
      doubflo f_TE   = (D.f[dirTE  ])[kte  ];
      doubflo f_BW   = (D.f[dirBW  ])[kbw  ];
      doubflo f_BE   = (D.f[dirBE  ])[kbe  ];
      doubflo f_TW   = (D.f[dirTW  ])[ktw  ];
      doubflo f_TN   = (D.f[dirTN  ])[ktn  ];
      doubflo f_BS   = (D.f[dirBS  ])[kbs  ];
      doubflo f_BN   = (D.f[dirBN  ])[kbn  ];
      doubflo f_TS   = (D.f[dirTS  ])[kts  ];
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_TNE  = (D.f[dirTNE ])[ktne ];
      doubflo f_TSW  = (D.f[dirTSW ])[ktsw ];
      doubflo f_TSE  = (D.f[dirTSE ])[ktse ];
      doubflo f_TNW  = (D.f[dirTNW ])[ktnw ];
      doubflo f_BNE  = (D.f[dirBNE ])[kbne ];
      doubflo f_BSW  = (D.f[dirBSW ])[kbsw ];
      doubflo f_BSE  = (D.f[dirBSE ])[kbse ];
      doubflo f_BNW  = (D.f[dirBNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3, drho, drho1;
      //////////////////////////////////////////////////////////////////////////
	  //Dichte
      drho1  =  f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW + 
                f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((D.f[dirZERO])[k1zero]); 
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 
      //////////////////////////////////////////////////////////////////////////
	  //Schallgeschwindigkeit
	  doubflo cs = one / sqrtf(three);
      //////////////////////////////////////////////////////////////////////////
	  doubflo rhoInterpol = drho1 * cs + (one - cs) * drho; 
	  //doubflo diffRho = (rhoBC[k] + one) / (rhoInterpol + one);
	  doubflo diffRhoToAdd = rhoBC[k] - rhoInterpol;
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //no velocity
	  //////////////////////////////////////////
      f_E    = f1_E   * cs + (one - cs) * f_E   ;
      f_W    = f1_W   * cs + (one - cs) * f_W   ;
      f_N    = f1_N   * cs + (one - cs) * f_N   ;
      f_S    = f1_S   * cs + (one - cs) * f_S   ;
      f_T    = f1_T   * cs + (one - cs) * f_T   ;
      f_B    = f1_B   * cs + (one - cs) * f_B   ;
      f_NE   = f1_NE  * cs + (one - cs) * f_NE  ;
      f_SW   = f1_SW  * cs + (one - cs) * f_SW  ;
      f_SE   = f1_SE  * cs + (one - cs) * f_SE  ;
      f_NW   = f1_NW  * cs + (one - cs) * f_NW  ;
      f_TE   = f1_TE  * cs + (one - cs) * f_TE  ;
      f_BW   = f1_BW  * cs + (one - cs) * f_BW  ;
      f_BE   = f1_BE  * cs + (one - cs) * f_BE  ;
      f_TW   = f1_TW  * cs + (one - cs) * f_TW  ;
      f_TN   = f1_TN  * cs + (one - cs) * f_TN  ;
      f_BS   = f1_BS  * cs + (one - cs) * f_BS  ;
      f_BN   = f1_BN  * cs + (one - cs) * f_BN  ;
      f_TS   = f1_TS  * cs + (one - cs) * f_TS  ;
      f_TNE  = f1_TNE * cs + (one - cs) * f_TNE ;
      f_TSW  = f1_TSW * cs + (one - cs) * f_TSW ;
      f_TSE  = f1_TSE * cs + (one - cs) * f_TSE ;
      f_TNW  = f1_TNW * cs + (one - cs) * f_TNW ;
      f_BNE  = f1_BNE * cs + (one - cs) * f_BNE ;
      f_BSW  = f1_BSW * cs + (one - cs) * f_BSW ;
      f_BSE  = f1_BSE * cs + (one - cs) * f_BSE ;
      f_BNW  = f1_BNW * cs + (one - cs) * f_BNW ;
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
	  f_E    = (f_E   + c2over27  * diffRhoToAdd);
      f_W    = (f_W   + c2over27  * diffRhoToAdd);
      f_N    = (f_N   + c2over27  * diffRhoToAdd);
      f_S    = (f_S   + c2over27  * diffRhoToAdd);
      f_T    = (f_T   + c2over27  * diffRhoToAdd);
      f_B    = (f_B   + c2over27  * diffRhoToAdd);
	  f_NE   = (f_NE  + c1over54  * diffRhoToAdd);
      f_SW   = (f_SW  + c1over54  * diffRhoToAdd);
      f_SE   = (f_SE  + c1over54  * diffRhoToAdd);
      f_NW   = (f_NW  + c1over54  * diffRhoToAdd);
      f_TE   = (f_TE  + c1over54  * diffRhoToAdd);
      f_BW   = (f_BW  + c1over54  * diffRhoToAdd);
      f_BE   = (f_BE  + c1over54  * diffRhoToAdd);
      f_TW   = (f_TW  + c1over54  * diffRhoToAdd);
      f_TN   = (f_TN  + c1over54  * diffRhoToAdd);
      f_BS   = (f_BS  + c1over54  * diffRhoToAdd);
      f_BN   = (f_BN  + c1over54  * diffRhoToAdd);
      f_TS   = (f_TS  + c1over54  * diffRhoToAdd);
      f_TNE  = (f_TNE + c1over216 * diffRhoToAdd);
      f_TSW  = (f_TSW + c1over216 * diffRhoToAdd);
      f_TSE  = (f_TSE + c1over216 * diffRhoToAdd);
      f_TNW  = (f_TNW + c1over216 * diffRhoToAdd);
      f_BNE  = (f_BNE + c1over216 * diffRhoToAdd);
      f_BSW  = (f_BSW + c1over216 * diffRhoToAdd);
      f_BSE  = (f_BSE + c1over216 * diffRhoToAdd);
      f_BNW  = (f_BNW + c1over216 * diffRhoToAdd);
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////
      //__syncthreads();
	  // -X
	  //(D.f[dirE   ])[ke   ] = f_E   ;
	  //(D.f[dirSE  ])[kse  ] = f_SE  ;
	  //(D.f[dirNE  ])[kne  ] = f_NE  ;
	  //(D.f[dirBE  ])[kbe  ] = f_BE  ;
	  //(D.f[dirTE  ])[kte  ] = f_TE  ;
	  //(D.f[dirTSE ])[ktse ] = f_TSE ;
	  //(D.f[dirTNE ])[ktne ] = f_TNE ;
	  //(D.f[dirBSE ])[kbse ] = f_BSE ;
	  //(D.f[dirBNE ])[kbne ] = f_BNE ;     
	  // X
	  (D.f[dirW   ])[kw   ] = f_W   ;
	  (D.f[dirSW  ])[ksw  ] = f_SW  ;
	  (D.f[dirNW  ])[knw  ] = f_NW  ;
	  (D.f[dirBW  ])[kbw  ] = f_BW  ;
	  (D.f[dirTW  ])[ktw  ] = f_TW  ;
	  (D.f[dirTSW ])[ktsw ] = f_TSW ;
	  (D.f[dirTNW ])[ktnw ] = f_TNW ;
	  (D.f[dirBSW ])[kbsw ] = f_BSW ;
	  (D.f[dirBNW ])[kbnw ] = f_BNW ;     
	  // Y
	  //(D.f[dirS   ])[ks   ] = f_S   ;
	  //(D.f[dirSE  ])[kse  ] = f_SE  ;
	  //(D.f[dirSW  ])[ksw  ] = f_SW  ;
	  //(D.f[dirTS  ])[kts  ] = f_TS  ;
	  //(D.f[dirBS  ])[kbs  ] = f_BS  ;
	  //(D.f[dirTSE ])[ktse ] = f_TSE ;
	  //(D.f[dirTSW ])[ktsw ] = f_TSW ;
	  //(D.f[dirBSE ])[kbse ] = f_BSE ;
	  //(D.f[dirBSW ])[kbsw ] = f_BSW ;     
	  // Z
	  //(D.f[dirB   ])[kb   ] = f_B   ;
	  //(D.f[dirBE  ])[kbe  ] = f_BE  ;
	  //(D.f[dirBW  ])[kbw  ] = f_BW  ;
	  //(D.f[dirBN  ])[kbn  ] = f_BN  ;
	  //(D.f[dirBS  ])[kbs  ] = f_BS  ;
	  //(D.f[dirBNE ])[kbne ] = f_BNE ;
	  //(D.f[dirBNW ])[kbnw ] = f_BNW ;
	  //(D.f[dirBSE ])[kbse ] = f_BSE ;
	  //(D.f[dirBSW ])[kbsw ] = f_BSW ;     
      //////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceIncompNEQ27( doubflo* rhoBC,
													doubflo* DD, 
													int* k_Q, 
													int* k_N, 
													int kQ, 
													doubflo om1, 
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat, 
													bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==true) //// ACHTUNG PREColl !!!!!!!!!!!!!!
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dirE   ])[k1e   ];
      f1_E    = (D.f[dirW   ])[k1w   ];
      f1_S    = (D.f[dirN   ])[k1n   ];
      f1_N    = (D.f[dirS   ])[k1s   ];
      f1_B    = (D.f[dirT   ])[k1t   ];
      f1_T    = (D.f[dirB   ])[k1b   ];
      f1_SW   = (D.f[dirNE  ])[k1ne  ];
      f1_NE   = (D.f[dirSW  ])[k1sw  ];
      f1_NW   = (D.f[dirSE  ])[k1se  ];
      f1_SE   = (D.f[dirNW  ])[k1nw  ];
      f1_BW   = (D.f[dirTE  ])[k1te  ];
      f1_TE   = (D.f[dirBW  ])[k1bw  ];
      f1_TW   = (D.f[dirBE  ])[k1be  ];
      f1_BE   = (D.f[dirTW  ])[k1tw  ];
      f1_BS   = (D.f[dirTN  ])[k1tn  ];
      f1_TN   = (D.f[dirBS  ])[k1bs  ];
      f1_TS   = (D.f[dirBN  ])[k1bn  ];
      f1_BN   = (D.f[dirTS  ])[k1ts  ];
      f1_ZERO = (D.f[dirZERO])[k1zero];
      f1_BSW  = (D.f[dirTNE ])[k1tne ];
      f1_BNE  = (D.f[dirTSW ])[k1tsw ];
      f1_BNW  = (D.f[dirTSE ])[k1tse ];
      f1_BSE  = (D.f[dirTNW ])[k1tnw ];
      f1_TSW  = (D.f[dirBNE ])[k1bne ];
      f1_TNE  = (D.f[dirBSW ])[k1bsw ];
      f1_TNW  = (D.f[dirBSE ])[k1bse ];
      f1_TSE  = (D.f[dirBNW ])[k1bnw ];

      //////////////////////////////////////////////////////////////////////////
      doubflo drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                          f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      doubflo vx1      =  ((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
						  ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
						  (f1_E - f1_W); 


      doubflo vx2    =   (-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
						 ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
						 (f1_N - f1_S); 

      doubflo vx3    =   ((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
						 (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
						 (f1_T - f1_B); 

      doubflo cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      f1_ZERO  -= c8over27*  (drho1-(drho1+one)*cusq);
      f1_E     -= c2over27*  (drho1+(drho1+one)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
      f1_W     -= c2over27*  (drho1+(drho1+one)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
      f1_N     -= c2over27*  (drho1+(drho1+one)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
      f1_S     -= c2over27*  (drho1+(drho1+one)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     -= c2over27*  (drho1+(drho1+one)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
      f1_B     -= c2over27*  (drho1+(drho1+one)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    -= c1over54*  (drho1+(drho1+one)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    -= c1over54*  (drho1+(drho1+one)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    -=  c1over54* (drho1+(drho1+one)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    -=  c1over54* (drho1+(drho1+one)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    -=  c1over54* (drho1+(drho1+one)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    -=  c1over54* (drho1+(drho1+one)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    -=  c1over54* (drho1+(drho1+one)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    -=  c1over54* (drho1+(drho1+one)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    -=  c1over54* (drho1+(drho1+one)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    -=  c1over54* (drho1+(drho1+one)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    -=  c1over54* (drho1+(drho1+one)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    -=  c1over54* (drho1+(drho1+one)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
	   
	  drho1 = rhoBC[k];

	  //if(vx1 < zero){
		 // vx1 *= 0.9;
	  //}
	  //if(vx2 < zero){
		 // vx2 *= c1o10;//0.9;
	  //}

      f1_ZERO  += c8over27*  (drho1-(drho1+one)*cusq);
      f1_E     += c2over27*  (drho1+(drho1+one)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
      f1_W     += c2over27*  (drho1+(drho1+one)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
      f1_N     += c2over27*  (drho1+(drho1+one)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
      f1_S     += c2over27*  (drho1+(drho1+one)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     += c2over27*  (drho1+(drho1+one)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
      f1_B     += c2over27*  (drho1+(drho1+one)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    += c1over54*  (drho1+(drho1+one)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    += c1over54*  (drho1+(drho1+one)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    +=  c1over54* (drho1+(drho1+one)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    +=  c1over54* (drho1+(drho1+one)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    +=  c1over54* (drho1+(drho1+one)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    +=  c1over54* (drho1+(drho1+one)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    +=  c1over54* (drho1+(drho1+one)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    +=  c1over54* (drho1+(drho1+one)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    +=  c1over54* (drho1+(drho1+one)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    +=  c1over54* (drho1+(drho1+one)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    +=  c1over54* (drho1+(drho1+one)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    +=  c1over54* (drho1+(drho1+one)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

	  //drho1 = (drho1 + rhoBC[k])/2.f;
	  //drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[dirE   ])[ke   ] = f1_W   ;  
      (D.f[dirW   ])[kw   ] = f1_E   ;	
      (D.f[dirN   ])[kn   ] = f1_S   ;	
      (D.f[dirS   ])[ks   ] = f1_N   ;	
      (D.f[dirT   ])[kt   ] = f1_B   ;	
      (D.f[dirB   ])[kb   ] = f1_T   ;	
      (D.f[dirNE  ])[kne  ] = f1_SW  ;	
      (D.f[dirSW  ])[ksw  ] = f1_NE  ;	
      (D.f[dirSE  ])[kse  ] = f1_NW  ;	
      (D.f[dirNW  ])[knw  ] = f1_SE  ;	
      (D.f[dirTE  ])[kte  ] = f1_BW  ;	
      (D.f[dirBW  ])[kbw  ] = f1_TE  ;	
      (D.f[dirBE  ])[kbe  ] = f1_TW  ;	
      (D.f[dirTW  ])[ktw  ] = f1_BE  ;	
      (D.f[dirTN  ])[ktn  ] = f1_BS  ;	
      (D.f[dirBS  ])[kbs  ] = f1_TN  ;	
      (D.f[dirBN  ])[kbn  ] = f1_TS  ;	
      (D.f[dirTS  ])[kts  ] = f1_BN  ;	
      (D.f[dirZERO])[kzero] = f1_ZERO;	
      (D.f[dirTNE ])[ktne ] = f1_BSW ;	
      (D.f[dirTSW ])[ktsw ] = f1_BNE ;	
      (D.f[dirTSE ])[ktse ] = f1_BNW ;	
      (D.f[dirTNW ])[ktnw ] = f1_BSE ;	
      (D.f[dirBNE ])[kbne ] = f1_TSW ;	
      (D.f[dirBSW ])[kbsw ] = f1_TNE ;	
      (D.f[dirBSE ])[kbse ] = f1_TNW ;	
      (D.f[dirBNW ])[kbnw ] = f1_TSE ;       
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceNEQ27(doubflo* rhoBC,
                                             doubflo* DD, 
                                             int* k_Q, 
                                             int* k_N, 
                                             int kQ, 
                                             doubflo om1, 
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int size_Mat, 
                                             bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==true) //// ACHTUNG PREColl !!!!!!!!!!!!!!
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dirE   ])[k1e   ];
      f1_E    = (D.f[dirW   ])[k1w   ];
      f1_S    = (D.f[dirN   ])[k1n   ];
      f1_N    = (D.f[dirS   ])[k1s   ];
      f1_B    = (D.f[dirT   ])[k1t   ];
      f1_T    = (D.f[dirB   ])[k1b   ];
      f1_SW   = (D.f[dirNE  ])[k1ne  ];
      f1_NE   = (D.f[dirSW  ])[k1sw  ];
      f1_NW   = (D.f[dirSE  ])[k1se  ];
      f1_SE   = (D.f[dirNW  ])[k1nw  ];
      f1_BW   = (D.f[dirTE  ])[k1te  ];
      f1_TE   = (D.f[dirBW  ])[k1bw  ];
      f1_TW   = (D.f[dirBE  ])[k1be  ];
      f1_BE   = (D.f[dirTW  ])[k1tw  ];
      f1_BS   = (D.f[dirTN  ])[k1tn  ];
      f1_TN   = (D.f[dirBS  ])[k1bs  ];
      f1_TS   = (D.f[dirBN  ])[k1bn  ];
      f1_BN   = (D.f[dirTS  ])[k1ts  ];
      f1_ZERO = (D.f[dirZERO])[k1zero];
      f1_BSW  = (D.f[dirTNE ])[k1tne ];
      f1_BNE  = (D.f[dirTSW ])[k1tsw ];
      f1_BNW  = (D.f[dirTSE ])[k1tse ];
      f1_BSE  = (D.f[dirTNW ])[k1tnw ];
      f1_TSW  = (D.f[dirBNE ])[k1bne ];
      f1_TNE  = (D.f[dirBSW ])[k1bsw ];
      f1_TNW  = (D.f[dirBSE ])[k1bse ];
      f1_TSE  = (D.f[dirBNW ])[k1bnw ];

      //////////////////////////////////////////////////////////////////////////
      doubflo drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                          f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      doubflo vx1      =  ((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
						  ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
						  (f1_E - f1_W); 


      doubflo vx2    =   (-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
						 ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
						 (f1_N - f1_S); 

      doubflo vx3    =   ((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
						 (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
						 (f1_T - f1_B); 

	  vx1 /= (drho1+one);
	  vx2 /= (drho1+one);
	  vx3 /= (drho1+one);

      doubflo cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      f1_ZERO  -= c8over27*  (drho1-(drho1+one)*cusq);
      f1_E     -= c2over27*  (drho1+(drho1+one)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
      f1_W     -= c2over27*  (drho1+(drho1+one)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
      f1_N     -= c2over27*  (drho1+(drho1+one)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
      f1_S     -= c2over27*  (drho1+(drho1+one)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     -= c2over27*  (drho1+(drho1+one)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
      f1_B     -= c2over27*  (drho1+(drho1+one)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    -= c1over54*  (drho1+(drho1+one)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    -= c1over54*  (drho1+(drho1+one)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    -=  c1over54* (drho1+(drho1+one)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    -=  c1over54* (drho1+(drho1+one)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    -=  c1over54* (drho1+(drho1+one)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    -=  c1over54* (drho1+(drho1+one)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    -=  c1over54* (drho1+(drho1+one)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    -=  c1over54* (drho1+(drho1+one)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    -=  c1over54* (drho1+(drho1+one)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    -=  c1over54* (drho1+(drho1+one)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    -=  c1over54* (drho1+(drho1+one)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    -=  c1over54* (drho1+(drho1+one)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   -=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   -=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
	   
	  drho1 = rhoBC[k];

      f1_ZERO  += c8over27*  (drho1-(drho1+one)*cusq);
      f1_E     += c2over27*  (drho1+(drho1+one)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
      f1_W     += c2over27*  (drho1+(drho1+one)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
      f1_N     += c2over27*  (drho1+(drho1+one)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
      f1_S     += c2over27*  (drho1+(drho1+one)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
      f1_T     += c2over27*  (drho1+(drho1+one)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
      f1_B     += c2over27*  (drho1+(drho1+one)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
      f1_NE    += c1over54*  (drho1+(drho1+one)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f1_SW    += c1over54*  (drho1+(drho1+one)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f1_SE    +=  c1over54* (drho1+(drho1+one)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f1_NW    +=  c1over54* (drho1+(drho1+one)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f1_TE    +=  c1over54* (drho1+(drho1+one)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f1_BW    +=  c1over54* (drho1+(drho1+one)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f1_BE    +=  c1over54* (drho1+(drho1+one)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f1_TW    +=  c1over54* (drho1+(drho1+one)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f1_TN    +=  c1over54* (drho1+(drho1+one)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f1_BS    +=  c1over54* (drho1+(drho1+one)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f1_BN    +=  c1over54* (drho1+(drho1+one)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f1_TS    +=  c1over54* (drho1+(drho1+one)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f1_TNE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f1_BSW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f1_BNE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f1_TSW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f1_TSE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f1_BNW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f1_BSE   +=  c1over216*(drho1+(drho1+one)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f1_TNW   +=  c1over216*(drho1+(drho1+one)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

	  //drho1 = (drho1 + rhoBC[k])/2.f;
	  //drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[dirE   ])[ke   ] = f1_W   ;  
      (D.f[dirW   ])[kw   ] = f1_E   ;	
      (D.f[dirN   ])[kn   ] = f1_S   ;	
      (D.f[dirS   ])[ks   ] = f1_N   ;	
      (D.f[dirT   ])[kt   ] = f1_B   ;	
      (D.f[dirB   ])[kb   ] = f1_T   ;	
      (D.f[dirNE  ])[kne  ] = f1_SW  ;	
      (D.f[dirSW  ])[ksw  ] = f1_NE  ;	
      (D.f[dirSE  ])[kse  ] = f1_NW  ;	
      (D.f[dirNW  ])[knw  ] = f1_SE  ;	
      (D.f[dirTE  ])[kte  ] = f1_BW  ;	
      (D.f[dirBW  ])[kbw  ] = f1_TE  ;	
      (D.f[dirBE  ])[kbe  ] = f1_TW  ;	
      (D.f[dirTW  ])[ktw  ] = f1_BE  ;	
      (D.f[dirTN  ])[ktn  ] = f1_BS  ;	
      (D.f[dirBS  ])[kbs  ] = f1_TN  ;	
      (D.f[dirBN  ])[kbn  ] = f1_TS  ;	
      (D.f[dirTS  ])[kts  ] = f1_BN  ;	
      (D.f[dirZERO])[kzero] = f1_ZERO;	
      (D.f[dirTNE ])[ktne ] = f1_BSW ;	
      (D.f[dirTSW ])[ktsw ] = f1_BNE ;	
      (D.f[dirTSE ])[ktse ] = f1_BNW ;	
      (D.f[dirTNW ])[ktnw ] = f1_BSE ;	
      (D.f[dirBNE ])[kbne ] = f1_TSW ;	
      (D.f[dirBSW ])[kbsw ] = f1_TNE ;	
      (D.f[dirBSE ])[kbse ] = f1_TNW ;	
      (D.f[dirBNW ])[kbnw ] = f1_TSE ;       
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
                                               doubflo* DD, 
                                               unsigned int size_Mat, 
                                               bool evenOrOdd) 
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
      if (evenOrOdd==true)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
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
      doubflo        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                   f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dirE   ])[k1e   ];
      f1_E    = (D.f[dirW   ])[k1w   ];
      f1_S    = (D.f[dirN   ])[k1n   ];
      f1_N    = (D.f[dirS   ])[k1s   ];
      f1_B    = (D.f[dirT   ])[k1t   ];
      f1_T    = (D.f[dirB   ])[k1b   ];
      f1_SW   = (D.f[dirNE  ])[k1ne  ];
      f1_NE   = (D.f[dirSW  ])[k1sw  ];
      f1_NW   = (D.f[dirSE  ])[k1se  ];
      f1_SE   = (D.f[dirNW  ])[k1nw  ];
      f1_BW   = (D.f[dirTE  ])[k1te  ];
      f1_TE   = (D.f[dirBW  ])[k1bw  ];
      f1_TW   = (D.f[dirBE  ])[k1be  ];
      f1_BE   = (D.f[dirTW  ])[k1tw  ];
      f1_BS   = (D.f[dirTN  ])[k1tn  ];
      f1_TN   = (D.f[dirBS  ])[k1bs  ];
      f1_TS   = (D.f[dirBN  ])[k1bn  ];
      f1_BN   = (D.f[dirTS  ])[k1ts  ];
      f1_ZERO = (D.f[dirZERO])[k1zero];
      f1_BSW  = (D.f[dirTNE ])[k1tne ];
      f1_BNE  = (D.f[dirTSW ])[k1tsw ];
      f1_BNW  = (D.f[dirTSE ])[k1tse ];
      f1_BSE  = (D.f[dirTNW ])[k1tnw ];
      f1_TSW  = (D.f[dirBNE ])[k1bne ];
      f1_TNE  = (D.f[dirBSW ])[k1bsw ];
      f1_TNW  = (D.f[dirBSE ])[k1bse ];
      f1_TSE  = (D.f[dirBNW ])[k1bnw ];

      doubflo drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                        f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

      __syncthreads();

      (D.f[dirE   ])[ke   ] = f1_W   -c2over27*drho1;
      (D.f[dirW   ])[kw   ] = f1_E   -c2over27*drho1;
      (D.f[dirN   ])[kn   ] = f1_S   -c2over27*drho1;
      (D.f[dirS   ])[ks   ] = f1_N   -c2over27*drho1;
      (D.f[dirT   ])[kt   ] = f1_B   -c2over27*drho1;
      (D.f[dirB   ])[kb   ] = f1_T   -c2over27*drho1;
      (D.f[dirNE  ])[kne  ] = f1_SW  -c1over54*drho1;
      (D.f[dirSW  ])[ksw  ] = f1_NE  -c1over54*drho1;
      (D.f[dirSE  ])[kse  ] = f1_NW  -c1over54*drho1;
      (D.f[dirNW  ])[knw  ] = f1_SE  -c1over54*drho1;
      (D.f[dirTE  ])[kte  ] = f1_BW  -c1over54*drho1;
      (D.f[dirBW  ])[kbw  ] = f1_TE  -c1over54*drho1;
      (D.f[dirBE  ])[kbe  ] = f1_TW  -c1over54*drho1;
      (D.f[dirTW  ])[ktw  ] = f1_BE  -c1over54*drho1;
      (D.f[dirTN  ])[ktn  ] = f1_BS  -c1over54*drho1;
      (D.f[dirBS  ])[kbs  ] = f1_TN  -c1over54*drho1;
      (D.f[dirBN  ])[kbn  ] = f1_TS  -c1over54*drho1;
      (D.f[dirTS  ])[kts  ] = f1_BN  -c1over54*drho1;
      (D.f[dirZERO])[kzero] = f1_ZERO-c8over27*drho1;
      (D.f[dirTNE ])[ktne ] = f1_BSW -c1over216*drho1;
      (D.f[dirTSW ])[ktsw ] = f1_BNE -c1over216*drho1;
      (D.f[dirTSE ])[ktse ] = f1_BNW -c1over216*drho1;
      (D.f[dirTNW ])[ktnw ] = f1_BSE -c1over216*drho1;
      (D.f[dirBNE ])[kbne ] = f1_TSW -c1over216*drho1;
      (D.f[dirBSW ])[kbsw ] = f1_TNE -c1over216*drho1;
      (D.f[dirBSE ])[kbse ] = f1_TNW -c1over216*drho1;
      (D.f[dirBNW ])[kbnw ] = f1_TSE -c1over216*drho1;       
   }
   __syncthreads();
}          
//////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDevice27(int inx,
                                           int iny,
                                           doubflo* rhoBC,
                                           doubflo* DD, 
                                           int* k_Q, 
                                           doubflo* QQ,
                                           unsigned int sizeQ,
                                           int kQ, 
                                           doubflo om1, 
                                           unsigned int* neighborX,
                                           unsigned int* neighborY,
                                           unsigned int* neighborZ,
                                           unsigned int size_Mat, 
                                           bool evenOrOdd)
{
   Distributions27 D;
   if (evenOrOdd==true)
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
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
   {
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *sizeQ];
      q_dirW   = &QQ[dirW   *sizeQ];
      q_dirN   = &QQ[dirN   *sizeQ];
      q_dirS   = &QQ[dirS   *sizeQ];
      q_dirT   = &QQ[dirT   *sizeQ];
      q_dirB   = &QQ[dirB   *sizeQ];
      q_dirNE  = &QQ[dirNE  *sizeQ];
      q_dirSW  = &QQ[dirSW  *sizeQ];
      q_dirSE  = &QQ[dirSE  *sizeQ];
      q_dirNW  = &QQ[dirNW  *sizeQ];
      q_dirTE  = &QQ[dirTE  *sizeQ];
      q_dirBW  = &QQ[dirBW  *sizeQ];
      q_dirBE  = &QQ[dirBE  *sizeQ];
      q_dirTW  = &QQ[dirTW  *sizeQ];
      q_dirTN  = &QQ[dirTN  *sizeQ];
      q_dirBS  = &QQ[dirBS  *sizeQ];
      q_dirBN  = &QQ[dirBN  *sizeQ];
      q_dirTS  = &QQ[dirTS  *sizeQ];
      q_dirTNE = &QQ[dirTNE *sizeQ];
      q_dirTSW = &QQ[dirTSW *sizeQ];
      q_dirTSE = &QQ[dirTSE *sizeQ];
      q_dirTNW = &QQ[dirTNW *sizeQ];
      q_dirBNE = &QQ[dirBNE *sizeQ];
      q_dirBSW = &QQ[dirBSW *sizeQ];
      q_dirBSE = &QQ[dirBSE *sizeQ];
      q_dirBNW = &QQ[dirBNW *sizeQ];
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      ////////////////////////////////////////////////////////////////////////////////
      doubflo q, vx1, vx2, vx3, drho;
      vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
                  (f_E - f_W); 


      vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
                  ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
                  (f_N - f_S); 

      vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
                  (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
                  (f_T - f_B); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      drho = rhoBC[k];
      //deltaRho = (rhoBC[k] + one) / (deltaRho + one);
      ////////////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirW])[kw]=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
         //(D.f[dirE])[ke]=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirE])[ke]=c2over27* (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq); 
         //(D.f[dirW])[kw]=c2over27* (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirS])[ks]=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
         //(D.f[dirN])[kn]=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirN])[kn]=c2over27* (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
         //(D.f[dirS])[ks]=c2over27* (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirB])[kb]=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
         //(D.f[dirT])[kt]=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirT])[kt]=c2over27* (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq); 
         //(D.f[dirB])[kb]=c2over27* (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirSW])[ksw]=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
         //(D.f[dirNE])[kne]=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirNE])[kne]=c1over54* (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
         //(D.f[dirSW])[ksw]=c1over54* (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirNW])[knw]=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
         //(D.f[dirSE])[kse]=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirSE])[kse]=c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
         //(D.f[dirNW])[knw]=c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBW])[kbw]=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
         //(D.f[dirTE])[kte]=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTE])[kte]=c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
         //(D.f[dirBW])[kbw]=c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTW])[ktw]=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
         //(D.f[dirBE])[kbe]=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBE])[kbe]=c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
         //(D.f[dirTW])[ktw]=c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBS])[kbs]=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
         //(D.f[dirTN])[ktn]=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTN])[ktn]=c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
         //(D.f[dirBS])[kbs]=c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTS])[kts]=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
         //(D.f[dirBN])[kbn]=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBN])[kbn]=c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
         //(D.f[dirTS])[kts]=c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBSW])[kbsw]=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
         //(D.f[dirTNE])[ktne]=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTNE])[ktne]=c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
         //(D.f[dirBSW])[kbsw]=c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTSW])[ktsw]=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
         //(D.f[dirBNE])[kbne]=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBNE])[kbne]=c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
         //(D.f[dirTSW])[ktsw]=c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBNW])[kbnw]=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
         //(D.f[dirTSE])[ktse]=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTSE])[ktse]=c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
         //(D.f[dirBNW])[kbnw]=c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTNW])[ktnw]=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
         //(D.f[dirBSE])[kbse]=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBSE])[kbse]=c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
         //(D.f[dirTNW])[ktnw]=c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceAntiBB27(   doubflo* rhoBC,
												   doubflo* vx,
												   doubflo* vy,
												   doubflo* vz,
												   doubflo* DD, 
												   int* k_Q, 
												   doubflo* QQ,
												   int kQ, 
												   doubflo om1, 
												   unsigned int* neighborX,
												   unsigned int* neighborY,
												   unsigned int* neighborZ,
												   unsigned int size_Mat, 
												   bool evenOrOdd)
{
   Distributions27 D;
   if (evenOrOdd==true)
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
      D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      D.f[dirBNW ] = &DD[dirTSE *size_Mat];
   }
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
   {
      doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
         *q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
         *q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
         *q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
         *q_dirBSE, *q_dirBNW; 
      q_dirE   = &QQ[dirE   *kQ];
      q_dirW   = &QQ[dirW   *kQ];
      q_dirN   = &QQ[dirN   *kQ];
      q_dirS   = &QQ[dirS   *kQ];
      q_dirT   = &QQ[dirT   *kQ];
      q_dirB   = &QQ[dirB   *kQ];
      q_dirNE  = &QQ[dirNE  *kQ];
      q_dirSW  = &QQ[dirSW  *kQ];
      q_dirSE  = &QQ[dirSE  *kQ];
      q_dirNW  = &QQ[dirNW  *kQ];
      q_dirTE  = &QQ[dirTE  *kQ];
      q_dirBW  = &QQ[dirBW  *kQ];
      q_dirBE  = &QQ[dirBE  *kQ];
      q_dirTW  = &QQ[dirTW  *kQ];
      q_dirTN  = &QQ[dirTN  *kQ];
      q_dirBS  = &QQ[dirBS  *kQ];
      q_dirBN  = &QQ[dirBN  *kQ];
      q_dirTS  = &QQ[dirTS  *kQ];
      q_dirTNE = &QQ[dirTNE *kQ];
      q_dirTSW = &QQ[dirTSW *kQ];
      q_dirTSE = &QQ[dirTSE *kQ];
      q_dirTNW = &QQ[dirTNW *kQ];
      q_dirBNE = &QQ[dirBNE *kQ];
      q_dirBSW = &QQ[dirBSW *kQ];
      q_dirBSE = &QQ[dirBSE *kQ];
      q_dirBNW = &QQ[dirBNW *kQ];
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
      doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
         f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW, f_ZERO;

      f_W    = (D.f[dirE   ])[ke   ];
      f_E    = (D.f[dirW   ])[kw   ];
      f_S    = (D.f[dirN   ])[kn   ];
      f_N    = (D.f[dirS   ])[ks   ];
      f_B    = (D.f[dirT   ])[kt   ];
      f_T    = (D.f[dirB   ])[kb   ];
      f_SW   = (D.f[dirNE  ])[kne  ];
      f_NE   = (D.f[dirSW  ])[ksw  ];
      f_NW   = (D.f[dirSE  ])[kse  ];
      f_SE   = (D.f[dirNW  ])[knw  ];
      f_BW   = (D.f[dirTE  ])[kte  ];
      f_TE   = (D.f[dirBW  ])[kbw  ];
      f_TW   = (D.f[dirBE  ])[kbe  ];
      f_BE   = (D.f[dirTW  ])[ktw  ];
      f_BS   = (D.f[dirTN  ])[ktn  ];
      f_TN   = (D.f[dirBS  ])[kbs  ];
      f_TS   = (D.f[dirBN  ])[kbn  ];
      f_BN   = (D.f[dirTS  ])[kts  ];
      f_BSW  = (D.f[dirTNE ])[ktne ];
      f_BNE  = (D.f[dirTSW ])[ktsw ];
      f_BNW  = (D.f[dirTSE ])[ktse ];
      f_BSE  = (D.f[dirTNW ])[ktnw ];
      f_TSW  = (D.f[dirBNE ])[kbne ];
      f_TNE  = (D.f[dirBSW ])[kbsw ];
      f_TNW  = (D.f[dirBSE ])[kbse ];
      f_TSE  = (D.f[dirBNW ])[kbnw ];
      f_ZERO = (D.f[dirZERO])[kzero];
      ////////////////////////////////////////////////////////////////////////////////
      //doubflo vx1, vx2, vx3, drho;
      //vx1    =  ((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //            ((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
      //            (f_E - f_W); 


      //vx2    =   (-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
      //            ((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
      //            (f_N - f_S); 

      //vx3    =   ((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
      //            (-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
      //            (f_T - f_B); 

      //doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      doubflo drho    = f_ZERO+f_E+f_W+f_N+f_S+f_T+f_B+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+
						f_TNE+f_TSW+f_TSE+f_TNW+f_BNE+f_BSW+f_BSE+f_BNW;
      drho = drho - rhoBC[k];
	  drho *= 0.01f;
      ////////////////////////////////////////////////////////////////////////////////
	  doubflo q;
      //deltaRho = (rhoBC[k] + one) / (deltaRho + one);
      ////////////////////////////////////////////////////////////////////////////////
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      q = q_dirE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirW])[kw]=f_W-c2over27*drho; 
      }

      q = q_dirW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirE])[ke]=f_E-c2over27*drho;
      }

      q = q_dirN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirS])[ks]=f_S-c2over27*drho; 
      }

      q = q_dirS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirN])[kn]=f_N-c2over27*drho; 
      }

      q = q_dirT[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirB])[kb]=f_B-c2over27*drho; 
      }

      q = q_dirB[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirT])[kt]=f_T-c2over27*drho;
      }

      q = q_dirNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirSW])[ksw]=f_SW-c1over54*drho;
      }

      q = q_dirSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirNE])[kne]=f_NE-c1over54*drho;
      }

      q = q_dirSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirNW])[knw]=f_NW-c1over54*drho;
      }

      q = q_dirNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirSE])[kse]=f_SE-c1over54*drho;
      }

      q = q_dirTE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBW])[kbw]=f_BW-c1over54*drho; 
      }

      q = q_dirBW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTE])[kte]=f_TE-c1over54*drho; 
      }

      q = q_dirBE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTW])[ktw]=f_TW-c1over54*drho;
      }

      q = q_dirTW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBE])[kbe]=f_BE-c1over54*drho;
      }

      q = q_dirTN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBS])[kbs]=f_BS-c1over54*drho;
      }

      q = q_dirBS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTN])[ktn]=f_TN-c1over54*drho;
      }

      q = q_dirBN[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTS])[kts]=f_TS-c1over54*drho;
      }

      q = q_dirTS[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBN])[kbn]=f_BN-c1over54*drho;
      }

      q = q_dirTNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBSW])[kbsw]=f_BSW-c1over216*drho;
      }

      q = q_dirBSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTNE])[ktne]=f_TNE-c1over216*drho;
      }

      q = q_dirBNE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTSW])[ktsw]=f_TSW-c1over216*drho;
      }

      q = q_dirTSW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBNE])[kbne]=f_BNE-c1over216*drho;
      }

      q = q_dirTSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBNW])[kbnw]=f_BNW-c1over216*drho;
      }

      q = q_dirBNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTSE])[ktse]=f_TSE-c1over216*drho;
      }

      q = q_dirBSE[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirTNW])[ktnw]=f_TNW-c1over216*drho;
      }

      q = q_dirTNW[k];
      if (q>=zero && q<=one)
      {
         (D.f[dirBSE])[kbse]=f_BSE-c1over216*drho;
      }
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceFixBackflow27( doubflo* rhoBC,
                                                      doubflo* DD, 
                                                      int* k_Q, 
                                                      int kQ, 
                                                      doubflo om1, 
                                                      unsigned int* neighborX,
                                                      unsigned int* neighborY,
                                                      unsigned int* neighborZ,
                                                      unsigned int size_Mat, 
                                                      bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      doubflo deltaRho;
      ////////////////////////////////////////////////////////////////////////////////
      deltaRho = rhoBC[k];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         (D.f[dirW])[kw]       = c2over27  * deltaRho;
         (D.f[dirE])[ke]       = c2over27  * deltaRho;
         (D.f[dirS])[ks]       = c2over27  * deltaRho;
         (D.f[dirN])[kn]       = c2over27  * deltaRho;
         (D.f[dirB])[kb]       = c2over27  * deltaRho;
         (D.f[dirT])[kt]       = c2over27  * deltaRho;
         (D.f[dirSW])[ksw]     = c1over54  * deltaRho;
         (D.f[dirNE])[kne]     = c1over54  * deltaRho;
         (D.f[dirNW])[knw]     = c1over54  * deltaRho;
         (D.f[dirSE])[kse]     = c1over54  * deltaRho;
         (D.f[dirBW])[kbw]     = c1over54  * deltaRho;
         (D.f[dirTE])[kte]     = c1over54  * deltaRho;
         (D.f[dirTW])[ktw]     = c1over54  * deltaRho;
         (D.f[dirBE])[kbe]     = c1over54  * deltaRho;
         (D.f[dirBS])[kbs]     = c1over54  * deltaRho;
         (D.f[dirTN])[ktn]     = c1over54  * deltaRho;
         (D.f[dirTS])[kts]     = c1over54  * deltaRho;
         (D.f[dirBN])[kbn]     = c1over54  * deltaRho;
         (D.f[dirBSW])[kbsw]   = c1over216 * deltaRho;
         (D.f[dirTNE])[ktne]   = c1over216 * deltaRho;
         (D.f[dirTSW])[ktsw]   = c1over216 * deltaRho;
         (D.f[dirBNE])[kbne]   = c1over216 * deltaRho;
         (D.f[dirBNW])[kbnw]   = c1over216 * deltaRho;
         (D.f[dirTSE])[ktse]   = c1over216 * deltaRho;
         (D.f[dirTNW])[ktnw]   = c1over216 * deltaRho;
         (D.f[dirBSE])[kbse]   = c1over216 * deltaRho;
         (D.f[dirZERO])[kzero] = c8over27  * deltaRho;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceDirDepBot27(  doubflo* rhoBC,
                                                     doubflo* DD, 
                                                     int* k_Q, 
                                                     int kQ, 
                                                     doubflo om1, 
                                                     unsigned int* neighborX,
                                                     unsigned int* neighborY,
                                                     unsigned int* neighborZ,
                                                     unsigned int size_Mat, 
                                                     bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      doubflo rho;
      ////////////////////////////////////////////////////////////////////////////////
      rho = rhoBC[k];
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 D;
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo f_E,f_W,f_N,f_S,f_T,f_NE,f_SW,f_SE,f_NW,f_TE,f_TW,f_TN,f_TS,f_ZERO,f_TNE,f_TSW,f_TSE,f_TNW;//,
            //f_B,f_BW,f_BE,f_BS,f_BN,f_BSW,f_BNE,f_BNW,f_BSE;

      f_E    = (D.f[dirE   ])[ke   ];
      f_W    = (D.f[dirW   ])[kw   ];
      f_N    = (D.f[dirN   ])[kn   ];
      f_S    = (D.f[dirS   ])[ks   ];
      f_T    = (D.f[dirT   ])[kt   ];
      f_NE   = (D.f[dirNE  ])[kne  ];
      f_SW   = (D.f[dirSW  ])[ksw  ];
      f_SE   = (D.f[dirSE  ])[kse  ];
      f_NW   = (D.f[dirNW  ])[knw  ];
      f_TE   = (D.f[dirTE  ])[kte  ];
      f_TW   = (D.f[dirTW  ])[ktw  ];
      f_TN   = (D.f[dirTN  ])[ktn  ];
      f_TS   = (D.f[dirTS  ])[kts  ];
      f_ZERO = (D.f[dirZERO])[kzero];
      f_TNE  = (D.f[dirTNE ])[ktne ];
      f_TSW  = (D.f[dirTSW ])[ktsw ];
      f_TSE  = (D.f[dirTSE ])[ktse ];
      f_TNW  = (D.f[dirTNW ])[ktnw ];
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

      //doubflo drho   =    f_ZERO+f_E+f_W+f_N+f_S+f_T+f_B+f_NE+f_SW+f_SE+f_NW+f_TE+f_BW+f_BE+f_TW+f_TN+f_BS+f_BN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW+f_BNE+f_BSW+f_BSE+f_BNW;
      //doubflo vx1     =  (f_E -f_W +f_NE-f_SW+f_SE-f_NW+f_TE-f_BW+f_BE-f_TW+ f_TNE-f_TSW+f_TSE-f_TNW+ f_BNE-f_BSW+f_BSE-f_BNW);
      //doubflo vx2     =  (f_N -f_S +f_NE-f_SW-f_SE+f_NW+f_TN-f_BS+f_BN-f_TS+ f_TNE-f_TSW-f_TSE+f_TNW+ f_BNE-f_BSW-f_BSE+f_BNW);
      //doubflo vx3     =  (f_T -f_B +f_TE-f_BW-f_BE+f_TW+f_TN-f_BS-f_BN+f_TS+ f_TNE+f_TSW+f_TSE+f_TNW- f_BNE-f_BSW-f_BSE-f_BNW);

      //doubflo cusq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

      //(D.f[dirZERO])[kzero] = c8over27*  (drho-cusq);
      //(D.f[dirE])[ke]    = c2over27*  (drho+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq);
      //(D.f[dirW])[kw]    = c2over27*  (drho+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq);
      //(D.f[dirN])[kn]     = c2over27*  (drho+three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq);
      //(D.f[dirS])[ks]    = c2over27*  (drho+three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq);
      //(D.f[dirT])[kt]    = c2over27*  (drho+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq);
      //(D.f[dirB])[kb]    = c2over27*  (drho+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq);
      //(D.f[dirNE])[kne]   = c1over54*  (drho+three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq);
      //(D.f[dirSW])[ksw]   = c1over54*  (drho+three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq);
      //(D.f[dirSE])[kse]   =  c1over54* (drho+three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq);
      //(D.f[dirNW])[knw]   =  c1over54* (drho+three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq);
      //(D.f[dirTE])[kte]   =  c1over54* (drho+three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq);
      //(D.f[dirBW])[kbw]   =  c1over54* (drho+three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq);
      //(D.f[dirBE])[kbe]   =  c1over54* (drho+three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq);
      //(D.f[dirTW])[ktw]   =  c1over54* (drho+three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq);
      //(D.f[dirTN])[ktn]   =  c1over54* (drho+three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq);
      //(D.f[dirBS])[kbs]   =  c1over54* (drho+three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq);
      //(D.f[dirBN])[kbn]   =  c1over54* (drho+three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq);
      //(D.f[dirTS])[kts]   =  c1over54* (drho+three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq);
      //(D.f[dirTNE])[ktne]  =  c1over216*(drho+three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq);
      //(D.f[dirBSW])[kbsw]  =  c1over216*(drho+three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq);
      //(D.f[dirBNE])[kbne]  =  c1over216*(drho+three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq);
      //(D.f[dirTSW])[ktsw]  =  c1over216*(drho+three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq);
      //(D.f[dirTSE])[ktse]  =  c1over216*(drho+three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq);
      //(D.f[dirBNW])[kbnw]  =  c1over216*(drho+three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq);
      //(D.f[dirBSE])[kbse]  =  c1over216*(drho+three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq);
      //(D.f[dirTNW])[ktnw]  =  c1over216*(drho+three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq);
      doubflo drho   =    f_ZERO+f_E+f_W+f_N+f_S+f_T+f_NE+f_SW+f_SE+f_NW+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      doubflo dTop   =    f_T+f_TE+f_TW+f_TN+f_TS+f_TNE+f_TSW+f_TSE+f_TNW;
      (D.f[dirB])[kb]     = (f_T+c2over27)*(rho-drho+one/six)/(dTop+one/six)-c2over27;
      (D.f[dirBW])[kbw]   = (f_TW+c1over54)*(rho-drho+one/six)/(dTop+one/six)-c1over54;
      (D.f[dirBE])[kbe]   = (f_TE+c1over54)*(rho-drho+one/six)/(dTop+one/six)-c1over54;
      (D.f[dirBS])[kbs]   = (f_TS+c1over54)*(rho-drho+one/six)/(dTop+one/six)-c1over54;
      (D.f[dirBN])[kbn]   = (f_TN+c1over54)*(rho-drho+one/six)/(dTop+one/six)-c1over54;
      (D.f[dirBSW])[kbsw] = (f_TSW+c1over216)*(rho-drho+one/six)/(dTop+one/six)-c1over216;
      (D.f[dirBNE])[kbne] = (f_TNE+c1over216)*(rho-drho+one/six)/(dTop+one/six)-c1over216;
      (D.f[dirBNW])[kbnw] = (f_TNW+c1over216)*(rho-drho+one/six)/(dTop+one/six)-c1over216;
      (D.f[dirBSE])[kbse] = (f_TSE+c1over216)*(rho-drho+one/six)/(dTop+one/six)-c1over216;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressNoRhoDevice27(  doubflo* rhoBC,
												 doubflo* DD, 
												 int* k_Q, 
												 int* k_N, 
												 int kQ, 
												 doubflo om1, 
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
												 bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==true)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo f1_E    = (D.f[dirE   ])[k1e   ];
      doubflo f1_W    = (D.f[dirW   ])[k1w   ];
      doubflo f1_N    = (D.f[dirN   ])[k1n   ];
      doubflo f1_S    = (D.f[dirS   ])[k1s   ];
      doubflo f1_T    = (D.f[dirT   ])[k1t   ];
      doubflo f1_B    = (D.f[dirB   ])[k1b   ];
      doubflo f1_NE   = (D.f[dirNE  ])[k1ne  ];
      doubflo f1_SW   = (D.f[dirSW  ])[k1sw  ];
      doubflo f1_SE   = (D.f[dirSE  ])[k1se  ];
      doubflo f1_NW   = (D.f[dirNW  ])[k1nw  ];
      doubflo f1_TE   = (D.f[dirTE  ])[k1te  ];
      doubflo f1_BW   = (D.f[dirBW  ])[k1bw  ];
      doubflo f1_BE   = (D.f[dirBE  ])[k1be  ];
      doubflo f1_TW   = (D.f[dirTW  ])[k1tw  ];
      doubflo f1_TN   = (D.f[dirTN  ])[k1tn  ];
      doubflo f1_BS   = (D.f[dirBS  ])[k1bs  ];
      doubflo f1_BN   = (D.f[dirBN  ])[k1bn  ];
      doubflo f1_TS   = (D.f[dirTS  ])[k1ts  ];
      doubflo f1_ZERO = (D.f[dirZERO])[k1zero];
      doubflo f1_TNE  = (D.f[dirTNE ])[k1tne ];
      doubflo f1_TSW  = (D.f[dirTSW ])[k1tsw ];
      doubflo f1_TSE  = (D.f[dirTSE ])[k1tse ];
      doubflo f1_TNW  = (D.f[dirTNW ])[k1tnw ];
      doubflo f1_BNE  = (D.f[dirBNE ])[k1bne ];
      doubflo f1_BSW  = (D.f[dirBSW ])[k1bsw ];
      doubflo f1_BSE  = (D.f[dirBSE ])[k1bse ];
      doubflo f1_BNW  = (D.f[dirBNW ])[k1bnw ];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo f_E    = (D.f[dirE   ])[ke   ];
      doubflo f_W    = (D.f[dirW   ])[kw   ];
      doubflo f_N    = (D.f[dirN   ])[kn   ];
      doubflo f_S    = (D.f[dirS   ])[ks   ];
      doubflo f_T    = (D.f[dirT   ])[kt   ];
      doubflo f_B    = (D.f[dirB   ])[kb   ];
      doubflo f_NE   = (D.f[dirNE  ])[kne  ];
      doubflo f_SW   = (D.f[dirSW  ])[ksw  ];
      doubflo f_SE   = (D.f[dirSE  ])[kse  ];
      doubflo f_NW   = (D.f[dirNW  ])[knw  ];
      doubflo f_TE   = (D.f[dirTE  ])[kte  ];
      doubflo f_BW   = (D.f[dirBW  ])[kbw  ];
      doubflo f_BE   = (D.f[dirBE  ])[kbe  ];
      doubflo f_TW   = (D.f[dirTW  ])[ktw  ];
      doubflo f_TN   = (D.f[dirTN  ])[ktn  ];
      doubflo f_BS   = (D.f[dirBS  ])[kbs  ];
      doubflo f_BN   = (D.f[dirBN  ])[kbn  ];
      doubflo f_TS   = (D.f[dirTS  ])[kts  ];
      doubflo f_ZERO = (D.f[dirZERO])[kzero];
      doubflo f_TNE  = (D.f[dirTNE ])[ktne ];
      doubflo f_TSW  = (D.f[dirTSW ])[ktsw ];
      doubflo f_TSE  = (D.f[dirTSE ])[ktse ];
      doubflo f_TNW  = (D.f[dirTNW ])[ktnw ];
      doubflo f_BNE  = (D.f[dirBNE ])[kbne ];
      doubflo f_BSW  = (D.f[dirBSW ])[kbsw ];
      doubflo f_BSE  = (D.f[dirBSE ])[kbse ];
      doubflo f_BNW  = (D.f[dirBNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////

      //doubflo vx1, vx2, vx3, drho;
      doubflo vx1, vx2, vx3, drho, drho1;
      //////////////////////////////////////////////////////////////////////////
	  //Dichte
      drho1  =  f1_TSE + f1_TNW + f1_TNE + f1_TSW + f1_BSE + f1_BNW + f1_BNE + f1_BSW +
                f1_BN + f1_TS + f1_TN + f1_BS + f1_BE + f1_TW + f1_TE + f1_BW + f1_SE + f1_NW + f1_NE + f1_SW + 
                f1_T + f1_B + f1_N + f1_S + f1_E + f1_W + ((D.f[dirZERO])[k1zero]); 
      drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
                f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
                f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 
      
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


      //doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

   //   //////////////////////////////////////////////////////////////////////////
	  ////doubflo omega = om1;
   //   doubflo cusq  = c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
   //   //////////////////////////////////////////////////////////////////////////
	  ////Täst MK
	  ////if(vx1 < zero) vx1 = zero;
   //   //////////////////////////////////////////////////////////////////////////
   //   doubflo fZERO = c8over27*  (drho1-(one + drho1)*(cusq))                                                           ;
   //   doubflo fE    = c2over27*  (drho1+(one + drho1)*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
   //   doubflo fW    = c2over27*  (drho1+(one + drho1)*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
   //   doubflo fN    = c2over27*  (drho1+(one + drho1)*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
   //   doubflo fS    = c2over27*  (drho1+(one + drho1)*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
   //   doubflo fT    = c2over27*  (drho1+(one + drho1)*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
   //   doubflo fB    = c2over27*  (drho1+(one + drho1)*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
   //   doubflo fNE   = c1over54*  (drho1+(one + drho1)*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
   //   doubflo fSW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
   //   doubflo fSE   = c1over54*  (drho1+(one + drho1)*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
   //   doubflo fNW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
   //   doubflo fTE   = c1over54*  (drho1+(one + drho1)*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
   //   doubflo fBW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
   //   doubflo fBE   = c1over54*  (drho1+(one + drho1)*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
   //   doubflo fTW   = c1over54*  (drho1+(one + drho1)*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
   //   doubflo fTN   = c1over54*  (drho1+(one + drho1)*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
   //   doubflo fBS   = c1over54*  (drho1+(one + drho1)*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
   //   doubflo fBN   = c1over54*  (drho1+(one + drho1)*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
   //   doubflo fTS   = c1over54*  (drho1+(one + drho1)*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
   //   doubflo fTNE  = c1over216* (drho1+(one + drho1)*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
   //   doubflo fBSW  = c1over216* (drho1+(one + drho1)*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
   //   doubflo fBNE  = c1over216* (drho1+(one + drho1)*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
   //   doubflo fTSW  = c1over216* (drho1+(one + drho1)*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
   //   doubflo fTSE  = c1over216* (drho1+(one + drho1)*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
   //   doubflo fBNW  = c1over216* (drho1+(one + drho1)*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
   //   doubflo fBSE  = c1over216* (drho1+(one + drho1)*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
   //   doubflo fTNW  = c1over216* (drho1+(one + drho1)*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

	  doubflo cs = one / sqrtf(three);
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //no velocity
	  //////////////////////////////////////////
      f_E    = f1_E   * cs + (one - cs) * f_E   ;
      f_W    = f1_W   * cs + (one - cs) * f_W   ;
      f_N    = f1_N   * cs + (one - cs) * f_N   ;
      f_S    = f1_S   * cs + (one - cs) * f_S   ;
      f_T    = f1_T   * cs + (one - cs) * f_T   ;
      f_B    = f1_B   * cs + (one - cs) * f_B   ;
      f_NE   = f1_NE  * cs + (one - cs) * f_NE  ;
      f_SW   = f1_SW  * cs + (one - cs) * f_SW  ;
      f_SE   = f1_SE  * cs + (one - cs) * f_SE  ;
      f_NW   = f1_NW  * cs + (one - cs) * f_NW  ;
      f_TE   = f1_TE  * cs + (one - cs) * f_TE  ;
      f_BW   = f1_BW  * cs + (one - cs) * f_BW  ;
      f_BE   = f1_BE  * cs + (one - cs) * f_BE  ;
      f_TW   = f1_TW  * cs + (one - cs) * f_TW  ;
      f_TN   = f1_TN  * cs + (one - cs) * f_TN  ;
      f_BS   = f1_BS  * cs + (one - cs) * f_BS  ;
      f_BN   = f1_BN  * cs + (one - cs) * f_BN  ;
      f_TS   = f1_TS  * cs + (one - cs) * f_TS  ;
      f_TNE  = f1_TNE * cs + (one - cs) * f_TNE ;
      f_TSW  = f1_TSW * cs + (one - cs) * f_TSW ;
      f_TSE  = f1_TSE * cs + (one - cs) * f_TSE ;
      f_TNW  = f1_TNW * cs + (one - cs) * f_TNW ;
      f_BNE  = f1_BNE * cs + (one - cs) * f_BNE ;
      f_BSW  = f1_BSW * cs + (one - cs) * f_BSW ;
      f_BSE  = f1_BSE * cs + (one - cs) * f_BSE ;
      f_BNW  = f1_BNW * cs + (one - cs) * f_BNW ;
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //with velocity
	  //if(true){//vx1 >= zero){
		 // doubflo csMvx = one / sqrtf(three) - vx1;
		 // //doubflo csMvy = one / sqrtf(three) - vx2;
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
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////
      //__syncthreads();
	  // -X
	  //(D.f[dirE   ])[ke   ] = f_E   ;
	  //(D.f[dirSE  ])[kse  ] = f_SE  ;
	  //(D.f[dirNE  ])[kne  ] = f_NE  ;
	  //(D.f[dirBE  ])[kbe  ] = f_BE  ;
	  //(D.f[dirTE  ])[kte  ] = f_TE  ;
	  //(D.f[dirTSE ])[ktse ] = f_TSE ;
	  //(D.f[dirTNE ])[ktne ] = f_TNE ;
	  //(D.f[dirBSE ])[kbse ] = f_BSE ;
	  //(D.f[dirBNE ])[kbne ] = f_BNE ;     
	  // X
	  //(D.f[dirW   ])[kw   ] = f_W   ;
	  //(D.f[dirSW  ])[ksw  ] = f_SW  ;
	  //(D.f[dirNW  ])[knw  ] = f_NW  ;
	  //(D.f[dirBW  ])[kbw  ] = f_BW  ;
	  //(D.f[dirTW  ])[ktw  ] = f_TW  ;
	  //(D.f[dirTSW ])[ktsw ] = f_TSW ;
	  //(D.f[dirTNW ])[ktnw ] = f_TNW ;
	  //(D.f[dirBSW ])[kbsw ] = f_BSW ;
	  //(D.f[dirBNW ])[kbnw ] = f_BNW ;     
	  // Y
	  //(D.f[dirS   ])[ks   ] = f_S   ;
	  //(D.f[dirSE  ])[kse  ] = f_SE  ;
	  //(D.f[dirSW  ])[ksw  ] = f_SW  ;
	  //(D.f[dirTS  ])[kts  ] = f_TS  ;
	  //(D.f[dirBS  ])[kbs  ] = f_BS  ;
	  //(D.f[dirTSE ])[ktse ] = f_TSE ;
	  //(D.f[dirTSW ])[ktsw ] = f_TSW ;
	  //(D.f[dirBSE ])[kbse ] = f_BSE ;
	  //(D.f[dirBSW ])[kbsw ] = f_BSW ;     
	  // Z
	  (D.f[dirB   ])[kb   ] = f_B   ;
	  (D.f[dirBE  ])[kbe  ] = f_BE  ;
	  (D.f[dirBW  ])[kbw  ] = f_BW  ;
	  (D.f[dirBN  ])[kbn  ] = f_BN  ;
	  (D.f[dirBS  ])[kbs  ] = f_BS  ;
	  (D.f[dirBNE ])[kbne ] = f_BNE ;
	  (D.f[dirBNW ])[kbnw ] = f_BNW ;
	  (D.f[dirBSE ])[kbse ] = f_BSE ;
	  (D.f[dirBSW ])[kbsw ] = f_BSW ;     
      //////////////////////////////////////////////////////////////////////////
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceOld27(doubflo* rhoBC,
                                             doubflo* DD, 
                                             int* k_Q, 
                                             int* k_N, 
                                             int kQ, 
                                             doubflo om1, 
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int size_Mat, 
                                             bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
                     f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dirE   ])[k1e   ];
      f1_E    = (D.f[dirW   ])[k1w   ];
      f1_S    = (D.f[dirN   ])[k1n   ];
      f1_N    = (D.f[dirS   ])[k1s   ];
      f1_B    = (D.f[dirT   ])[k1t   ];
      f1_T    = (D.f[dirB   ])[k1b   ];
      f1_SW   = (D.f[dirNE  ])[k1ne  ];
      f1_NE   = (D.f[dirSW  ])[k1sw  ];
      f1_NW   = (D.f[dirSE  ])[k1se  ];
      f1_SE   = (D.f[dirNW  ])[k1nw  ];
      f1_BW   = (D.f[dirTE  ])[k1te  ];
      f1_TE   = (D.f[dirBW  ])[k1bw  ];
      f1_TW   = (D.f[dirBE  ])[k1be  ];
      f1_BE   = (D.f[dirTW  ])[k1tw  ];
      f1_BS   = (D.f[dirTN  ])[k1tn  ];
      f1_TN   = (D.f[dirBS  ])[k1bs  ];
      f1_TS   = (D.f[dirBN  ])[k1bn  ];
      f1_BN   = (D.f[dirTS  ])[k1ts  ];
      f1_ZERO = (D.f[dirZERO])[k1zero];
      f1_BSW  = (D.f[dirTNE ])[k1tne ];
      f1_BNE  = (D.f[dirTSW ])[k1tsw ];
      f1_BNW  = (D.f[dirTSE ])[k1tse ];
      f1_BSE  = (D.f[dirTNW ])[k1tnw ];
      f1_TSW  = (D.f[dirBNE ])[k1bne ];
      f1_TNE  = (D.f[dirBSW ])[k1bsw ];
      f1_TNW  = (D.f[dirBSE ])[k1bse ];
      f1_TSE  = (D.f[dirBNW ])[k1bnw ];

      //////////////////////////////////////////////////////////////////////////
      doubflo drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
                          f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

	  //drho1 = (drho1 + rhoBC[k])/2.f;
	  drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

      __syncthreads();

      (D.f[dirE   ])[ke   ] = f1_W   -c2over27*drho1;   //  c1o100;  // zero;  //
      (D.f[dirW   ])[kw   ] = f1_E   -c2over27*drho1;	//  c1o100;  // zero;  //
      (D.f[dirN   ])[kn   ] = f1_S   -c2over27*drho1;	//  c1o100;  // zero;  //
      (D.f[dirS   ])[ks   ] = f1_N   -c2over27*drho1;	//  c1o100;  // zero;  //
      (D.f[dirT   ])[kt   ] = f1_B   -c2over27*drho1;	//  c1o100;  // zero;  //
      (D.f[dirB   ])[kb   ] = f1_T   -c2over27*drho1;	//  c1o100;  // zero;  //
      (D.f[dirNE  ])[kne  ] = f1_SW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirSW  ])[ksw  ] = f1_NE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirSE  ])[kse  ] = f1_NW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirNW  ])[knw  ] = f1_SE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTE  ])[kte  ] = f1_BW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBW  ])[kbw  ] = f1_TE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBE  ])[kbe  ] = f1_TW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTW  ])[ktw  ] = f1_BE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTN  ])[ktn  ] = f1_BS  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBS  ])[kbs  ] = f1_TN  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBN  ])[kbn  ] = f1_TS  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTS  ])[kts  ] = f1_BN  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirZERO])[kzero] = f1_ZERO-c8over27*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTNE ])[ktne ] = f1_BSW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTSW ])[ktsw ] = f1_BNE -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTSE ])[ktse ] = f1_BNW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTNW ])[ktnw ] = f1_BSE -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBNE ])[kbne ] = f1_TSW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBSW ])[kbsw ] = f1_TNE -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBSE ])[kbse ] = f1_TNW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBNW ])[kbnw ] = f1_TSE -c1over216*drho1;  //  c1o100;  // zero;  //      
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceEQZ27(doubflo* rhoBC,
                                             doubflo* DD, 
                                             int* k_Q, 
                                             int* k_N,
											 doubflo* kTestRE,
                                             int kQ, 
                                             doubflo om1, 
                                             unsigned int* neighborX,
                                             unsigned int* neighborY,
                                             unsigned int* neighborZ,
                                             unsigned int size_Mat, 
                                             bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==true)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      ////////////////////////////////////////////////////////////////////////////////
      Distributions27 kDistTest;
         kDistTest.f[dirE   ] = &kTestRE[dirE   *kQ];
         kDistTest.f[dirW   ] = &kTestRE[dirW   *kQ];
         kDistTest.f[dirN   ] = &kTestRE[dirN   *kQ];
         kDistTest.f[dirS   ] = &kTestRE[dirS   *kQ];
         kDistTest.f[dirT   ] = &kTestRE[dirT   *kQ];
         kDistTest.f[dirB   ] = &kTestRE[dirB   *kQ];
         kDistTest.f[dirNE  ] = &kTestRE[dirNE  *kQ];
         kDistTest.f[dirSW  ] = &kTestRE[dirSW  *kQ];
         kDistTest.f[dirSE  ] = &kTestRE[dirSE  *kQ];
         kDistTest.f[dirNW  ] = &kTestRE[dirNW  *kQ];
         kDistTest.f[dirTE  ] = &kTestRE[dirTE  *kQ];
         kDistTest.f[dirBW  ] = &kTestRE[dirBW  *kQ];
         kDistTest.f[dirBE  ] = &kTestRE[dirBE  *kQ];
         kDistTest.f[dirTW  ] = &kTestRE[dirTW  *kQ];
         kDistTest.f[dirTN  ] = &kTestRE[dirTN  *kQ];
         kDistTest.f[dirBS  ] = &kTestRE[dirBS  *kQ];
         kDistTest.f[dirBN  ] = &kTestRE[dirBN  *kQ];
         kDistTest.f[dirTS  ] = &kTestRE[dirTS  *kQ];
         kDistTest.f[dirZERO] = &kTestRE[dirZERO*kQ];
         kDistTest.f[dirTNE ] = &kTestRE[dirTNE *kQ];
         kDistTest.f[dirTSW ] = &kTestRE[dirTSW *kQ];
         kDistTest.f[dirTSE ] = &kTestRE[dirTSE *kQ];
         kDistTest.f[dirTNW ] = &kTestRE[dirTNW *kQ];
         kDistTest.f[dirBNE ] = &kTestRE[dirBNE *kQ];
         kDistTest.f[dirBSW ] = &kTestRE[dirBSW *kQ];
         kDistTest.f[dirBSE ] = &kTestRE[dirBSE *kQ];
         kDistTest.f[dirBNW ] = &kTestRE[dirBNW *kQ];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   //doubflo f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   //f1_W    = (D.f[dirE   ])[k1e   ];
   //   //f1_E    = (D.f[dirW   ])[k1w   ];
   //   //f1_S    = (D.f[dirN   ])[k1n   ];
   //   //f1_N    = (D.f[dirS   ])[k1s   ];
   //   //f1_B    = (D.f[dirT   ])[k1t   ];
   //   //f1_T    = (D.f[dirB   ])[k1b   ];
   //   //f1_SW   = (D.f[dirNE  ])[k1ne  ];
   //   //f1_NE   = (D.f[dirSW  ])[k1sw  ];
   //   //f1_NW   = (D.f[dirSE  ])[k1se  ];
   //   //f1_SE   = (D.f[dirNW  ])[k1nw  ];
   //   //f1_BW   = (D.f[dirTE  ])[k1te  ];
   //   //f1_TE   = (D.f[dirBW  ])[k1bw  ];
   //   //f1_TW   = (D.f[dirBE  ])[k1be  ];
   //   //f1_BE   = (D.f[dirTW  ])[k1tw  ];
   //   //f1_BS   = (D.f[dirTN  ])[k1tn  ];
   //   //f1_TN   = (D.f[dirBS  ])[k1bs  ];
   //   //f1_TS   = (D.f[dirBN  ])[k1bn  ];
   //   //f1_BN   = (D.f[dirTS  ])[k1ts  ];
   //   //f1_ZERO = (D.f[dirZERO])[k1zero];
   //   //f1_BSW  = (D.f[dirTNE ])[k1tne ];
   //   //f1_BNE  = (D.f[dirTSW ])[k1tsw ];
   //   //f1_BNW  = (D.f[dirTSE ])[k1tse ];
   //   //f1_BSE  = (D.f[dirTNW ])[k1tnw ];
   //   //f1_TSW  = (D.f[dirBNE ])[k1bne ];
   //   //f1_TNE  = (D.f[dirBSW ])[k1bsw ];
   //   //f1_TNW  = (D.f[dirBSE ])[k1bse ];
   //   //f1_TSE  = (D.f[dirBNW ])[k1bnw ];
   //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //   doubflo f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;
   //   f1_E    = (D.f[dirE   ])[k1e   ];
   //   f1_W    = (D.f[dirW   ])[k1w   ];
   //   f1_N    = (D.f[dirN   ])[k1n   ];
   //   f1_S    = (D.f[dirS   ])[k1s   ];
   //   f1_T    = (D.f[dirT   ])[k1t   ];
   //   f1_B    = (D.f[dirB   ])[k1b   ];
   //   f1_NE   = (D.f[dirNE  ])[k1ne  ];
   //   f1_SW   = (D.f[dirSW  ])[k1sw  ];
   //   f1_SE   = (D.f[dirSE  ])[k1se  ];
   //   f1_NW   = (D.f[dirNW  ])[k1nw  ];
   //   f1_TE   = (D.f[dirTE  ])[k1te  ];
   //   f1_BW   = (D.f[dirBW  ])[k1bw  ];
   //   f1_BE   = (D.f[dirBE  ])[k1be  ];
   //   f1_TW   = (D.f[dirTW  ])[k1tw  ];
   //   f1_TN   = (D.f[dirTN  ])[k1tn  ];
   //   f1_BS   = (D.f[dirBS  ])[k1bs  ];
   //   f1_BN   = (D.f[dirBN  ])[k1bn  ];
   //   f1_TS   = (D.f[dirTS  ])[k1ts  ];
   //   f1_ZERO = (D.f[dirZERO])[k1zero];
   //   f1_TNE  = (D.f[dirTNE ])[k1tne ];
   //   f1_TSW  = (D.f[dirTSW ])[k1tsw ];
   //   f1_TSE  = (D.f[dirTSE ])[k1tse ];
   //   f1_TNW  = (D.f[dirTNW ])[k1tnw ];
   //   f1_BNE  = (D.f[dirBNE ])[k1bne ];
   //   f1_BSW  = (D.f[dirBSW ])[k1bsw ];
   //   f1_BSE  = (D.f[dirBSE ])[k1bse ];
   //   f1_BNW  = (D.f[dirBNW ])[k1bnw ];
   //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //   //////////////////////////////////////////////////////////////////////////
   //   doubflo drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+ f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;
	  //doubflo vx1      = (((f1_TNE-f1_BSW)+(f1_BSE-f1_TNW)+(f1_BNE-f1_TSW)+(f1_TSE-f1_BNW)) + (((f1_NE-f1_SW)+(f1_TE-f1_BW))+((f1_SE-f1_NW)+(f1_BE-f1_TW))) + (f1_E-f1_W)) / (one + drho1);
	  //doubflo vx2      = (((f1_TNE-f1_BSW)+(f1_TNW-f1_BSE)+(f1_BNE-f1_TSW)+(f1_BNW-f1_TSE)) + (((f1_NE-f1_SW)+(f1_TN-f1_BS))+((f1_BN-f1_TS)+(f1_NW-f1_SE))) + (f1_N-f1_S)) / (one + drho1);
	  //doubflo vx3      = (((f1_TNE-f1_BSW)+(f1_TNW-f1_BSE)+(f1_TSW-f1_BNE)+(f1_TSE-f1_BNW)) + (((f1_TE-f1_BW)+(f1_TN-f1_BS))+((f1_TW-f1_BE)+(f1_TS-f1_BN))) + (f1_T-f1_B)) / (one + drho1);
   //   //////////////////////////////////////////////////////////////////////////
	  ////doubflo omega = om1;
   //   doubflo cusq  = c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
   //   //////////////////////////////////////////////////////////////////////////
	  ////Täst MK
	  ////if(vx1 < zero) vx1 = zero;
   //   //////////////////////////////////////////////////////////////////////////
	  ////becomes higher with neighbor source and lower with local source
   //   //doubflo fZERO = c8over27*  (rhoBC[k]-(one + rhoBC[k])*(cusq))                                                           ;
   //   //doubflo fE    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1        )+c9over2*( vx1        )*( vx1        )-cusq));
   //   //doubflo fW    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cusq));
   //   //doubflo fN    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(    vx2     )+c9over2*(     vx2    )*(     vx2    )-cusq));
   //   //doubflo fS    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(   -vx2     )+c9over2*(    -vx2    )*(    -vx2    )-cusq));
   //   //doubflo fT    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(         vx3)+c9over2*(         vx3)*(         vx3)-cusq));
   //   //doubflo fB    = c2over27*  (rhoBC[k]+(one + rhoBC[k])*(three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cusq));
   //   //doubflo fNE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2    )+c9over2*( vx1+vx2    )*( vx1+vx2    )-cusq));
   //   //doubflo fSW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2    )+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
   //   //doubflo fSE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2    )+c9over2*( vx1-vx2    )*( vx1-vx2    )-cusq));
   //   //doubflo fNW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2    )+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
   //   //doubflo fTE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1    +vx3)+c9over2*( vx1    +vx3)*( vx1    +vx3)-cusq));
   //   //doubflo fBW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1    -vx3)+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
   //   //doubflo fBE   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*( vx1    -vx3)+c9over2*( vx1    -vx3)*( vx1    -vx3)-cusq));
   //   //doubflo fTW   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1    +vx3)+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
   //   //doubflo fTN   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(     vx2+vx3)+c9over2*(     vx2+vx3)*(     vx2+vx3)-cusq));
   //   //doubflo fBS   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(    -vx2-vx3)+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
   //   //doubflo fBN   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(     vx2-vx3)+c9over2*(     vx2-vx3)*(     vx2-vx3)-cusq));
   //   //doubflo fTS   = c1over54*  (rhoBC[k]+(one + rhoBC[k])*(three*(    -vx2+vx3)+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
   //   //doubflo fTNE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2+vx3)+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
   //   //doubflo fBSW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2-vx3)+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
   //   //doubflo fBNE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1+vx2-vx3)+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
   //   //doubflo fTSW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1-vx2+vx3)+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
   //   //doubflo fTSE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2+vx3)+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
   //   //doubflo fBNW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2-vx3)+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
   //   //doubflo fBSE  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*( vx1-vx2-vx3)+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
   //   //doubflo fTNW  = c1over216* (rhoBC[k]+(one + rhoBC[k])*(three*(-vx1+vx2+vx3)+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));
   //   //////////////////////////////////////////////////////////////////////////
	  //// based on VirtualFluids (kucher + fard)
   //   doubflo fZERO = c8over27  * rhoBC[k] * (one                                                                      - cusq);
   //   doubflo fE    = c2over27  * rhoBC[k] * (one + three * ( vx1        ) + c9over2 * ( vx1        ) * ( vx1        ) - cusq);
   //   doubflo fW    = c2over27  * rhoBC[k] * (one + three * (-vx1        ) + c9over2 * (-vx1        ) * (-vx1        ) - cusq);
   //   doubflo fN    = c2over27  * rhoBC[k] * (one + three * (     vx2    ) + c9over2 * (     vx2    ) * (     vx2    ) - cusq);
   //   doubflo fS    = c2over27  * rhoBC[k] * (one + three * (    -vx2    ) + c9over2 * (    -vx2    ) * (    -vx2    ) - cusq);
   //   doubflo fT    = c2over27  * rhoBC[k] * (one + three * (         vx3) + c9over2 * (         vx3) * (         vx3) - cusq);
   //   doubflo fB    = c2over27  * rhoBC[k] * (one + three * (        -vx3) + c9over2 * (        -vx3) * (        -vx3) - cusq);
   //   doubflo fNE   = c1over54  * rhoBC[k] * (one + three * ( vx1+vx2    ) + c9over2 * ( vx1+vx2    ) * ( vx1+vx2    ) - cusq);
   //   doubflo fSW   = c1over54  * rhoBC[k] * (one + three * (-vx1-vx2    ) + c9over2 * (-vx1-vx2    ) * (-vx1-vx2    ) - cusq);
   //   doubflo fSE   = c1over54  * rhoBC[k] * (one + three * ( vx1-vx2    ) + c9over2 * ( vx1-vx2    ) * ( vx1-vx2    ) - cusq);
   //   doubflo fNW   = c1over54  * rhoBC[k] * (one + three * (-vx1+vx2    ) + c9over2 * (-vx1+vx2    ) * (-vx1+vx2    ) - cusq);
   //   doubflo fTE   = c1over54  * rhoBC[k] * (one + three * ( vx1    +vx3) + c9over2 * ( vx1    +vx3) * ( vx1    +vx3) - cusq);
   //   doubflo fBW   = c1over54  * rhoBC[k] * (one + three * (-vx1    -vx3) + c9over2 * (-vx1    -vx3) * (-vx1    -vx3) - cusq);
   //   doubflo fBE   = c1over54  * rhoBC[k] * (one + three * ( vx1    -vx3) + c9over2 * ( vx1    -vx3) * ( vx1    -vx3) - cusq);
   //   doubflo fTW   = c1over54  * rhoBC[k] * (one + three * (-vx1    +vx3) + c9over2 * (-vx1    +vx3) * (-vx1    +vx3) - cusq);
   //   doubflo fTN   = c1over54  * rhoBC[k] * (one + three * (     vx2+vx3) + c9over2 * (     vx2+vx3) * (     vx2+vx3) - cusq);
   //   doubflo fBS   = c1over54  * rhoBC[k] * (one + three * (    -vx2-vx3) + c9over2 * (    -vx2-vx3) * (    -vx2-vx3) - cusq);
   //   doubflo fBN   = c1over54  * rhoBC[k] * (one + three * (     vx2-vx3) + c9over2 * (     vx2-vx3) * (     vx2-vx3) - cusq);
   //   doubflo fTS   = c1over54  * rhoBC[k] * (one + three * (    -vx2+vx3) + c9over2 * (    -vx2+vx3) * (    -vx2+vx3) - cusq);
   //   doubflo fTNE  = c1over216 * rhoBC[k] * (one + three * ( vx1+vx2+vx3) + c9over2 * ( vx1+vx2+vx3) * ( vx1+vx2+vx3) - cusq);
   //   doubflo fBSW  = c1over216 * rhoBC[k] * (one + three * (-vx1-vx2-vx3) + c9over2 * (-vx1-vx2-vx3) * (-vx1-vx2-vx3) - cusq);
   //   doubflo fBNE  = c1over216 * rhoBC[k] * (one + three * ( vx1+vx2-vx3) + c9over2 * ( vx1+vx2-vx3) * ( vx1+vx2-vx3) - cusq);
   //   doubflo fTSW  = c1over216 * rhoBC[k] * (one + three * (-vx1-vx2+vx3) + c9over2 * (-vx1-vx2+vx3) * (-vx1-vx2+vx3) - cusq);
   //   doubflo fTSE  = c1over216 * rhoBC[k] * (one + three * ( vx1-vx2+vx3) + c9over2 * ( vx1-vx2+vx3) * ( vx1-vx2+vx3) - cusq);
   //   doubflo fBNW  = c1over216 * rhoBC[k] * (one + three * (-vx1+vx2-vx3) + c9over2 * (-vx1+vx2-vx3) * (-vx1+vx2-vx3) - cusq);
   //   doubflo fBSE  = c1over216 * rhoBC[k] * (one + three * ( vx1-vx2-vx3) + c9over2 * ( vx1-vx2-vx3) * ( vx1-vx2-vx3) - cusq);
   //   doubflo fTNW  = c1over216 * rhoBC[k] * (one + three * (-vx1+vx2+vx3) + c9over2 * (-vx1+vx2+vx3) * (-vx1+vx2+vx3) - cusq);
   ////   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////test
   ////   doubflo fZERO = c8over27  * ((drho1 + rhoBC[k]) / two) * (one                                                                      - cusq);
   ////   doubflo fE    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1        ) + c9over2 * ( vx1        ) * ( vx1        ) - cusq);
   ////   doubflo fW    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1        ) + c9over2 * (-vx1        ) * (-vx1        ) - cusq);
   ////   doubflo fN    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2    ) + c9over2 * (     vx2    ) * (     vx2    ) - cusq);
   ////   doubflo fS    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2    ) + c9over2 * (    -vx2    ) * (    -vx2    ) - cusq);
   ////   doubflo fT    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (         vx3) + c9over2 * (         vx3) * (         vx3) - cusq);
   ////   doubflo fB    = c2over27  * ((drho1 + rhoBC[k]) / two) * (one + three * (        -vx3) + c9over2 * (        -vx3) * (        -vx3) - cusq);
   ////   doubflo fNE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2    ) + c9over2 * ( vx1+vx2    ) * ( vx1+vx2    ) - cusq);
   ////   doubflo fSW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2    ) + c9over2 * (-vx1-vx2    ) * (-vx1-vx2    ) - cusq);
   ////   doubflo fSE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2    ) + c9over2 * ( vx1-vx2    ) * ( vx1-vx2    ) - cusq);
   ////   doubflo fNW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2    ) + c9over2 * (-vx1+vx2    ) * (-vx1+vx2    ) - cusq);
   ////   doubflo fTE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1    +vx3) + c9over2 * ( vx1    +vx3) * ( vx1    +vx3) - cusq);
   ////   doubflo fBW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1    -vx3) + c9over2 * (-vx1    -vx3) * (-vx1    -vx3) - cusq);
   ////   doubflo fBE   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1    -vx3) + c9over2 * ( vx1    -vx3) * ( vx1    -vx3) - cusq);
   ////   doubflo fTW   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1    +vx3) + c9over2 * (-vx1    +vx3) * (-vx1    +vx3) - cusq);
   ////   doubflo fTN   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2+vx3) + c9over2 * (     vx2+vx3) * (     vx2+vx3) - cusq);
   ////   doubflo fBS   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2-vx3) + c9over2 * (    -vx2-vx3) * (    -vx2-vx3) - cusq);
   ////   doubflo fBN   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (     vx2-vx3) + c9over2 * (     vx2-vx3) * (     vx2-vx3) - cusq);
   ////   doubflo fTS   = c1over54  * ((drho1 + rhoBC[k]) / two) * (one + three * (    -vx2+vx3) + c9over2 * (    -vx2+vx3) * (    -vx2+vx3) - cusq);
   ////   doubflo fTNE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2+vx3) + c9over2 * ( vx1+vx2+vx3) * ( vx1+vx2+vx3) - cusq);
   ////   doubflo fBSW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2-vx3) + c9over2 * (-vx1-vx2-vx3) * (-vx1-vx2-vx3) - cusq);
   ////   doubflo fBNE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1+vx2-vx3) + c9over2 * ( vx1+vx2-vx3) * ( vx1+vx2-vx3) - cusq);
   ////   doubflo fTSW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1-vx2+vx3) + c9over2 * (-vx1-vx2+vx3) * (-vx1-vx2+vx3) - cusq);
   ////   doubflo fTSE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2+vx3) + c9over2 * ( vx1-vx2+vx3) * ( vx1-vx2+vx3) - cusq);
   ////   doubflo fBNW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2-vx3) + c9over2 * (-vx1+vx2-vx3) * (-vx1+vx2-vx3) - cusq);
   ////   doubflo fBSE  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * ( vx1-vx2-vx3) + c9over2 * ( vx1-vx2-vx3) * ( vx1-vx2-vx3) - cusq);
   ////   doubflo fTNW  = c1over216 * ((drho1 + rhoBC[k]) / two) * (one + three * (-vx1+vx2+vx3) + c9over2 * (-vx1+vx2+vx3) * (-vx1+vx2+vx3) - cusq);

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // based on BGK Plus Comp
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//double mfabb = (D.f[dirE   ])[k1e   ];
			//double mfcbb = (D.f[dirW   ])[k1w   ];
			//double mfbab = (D.f[dirN   ])[k1n   ];
			//double mfbcb = (D.f[dirS   ])[k1s   ];
			//double mfbba = (D.f[dirT   ])[k1t   ];
			//double mfbbc = (D.f[dirB   ])[k1b   ];
			//double mfaab = (D.f[dirNE  ])[k1ne  ];
			//double mfccb = (D.f[dirSW  ])[k1sw  ];
			//double mfacb = (D.f[dirSE  ])[k1se  ];
			//double mfcab = (D.f[dirNW  ])[k1nw  ];
			//double mfaba = (D.f[dirTE  ])[k1te  ];
			//double mfcbc = (D.f[dirBW  ])[k1bw  ];
			//double mfabc = (D.f[dirBE  ])[k1be  ];
			//double mfcba = (D.f[dirTW  ])[k1tw  ];
			//double mfbaa = (D.f[dirTN  ])[k1tn  ];
			//double mfbcc = (D.f[dirBS  ])[k1bs  ];
			//double mfbac = (D.f[dirBN  ])[k1bn  ];
			//double mfbca = (D.f[dirTS  ])[k1ts  ];
			//double mfbbb = (D.f[dirZERO])[k1zero];
			//double mfaaa = (D.f[dirTNE ])[k1tne ];
			//double mfcca = (D.f[dirTSW ])[k1tsw ];
			//double mfaca = (D.f[dirTSE ])[k1tse ];
			//double mfcaa = (D.f[dirTNW ])[k1tnw ];
			//double mfaac = (D.f[dirBNE ])[k1bne ];
			//double mfccc = (D.f[dirBSW ])[k1bsw ];
			//double mfacc = (D.f[dirBSE ])[k1bse ];
			//double mfcac = (D.f[dirBNW ])[k1bnw ];
			doubflo mfabb = (D.f[dirE   ])[k1e   ];
			doubflo mfcbb = (D.f[dirW   ])[k1w   ];
			doubflo mfbab = (D.f[dirN   ])[k1n   ];
			doubflo mfbcb = (D.f[dirS   ])[k1s   ];
			doubflo mfbba = (D.f[dirT   ])[k1t   ];
			doubflo mfbbc = (D.f[dirB   ])[k1b   ];
			doubflo mfaab = (D.f[dirNE  ])[k1ne  ];
			doubflo mfccb = (D.f[dirSW  ])[k1sw  ];
			doubflo mfacb = (D.f[dirSE  ])[k1se  ];
			doubflo mfcab = (D.f[dirNW  ])[k1nw  ];
			doubflo mfaba = (D.f[dirTE  ])[k1te  ];
			doubflo mfcbc = (D.f[dirBW  ])[k1bw  ];
			doubflo mfabc = (D.f[dirBE  ])[k1be  ];
			doubflo mfcba = (D.f[dirTW  ])[k1tw  ];
			doubflo mfbaa = (D.f[dirTN  ])[k1tn  ];
			doubflo mfbcc = (D.f[dirBS  ])[k1bs  ];
			doubflo mfbac = (D.f[dirBN  ])[k1bn  ];
			doubflo mfbca = (D.f[dirTS  ])[k1ts  ];
			doubflo mfbbb = (D.f[dirZERO])[k1zero];
			doubflo mfaaa = (D.f[dirTNE ])[k1tne ];
			doubflo mfcca = (D.f[dirTSW ])[k1tsw ];
			doubflo mfaca = (D.f[dirTSE ])[k1tse ];
			doubflo mfcaa = (D.f[dirTNW ])[k1tnw ];
			doubflo mfaac = (D.f[dirBNE ])[k1bne ];
			doubflo mfccc = (D.f[dirBSW ])[k1bsw ];
			doubflo mfacc = (D.f[dirBSE ])[k1bse ];
			doubflo mfcac = (D.f[dirBNW ])[k1bnw ];

			//doubflo mfcbb = (D.f[dirE   ])[ke   ];
			//doubflo mfabb = (D.f[dirW   ])[kw   ];
			//doubflo mfbcb = (D.f[dirN   ])[kn   ];
			//doubflo mfbab = (D.f[dirS   ])[ks   ];
			//doubflo mfbbc = (D.f[dirT   ])[kt   ];
			//doubflo mfbba = (D.f[dirB   ])[kb   ];
			//doubflo mfccb = (D.f[dirNE  ])[kne  ];
			//doubflo mfaab = (D.f[dirSW  ])[ksw  ];
			//doubflo mfcab = (D.f[dirSE  ])[kse  ];
			//doubflo mfacb = (D.f[dirNW  ])[knw  ];
			//doubflo mfcbc = (D.f[dirTE  ])[kte  ];
			//doubflo mfaba = (D.f[dirBW  ])[kbw  ];
			//doubflo mfcba = (D.f[dirBE  ])[kbe  ];
			//doubflo mfabc = (D.f[dirTW  ])[ktw  ];
			//doubflo mfbcc = (D.f[dirTN  ])[ktn  ];
			//doubflo mfbaa = (D.f[dirBS  ])[kbs  ];
			//doubflo mfbca = (D.f[dirBN  ])[kbn  ];
			//doubflo mfbac = (D.f[dirTS  ])[kts  ];
			//doubflo mfbbb = (D.f[dirZERO])[kzero];
			//doubflo mfccc = (D.f[dirTNE ])[ktne ];
			//doubflo mfaac = (D.f[dirTSW ])[ktsw ];
			//doubflo mfcac = (D.f[dirTSE ])[ktse ];
			//doubflo mfacc = (D.f[dirTNW ])[ktnw ];
			//doubflo mfcca = (D.f[dirBNE ])[kbne ];
			//doubflo mfaaa = (D.f[dirBSW ])[kbsw ];
			//doubflo mfcaa = (D.f[dirBSE ])[kbse ];
			//doubflo mfaca = (D.f[dirBNW ])[kbnw ];
			////////////////////////////////////////////////////////////////////////////////////
			//doubflo rho   = (((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
			//				(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
			//				((mfabb+mfcbb) + (mfbab+mfbcb)) + (mfbba+mfbbc)) + mfbbb) + one;//!!!!Achtung + one
			////////////////////////////////////////////////////////////////////////////////////
			doubflo rho = rhoBC[k];
			////////////////////////////////////////////////////////////////////////////////////
			doubflo OoRho = one / (rho * 1.5f);
			////////////////////////////////////////////////////////////////////////////////////
			doubflo vvx    = ((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) * OoRho;
			doubflo vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) * OoRho;
			doubflo vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
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
			//	(kDistTest.f[dirE   ])[k] = mfabb;
			//	(kDistTest.f[dirW   ])[k] = mfcbb;
			//	(kDistTest.f[dirN   ])[k] = mfbab;
			//	(kDistTest.f[dirS   ])[k] = mfbcb;
			//	(kDistTest.f[dirT   ])[k] = mfbba;
			//	(kDistTest.f[dirB   ])[k] = mfbbc;
			//	(kDistTest.f[dirNE  ])[k] = mfaab;
			//	(kDistTest.f[dirSW  ])[k] = mfccb;
			//	(kDistTest.f[dirSE  ])[k] = mfacb;
			//	(kDistTest.f[dirNW  ])[k] = mfcab;
			//	(kDistTest.f[dirTE  ])[k] = mfaba;
			//	(kDistTest.f[dirBW  ])[k] = mfcbc;
			//	(kDistTest.f[dirBE  ])[k] = mfabc;
			//	(kDistTest.f[dirTW  ])[k] = mfcba;
			//	(kDistTest.f[dirTN  ])[k] = mfbaa;
			//	(kDistTest.f[dirBS  ])[k] = mfbcc;
			//	(kDistTest.f[dirBN  ])[k] = mfbac;
			//	(kDistTest.f[dirTS  ])[k] = mfbca;
			//	(kDistTest.f[dirZERO])[k] = KQK;
			//	(kDistTest.f[dirTNE ])[k] = mfaaa;
			//	(kDistTest.f[dirTSW ])[k] = mfcca;
			//	(kDistTest.f[dirTSE ])[k] = mfaca;
			//	(kDistTest.f[dirTNW ])[k] = mfcaa;
			//	(kDistTest.f[dirBNE ])[k] = mfaac;
			//	(kDistTest.f[dirBSW ])[k] = mfccc;
			//	(kDistTest.f[dirBSE ])[k] = mfacc;
			//	(kDistTest.f[dirBNW ])[k] = mfcac;
			//}else{
			//	(kDistTest.f[dirE   ])[k] = zero;
			//	(kDistTest.f[dirW   ])[k] = zero;
			//	(kDistTest.f[dirN   ])[k] = zero;
			//	(kDistTest.f[dirS   ])[k] = zero;
			//	(kDistTest.f[dirT   ])[k] = zero;
			//	(kDistTest.f[dirB   ])[k] = zero;
			//	(kDistTest.f[dirNE  ])[k] = zero;
			//	(kDistTest.f[dirSW  ])[k] = zero;
			//	(kDistTest.f[dirSE  ])[k] = zero;
			//	(kDistTest.f[dirNW  ])[k] = zero;
			//	(kDistTest.f[dirTE  ])[k] = zero;
			//	(kDistTest.f[dirBW  ])[k] = zero;
			//	(kDistTest.f[dirBE  ])[k] = zero;
			//	(kDistTest.f[dirTW  ])[k] = zero;
			//	(kDistTest.f[dirTN  ])[k] = zero;
			//	(kDistTest.f[dirBS  ])[k] = zero;
			//	(kDistTest.f[dirBN  ])[k] = zero;
			//	(kDistTest.f[dirTS  ])[k] = zero;
			//	(kDistTest.f[dirZERO])[k] = zero;
			//	(kDistTest.f[dirTNE ])[k] = zero;
			//	(kDistTest.f[dirTSW ])[k] = zero;
			//	(kDistTest.f[dirTSE ])[k] = zero;
			//	(kDistTest.f[dirTNW ])[k] = zero;
			//	(kDistTest.f[dirBNE ])[k] = zero;
			//	(kDistTest.f[dirBSW ])[k] = zero;
			//	(kDistTest.f[dirBSE ])[k] = zero;
			//	(kDistTest.f[dirBNW ])[k] = zero;
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
            doubflo XXb    = -c2o3 + vx2;
            doubflo XXc    = -c1o2 * (XXb + one + vvx);
            doubflo XXa    = XXc + vvx;
            doubflo YYb    = -c2o3 + vy2;
            doubflo YYc    = -c1o2 * (YYb + one + vvy);
            doubflo YYa    = YYc + vvy;
            doubflo ZZb    = -c2o3 + vz2;
            doubflo ZZc    = -c1o2 * (ZZb + one + vvz);
            doubflo ZZa    = ZZc + vvz;
			//////////////////////////////////////////////////////////////////////////////////
			//unkonditioniert
            mfcbb = -(rhoBC[k] + one) * XXc * YYb * ZZb - c2over27; 
			mfabb = -(rhoBC[k] + one) * XXa * YYb * ZZb - c2over27;
			mfbcb = -(rhoBC[k] + one) * XXb * YYc * ZZb - c2over27;
			mfbab = -(rhoBC[k] + one) * XXb * YYa * ZZb - c2over27;
			mfbbc = -(rhoBC[k] + one) * XXb * YYb * ZZc - c2over27;
			mfbba = -(rhoBC[k] + one) * XXb * YYb * ZZa - c2over27;
			mfccb = -(rhoBC[k] + one) * XXc * YYc * ZZb - c1over54;
			mfaab = -(rhoBC[k] + one) * XXa * YYa * ZZb - c1over54;
			mfcab = -(rhoBC[k] + one) * XXc * YYa * ZZb - c1over54;
			mfacb = -(rhoBC[k] + one) * XXa * YYc * ZZb - c1over54;
			mfcbc = -(rhoBC[k] + one) * XXc * YYb * ZZc - c1over54;
			mfaba = -(rhoBC[k] + one) * XXa * YYb * ZZa - c1over54;
			mfcba = -(rhoBC[k] + one) * XXc * YYb * ZZa - c1over54;
			mfabc = -(rhoBC[k] + one) * XXa * YYb * ZZc - c1over54;
			mfbcc = -(rhoBC[k] + one) * XXb * YYc * ZZc - c1over54;
			mfbaa = -(rhoBC[k] + one) * XXb * YYa * ZZa - c1over54;
			mfbca = -(rhoBC[k] + one) * XXb * YYc * ZZa - c1over54;
			mfbac = -(rhoBC[k] + one) * XXb * YYa * ZZc - c1over54;
			mfbbb = -(rhoBC[k] + one) * XXb * YYb * ZZb - c8over27;
			mfccc = -(rhoBC[k] + one) * XXc * YYc * ZZc - c1over216;
			mfaac = -(rhoBC[k] + one) * XXa * YYa * ZZc - c1over216;
			mfcac = -(rhoBC[k] + one) * XXc * YYa * ZZc - c1over216;
			mfacc = -(rhoBC[k] + one) * XXa * YYc * ZZc - c1over216;
			mfcca = -(rhoBC[k] + one) * XXc * YYc * ZZa - c1over216;
			mfaaa = -(rhoBC[k] + one) * XXa * YYa * ZZa - c1over216;
			mfcaa = -(rhoBC[k] + one) * XXc * YYa * ZZa - c1over216;
			mfaca = -(rhoBC[k] + one) * XXa * YYc * ZZa - c1over216;
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
      //if (evenOrOdd==true)
      //{
      //   D.f[dirE   ] = &DD[dirE   *size_Mat];
      //   D.f[dirW   ] = &DD[dirW   *size_Mat];
      //   D.f[dirN   ] = &DD[dirN   *size_Mat];
      //   D.f[dirS   ] = &DD[dirS   *size_Mat];
      //   D.f[dirT   ] = &DD[dirT   *size_Mat];
      //   D.f[dirB   ] = &DD[dirB   *size_Mat];
      //   D.f[dirNE  ] = &DD[dirNE  *size_Mat];
      //   D.f[dirSW  ] = &DD[dirSW  *size_Mat];
      //   D.f[dirSE  ] = &DD[dirSE  *size_Mat];
      //   D.f[dirNW  ] = &DD[dirNW  *size_Mat];
      //   D.f[dirTE  ] = &DD[dirTE  *size_Mat];
      //   D.f[dirBW  ] = &DD[dirBW  *size_Mat];
      //   D.f[dirBE  ] = &DD[dirBE  *size_Mat];
      //   D.f[dirTW  ] = &DD[dirTW  *size_Mat];
      //   D.f[dirTN  ] = &DD[dirTN  *size_Mat];
      //   D.f[dirBS  ] = &DD[dirBS  *size_Mat];
      //   D.f[dirBN  ] = &DD[dirBN  *size_Mat];
      //   D.f[dirTS  ] = &DD[dirTS  *size_Mat];
      //   D.f[dirZERO] = &DD[dirZERO*size_Mat];
      //   D.f[dirTNE ] = &DD[dirTNE *size_Mat];
      //   D.f[dirTSW ] = &DD[dirTSW *size_Mat];
      //   D.f[dirTSE ] = &DD[dirTSE *size_Mat];
      //   D.f[dirTNW ] = &DD[dirTNW *size_Mat];
      //   D.f[dirBNE ] = &DD[dirBNE *size_Mat];
      //   D.f[dirBSW ] = &DD[dirBSW *size_Mat];
      //   D.f[dirBSE ] = &DD[dirBSE *size_Mat];
      //   D.f[dirBNW ] = &DD[dirBNW *size_Mat];
      //} 
      //else
      //{
      //   D.f[dirW   ] = &DD[dirE   *size_Mat];
      //   D.f[dirE   ] = &DD[dirW   *size_Mat];
      //   D.f[dirS   ] = &DD[dirN   *size_Mat];
      //   D.f[dirN   ] = &DD[dirS   *size_Mat];
      //   D.f[dirB   ] = &DD[dirT   *size_Mat];
      //   D.f[dirT   ] = &DD[dirB   *size_Mat];
      //   D.f[dirSW  ] = &DD[dirNE  *size_Mat];
      //   D.f[dirNE  ] = &DD[dirSW  *size_Mat];
      //   D.f[dirNW  ] = &DD[dirSE  *size_Mat];
      //   D.f[dirSE  ] = &DD[dirNW  *size_Mat];
      //   D.f[dirBW  ] = &DD[dirTE  *size_Mat];
      //   D.f[dirTE  ] = &DD[dirBW  *size_Mat];
      //   D.f[dirTW  ] = &DD[dirBE  *size_Mat];
      //   D.f[dirBE  ] = &DD[dirTW  *size_Mat];
      //   D.f[dirBS  ] = &DD[dirTN  *size_Mat];
      //   D.f[dirTN  ] = &DD[dirBS  *size_Mat];
      //   D.f[dirTS  ] = &DD[dirBN  *size_Mat];
      //   D.f[dirBN  ] = &DD[dirTS  *size_Mat];
      //   D.f[dirZERO] = &DD[dirZERO*size_Mat];
      //   D.f[dirTNE ] = &DD[dirBSW *size_Mat];
      //   D.f[dirTSW ] = &DD[dirBNE *size_Mat];
      //   D.f[dirTSE ] = &DD[dirBNW *size_Mat];
      //   D.f[dirTNW ] = &DD[dirBSE *size_Mat];
      //   D.f[dirBNE ] = &DD[dirTSW *size_Mat];
      //   D.f[dirBSW ] = &DD[dirTNE *size_Mat];
      //   D.f[dirBSE ] = &DD[dirTNW *size_Mat];
      //   D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      //}
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //__syncthreads();

			(D.f[dirE   ])[ke   ] = mfabb;//mfcbb;
			(D.f[dirW   ])[kw   ] = mfcbb;//mfabb;
			(D.f[dirN   ])[kn   ] = mfbab;//mfbcb;
			(D.f[dirS   ])[ks   ] = mfbcb;//mfbab;
			(D.f[dirT   ])[kt   ] = mfbba;//mfbbc;
			(D.f[dirB   ])[kb   ] = mfbbc;//mfbba;
			(D.f[dirNE  ])[kne  ] = mfaab;//mfccb;
			(D.f[dirSW  ])[ksw  ] = mfccb;//mfaab;
			(D.f[dirSE  ])[kse  ] = mfacb;//mfcab;
			(D.f[dirNW  ])[knw  ] = mfcab;//mfacb;
			(D.f[dirTE  ])[kte  ] = mfaba;//mfcbc;
			(D.f[dirBW  ])[kbw  ] = mfcbc;//mfaba;
			(D.f[dirBE  ])[kbe  ] = mfabc;//mfcba;
			(D.f[dirTW  ])[ktw  ] = mfcba;//mfabc;
			(D.f[dirTN  ])[ktn  ] = mfbaa;//mfbcc;
			(D.f[dirBS  ])[kbs  ] = mfbcc;//mfbaa;
			(D.f[dirBN  ])[kbn  ] = mfbac;//mfbca;
			(D.f[dirTS  ])[kts  ] = mfbca;//mfbac;
			(D.f[dirZERO])[kzero] = mfbbb;//mfbbb;
			(D.f[dirTNE ])[ktne ] = mfaaa;//mfccc;
			(D.f[dirTSW ])[ktsw ] = mfcca;//mfaac;
			(D.f[dirTSE ])[ktse ] = mfaca;//mfcac;
			(D.f[dirTNW ])[ktnw ] = mfcaa;//mfacc;
			(D.f[dirBNE ])[kbne ] = mfaac;//mfcca;
			(D.f[dirBSW ])[kbsw ] = mfccc;//mfaaa;
			(D.f[dirBSE ])[kbse ] = mfacc;//mfcaa;
			(D.f[dirBNW ])[kbnw ] = mfcac;//mfaca;
			//(D.f[dirE   ])[ke   ] = mfcbb;
			//(D.f[dirW   ])[kw   ] = mfabb;
			//(D.f[dirN   ])[kn   ] = mfbcb;
			//(D.f[dirS   ])[ks   ] = mfbab;
			//(D.f[dirT   ])[kt   ] = mfbbc;
			//(D.f[dirB   ])[kb   ] = mfbba;
			//(D.f[dirNE  ])[kne  ] = mfccb;
			//(D.f[dirSW  ])[ksw  ] = mfaab;
			//(D.f[dirSE  ])[kse  ] = mfcab;
			//(D.f[dirNW  ])[knw  ] = mfacb;
			//(D.f[dirTE  ])[kte  ] = mfcbc;
			//(D.f[dirBW  ])[kbw  ] = mfaba;
			//(D.f[dirBE  ])[kbe  ] = mfcba;
			//(D.f[dirTW  ])[ktw  ] = mfabc;
			//(D.f[dirTN  ])[ktn  ] = mfbcc;
			//(D.f[dirBS  ])[kbs  ] = mfbaa;
			//(D.f[dirBN  ])[kbn  ] = mfbca;
			//(D.f[dirTS  ])[kts  ] = mfbac;
			//(D.f[dirZERO])[kzero] = mfbbb;
			//(D.f[dirTNE ])[ktne ] = mfccc;
			//(D.f[dirTSW ])[ktsw ] = mfaac;
			//(D.f[dirTSE ])[ktse ] = mfcac;
			//(D.f[dirTNW ])[ktnw ] = mfacc;
			//(D.f[dirBNE ])[kbne ] = mfcca;
			//(D.f[dirBSW ])[kbsw ] = mfaaa;
			//(D.f[dirBSE ])[kbse ] = mfcaa;
			//(D.f[dirBNW ])[kbnw ] = mfaca;

      //(D.f[dirE   ])[ke   ] = fE ;  //f1_E ;   //fW;    //fE ;  
      //(D.f[dirW   ])[kw   ] = fW ;  //f1_W ;   //fE;    //fW ;  
      //(D.f[dirN   ])[kn   ] = fN ;  //f1_N ;   //fS;    //fN ;  
      //(D.f[dirS   ])[ks   ] = fS ;  //f1_S ;   //fN;    //fS ;  
      //(D.f[dirT   ])[kt   ] = fT ;  //f1_T ;   //fB;    //fT ;  
      //(D.f[dirB   ])[kb   ] = fB ;  //f1_B ;   //fT;    //fB ;  
      //(D.f[dirNE  ])[kne  ] = fNE;  //f1_NE;   //fSW;   //fNE;  
      //(D.f[dirSW  ])[ksw  ] = fSW;  //f1_SW;   //fNE;   //fSW;  
      //(D.f[dirSE  ])[kse  ] = fSE;  //f1_SE;   //fNW;   //fSE;  
      //(D.f[dirNW  ])[knw  ] = fNW;  //f1_NW;   //fSE;   //fNW;  
      //(D.f[dirTE  ])[kte  ] = fTE;  //f1_TE;   //fBW;   //fTE;  
      //(D.f[dirBW  ])[kbw  ] = fBW;  //f1_BW;   //fTE;   //fBW;  
      //(D.f[dirBE  ])[kbe  ] = fBE;  //f1_BE;   //fTW;   //fBE;  
      //(D.f[dirTW  ])[ktw  ] = fTW;  //f1_TW;   //fBE;   //fTW;  
      //(D.f[dirTN  ])[ktn  ] = fTN;  //f1_TN;   //fBS;   //fTN;  
      //(D.f[dirBS  ])[kbs  ] = fBS;  //f1_BS;   //fTN;   //fBS;  
      //(D.f[dirBN  ])[kbn  ] = fBN;  //f1_BN;   //fTS;   //fBN;  
      //(D.f[dirTS  ])[kts  ] = fTS;  //f1_TS;   //fBN;   //fTS;  
      //(D.f[dirZERO])[kzero] = fZERO;//f1_ZERO; //fZERO; //fZERO;
      //(D.f[dirTNE ])[ktne ] = fTNE; //f1_TNE;  //fBSW;  //fTNE; 
      //(D.f[dirBSW ])[kbsw ] = fBSW; //f1_BSW;  //fTNE;  //fBSW; 
      //(D.f[dirBNE ])[kbne ] = fBNE; //f1_BNE;  //fTSW;  //fBNE; 
      //(D.f[dirTSW ])[ktsw ] = fTSW; //f1_TSW;  //fBNE;  //fTSW; 
      //(D.f[dirTSE ])[ktse ] = fTSE; //f1_TSE;  //fBNW;  //fTSE; 
      //(D.f[dirBNW ])[kbnw ] = fBNW; //f1_BNW;  //fTSE;  //fBNW; 
      //(D.f[dirBSE ])[kbse ] = fBSE; //f1_BSE;  //fTNW;  //fBSE; 
      //(D.f[dirTNW ])[ktnw ] = fTNW; //f1_TNW;  //fBSE;  //fTNW; 
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceZero27(	 doubflo* DD, 
												 int* k_Q, 
												 int kQ, 
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
												 bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //__syncthreads();
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      (D.f[dirE   ])[ke   ] =zero;
      (D.f[dirW   ])[kw   ] =zero;
      (D.f[dirN   ])[kn   ] =zero;
      (D.f[dirS   ])[ks   ] =zero;
      (D.f[dirT   ])[kt   ] =zero;
      (D.f[dirB   ])[kb   ] =zero;
      (D.f[dirNE  ])[kne  ] =zero;
      (D.f[dirSW  ])[ksw  ] =zero;
      (D.f[dirSE  ])[kse  ] =zero;
      (D.f[dirNW  ])[knw  ] =zero;
      (D.f[dirTE  ])[kte  ] =zero;
      (D.f[dirBW  ])[kbw  ] =zero;
      (D.f[dirBE  ])[kbe  ] =zero;
      (D.f[dirTW  ])[ktw  ] =zero;
      (D.f[dirTN  ])[ktn  ] =zero;
      (D.f[dirBS  ])[kbs  ] =zero;
      (D.f[dirBN  ])[kbn  ] =zero;
      (D.f[dirTS  ])[kts  ] =zero;
      (D.f[dirZERO])[kzero] =zero;
      (D.f[dirTNE ])[ktne ] =zero;
      (D.f[dirTSW ])[ktsw ] =zero;
      (D.f[dirTSE ])[ktse ] =zero;
      (D.f[dirTNW ])[ktnw ] =zero;
      (D.f[dirBNE ])[kbne ] =zero;
      (D.f[dirBSW ])[kbsw ] =zero;
      (D.f[dirBSE ])[kbse ] =zero;
      (D.f[dirBNW ])[kbnw ] =zero;
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDeviceFake27(	 doubflo* rhoBC,
												 doubflo* DD, 
												 int* k_Q, 
												 int* k_N, 
												 int kQ, 
												 doubflo om1, 
												 unsigned int* neighborX,
												 unsigned int* neighborY,
												 unsigned int* neighborZ,
												 unsigned int size_Mat, 
												 bool evenOrOdd)
{
   ////////////////////////////////////////////////////////////////////////////////
   const unsigned  x = threadIdx.x;  // Globaler x-Index 
   const unsigned  y = blockIdx.x;   // Globaler y-Index 
   const unsigned  z = blockIdx.y;   // Globaler z-Index 

   const unsigned nx = blockDim.x;
   const unsigned ny = gridDim.x;

   const unsigned k = nx*(ny*z + y) + x;
   //////////////////////////////////////////////////////////////////////////

   if(k<kQ)
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
      if (evenOrOdd==false)
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
         D.f[dirTNE ] = &DD[dirBSW *size_Mat];
         D.f[dirTSW ] = &DD[dirBNE *size_Mat];
         D.f[dirTSE ] = &DD[dirBNW *size_Mat];
         D.f[dirTNW ] = &DD[dirBSE *size_Mat];
         D.f[dirBNE ] = &DD[dirTSW *size_Mat];
         D.f[dirBSW ] = &DD[dirTNE *size_Mat];
         D.f[dirBSE ] = &DD[dirTNW *size_Mat];
         D.f[dirBNW ] = &DD[dirTSE *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      doubflo        f1_E,f1_W,f1_N,f1_S,f1_T,f1_B,f1_NE,f1_SW,f1_SE,f1_NW,f1_TE,f1_BW,f1_BE,f1_TW,f1_TN,f1_BS,f1_BN,f1_TS,f1_ZERO,
         f1_TNE,f1_TSW,f1_TSE,f1_TNW,f1_BNE,f1_BSW,f1_BSE,f1_BNW;

      f1_W    = (D.f[dirE   ])[k1e   ];
      f1_E    = (D.f[dirW   ])[k1w   ];
      f1_S    = (D.f[dirN   ])[k1n   ];
      f1_N    = (D.f[dirS   ])[k1s   ];
      f1_B    = (D.f[dirT   ])[k1t   ];
      f1_T    = (D.f[dirB   ])[k1b   ];
      f1_SW   = (D.f[dirNE  ])[k1ne  ];
      f1_NE   = (D.f[dirSW  ])[k1sw  ];
      f1_NW   = (D.f[dirSE  ])[k1se  ];
      f1_SE   = (D.f[dirNW  ])[k1nw  ];
      f1_BW   = (D.f[dirTE  ])[k1te  ];
      f1_TE   = (D.f[dirBW  ])[k1bw  ];
      f1_TW   = (D.f[dirBE  ])[k1be  ];
      f1_BE   = (D.f[dirTW  ])[k1tw  ];
      f1_BS   = (D.f[dirTN  ])[k1tn  ];
      f1_TN   = (D.f[dirBS  ])[k1bs  ];
      f1_TS   = (D.f[dirBN  ])[k1bn  ];
      f1_BN   = (D.f[dirTS  ])[k1ts  ];
      f1_ZERO = (D.f[dirZERO])[k1zero];
      f1_BSW  = (D.f[dirTNE ])[k1tne ];
      f1_BNE  = (D.f[dirTSW ])[k1tsw ];
      f1_BNW  = (D.f[dirTSE ])[k1tse ];
      f1_BSE  = (D.f[dirTNW ])[k1tnw ];
      f1_TSW  = (D.f[dirBNE ])[k1bne ];
      f1_TNE  = (D.f[dirBSW ])[k1bsw ];
      f1_TNW  = (D.f[dirBSE ])[k1bse ];
      f1_TSE  = (D.f[dirBNW ])[k1bnw ];

      ////////////////////////////////////////////////////////////////////////////////
      doubflo vx1, vx2, vx3;
      vx1    =  ((f1_TSE - f1_BNW) - (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                  ((f1_BE - f1_TW)   + (f1_TE - f1_BW))   + ((f1_SE - f1_NW)   + (f1_NE - f1_SW)) +
                  (f1_E - f1_W); 


      vx2    =   (-(f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) - (f1_TSW - f1_BNE)) +
                  ((f1_BN - f1_TS)   + (f1_TN - f1_BS))    + (-(f1_SE - f1_NW)  + (f1_NE - f1_SW)) +
                  (f1_N - f1_S); 

      vx3    =   ((f1_TSE - f1_BNW) + (f1_TNW - f1_BSE)) + ((f1_TNE - f1_BSW) + (f1_TSW - f1_BNE)) +
                  (-(f1_BN - f1_TS)  + (f1_TN - f1_BS))   + ((f1_TE - f1_BW)   - (f1_BE - f1_TW)) +
                  (f1_T - f1_B); 

      doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);
      //////////////////////////////////////////////////////////////////////////
      doubflo drho1    =  f1_ZERO+f1_E+f1_W+f1_N+f1_S+f1_T+f1_B+f1_NE+f1_SW+f1_SE+f1_NW+f1_TE+f1_BW+f1_BE+f1_TW+f1_TN+f1_BS+f1_BN+f1_TS+
         f1_TNE+f1_TSW+f1_TSE+f1_TNW+f1_BNE+f1_BSW+f1_BSE+f1_BNW;

	  //drho1 = (drho1 + rhoBC[k])/2.f;
	  drho1 = drho1 - rhoBC[k];

      __syncthreads();

      (D.f[dirE   ])[ke   ] = c2over27* (rhoBC[k]+three*(-vx1        )+c9over2*(-vx1        )*(-vx1        )-cu_sq);
      (D.f[dirW   ])[kw   ] = c2over27* (rhoBC[k]+three*( vx1        )+c9over2*( vx1        )*( vx1        )-cu_sq);
      (D.f[dirN   ])[kn   ] = c2over27* (rhoBC[k]+three*(    -vx2    )+c9over2*(    -vx2    )*(    -vx2    )-cu_sq);
      (D.f[dirS   ])[ks   ] = c2over27* (rhoBC[k]+three*(     vx2    )+c9over2*(     vx2    )*(     vx2    )-cu_sq);
      (D.f[dirT   ])[kt   ] = c2over27* (rhoBC[k]+three*(        -vx3)+c9over2*(        -vx3)*(        -vx3)-cu_sq);
      (D.f[dirB   ])[kb   ] = c2over27* (rhoBC[k]+three*(         vx3)+c9over2*(         vx3)*(         vx3)-cu_sq);
      (D.f[dirNE  ])[kne  ] = f1_SW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirSW  ])[ksw  ] = f1_NE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirSE  ])[kse  ] = f1_NW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirNW  ])[knw  ] = f1_SE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTE  ])[kte  ] = f1_BW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBW  ])[kbw  ] = f1_TE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBE  ])[kbe  ] = f1_TW  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTW  ])[ktw  ] = f1_BE  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTN  ])[ktn  ] = f1_BS  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBS  ])[kbs  ] = f1_TN  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBN  ])[kbn  ] = f1_TS  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTS  ])[kts  ] = f1_BN  -c1over54*drho1;	//  c1o100;  // zero;  //
      (D.f[dirZERO])[kzero] = f1_ZERO-c8over27*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTNE ])[ktne ] = f1_BSW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTSW ])[ktsw ] = f1_BNE -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTSE ])[ktse ] = f1_BNW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirTNW ])[ktnw ] = f1_BSE -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBNE ])[kbne ] = f1_TSW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBSW ])[kbsw ] = f1_TNE -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBSE ])[kbse ] = f1_TNW -c1over216*drho1;	//  c1o100;  // zero;  //
      (D.f[dirBNW ])[kbnw ] = f1_TSE -c1over216*drho1;  //  c1o100;  // zero;  //      
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







































//////////////////////////////////////////////////////////////////////////
extern "C" __global__ void QPressDevice27_IntBB(doubflo* rho,
												doubflo* DD, 
												int* k_Q, 
												doubflo* QQ,
												unsigned int sizeQ,
												int kQ, 
												doubflo om1, 
												unsigned int* neighborX,
												unsigned int* neighborY,
												unsigned int* neighborZ,
												unsigned int size_Mat, 
												bool evenOrOdd)
{
	Distributions27 D;
	if (evenOrOdd==true)
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
		D.f[dirTNE ] = &DD[dirBSW *size_Mat];
		D.f[dirTSW ] = &DD[dirBNE *size_Mat];
		D.f[dirTSE ] = &DD[dirBNW *size_Mat];
		D.f[dirTNW ] = &DD[dirBSE *size_Mat];
		D.f[dirBNE ] = &DD[dirTSW *size_Mat];
		D.f[dirBSW ] = &DD[dirTNE *size_Mat];
		D.f[dirBSE ] = &DD[dirTNW *size_Mat];
		D.f[dirBNW ] = &DD[dirTSE *size_Mat];
	}
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////

	if(k<kQ)
	{
		////////////////////////////////////////////////////////////////////////////////
		//doubflo VeloX = vx[k];
		//doubflo VeloY = vy[k];
		//doubflo VeloZ = vz[k]; //(16.0*(u0*2.0)*bbx*bby*(grid_nx-bbx)*(grid_ny-bby))/(grid_nx*grid_nx*grid_ny*grid_ny)
		////////////////////////////////////////////////////////////////////////////////
		doubflo *q_dirE,   *q_dirW,   *q_dirN,   *q_dirS,   *q_dirT,   *q_dirB, 
			*q_dirNE,  *q_dirSW,  *q_dirSE,  *q_dirNW,  *q_dirTE,  *q_dirBW,
			*q_dirBE,  *q_dirTW,  *q_dirTN,  *q_dirBS,  *q_dirBN,  *q_dirTS,
			*q_dirTNE, *q_dirTSW, *q_dirTSE, *q_dirTNW, *q_dirBNE, *q_dirBSW,
			*q_dirBSE, *q_dirBNW; 
		q_dirE   = &QQ[dirE   *sizeQ];
		q_dirW   = &QQ[dirW   *sizeQ];
		q_dirN   = &QQ[dirN   *sizeQ];
		q_dirS   = &QQ[dirS   *sizeQ];
		q_dirT   = &QQ[dirT   *sizeQ];
		q_dirB   = &QQ[dirB   *sizeQ];
		q_dirNE  = &QQ[dirNE  *sizeQ];
		q_dirSW  = &QQ[dirSW  *sizeQ];
		q_dirSE  = &QQ[dirSE  *sizeQ];
		q_dirNW  = &QQ[dirNW  *sizeQ];
		q_dirTE  = &QQ[dirTE  *sizeQ];
		q_dirBW  = &QQ[dirBW  *sizeQ];
		q_dirBE  = &QQ[dirBE  *sizeQ];
		q_dirTW  = &QQ[dirTW  *sizeQ];
		q_dirTN  = &QQ[dirTN  *sizeQ];
		q_dirBS  = &QQ[dirBS  *sizeQ];
		q_dirBN  = &QQ[dirBN  *sizeQ];
		q_dirTS  = &QQ[dirTS  *sizeQ];
		q_dirTNE = &QQ[dirTNE *sizeQ];
		q_dirTSW = &QQ[dirTSW *sizeQ];
		q_dirTSE = &QQ[dirTSE *sizeQ];
		q_dirTNW = &QQ[dirTNW *sizeQ];
		q_dirBNE = &QQ[dirBNE *sizeQ];
		q_dirBSW = &QQ[dirBSW *sizeQ];
		q_dirBSE = &QQ[dirBSE *sizeQ];
		q_dirBNW = &QQ[dirBNW *sizeQ];
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
		doubflo f_E,  f_W,  f_N,  f_S,  f_T,  f_B,   f_NE,  f_SW,  f_SE,  f_NW,  f_TE,  f_BW,  f_BE,
			f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

		f_W    = (D.f[dirE   ])[ke   ];
		f_E    = (D.f[dirW   ])[kw   ];
		f_S    = (D.f[dirN   ])[kn   ];
		f_N    = (D.f[dirS   ])[ks   ];
		f_B    = (D.f[dirT   ])[kt   ];
		f_T    = (D.f[dirB   ])[kb   ];
		f_SW   = (D.f[dirNE  ])[kne  ];
		f_NE   = (D.f[dirSW  ])[ksw  ];
		f_NW   = (D.f[dirSE  ])[kse  ];
		f_SE   = (D.f[dirNW  ])[knw  ];
		f_BW   = (D.f[dirTE  ])[kte  ];
		f_TE   = (D.f[dirBW  ])[kbw  ];
		f_TW   = (D.f[dirBE  ])[kbe  ];
		f_BE   = (D.f[dirTW  ])[ktw  ];
		f_BS   = (D.f[dirTN  ])[ktn  ];
		f_TN   = (D.f[dirBS  ])[kbs  ];
		f_TS   = (D.f[dirBN  ])[kbn  ];
		f_BN   = (D.f[dirTS  ])[kts  ];
		f_BSW  = (D.f[dirTNE ])[ktne ];
		f_BNE  = (D.f[dirTSW ])[ktsw ];
		f_BNW  = (D.f[dirTSE ])[ktse ];
		f_BSE  = (D.f[dirTNW ])[ktnw ];
		f_TSW  = (D.f[dirBNE ])[kbne ];
		f_TNE  = (D.f[dirBSW ])[kbsw ];
		f_TNW  = (D.f[dirBSE ])[kbse ];
		f_TSE  = (D.f[dirBNW ])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////
		doubflo vx1, vx2, vx3, drho, feq, q;
		drho   =  f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
			f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW + 
			f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[dirZERO])[kzero]); 

		vx1    = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
			((f_BE - f_TW)   + (f_TE - f_BW))   + ((f_SE - f_NW)   + (f_NE - f_SW)) +
			(f_E - f_W))/(one+drho); 


		vx2    =  ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
			((f_BN - f_TS)   + (f_TN - f_BS))    + (-(f_SE - f_NW)  + (f_NE - f_SW)) +
			(f_N - f_S))/(one+drho); 

		vx3    =  (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
			(-(f_BN - f_TS)  + (f_TN - f_BS))   + ((f_TE - f_BW)   - (f_BE - f_TW)) +
			(f_T - f_B))/(one+drho); 

		doubflo cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3);

		//////////////////////////////////////////////////////////////////////////
		if (evenOrOdd==false)
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
			D.f[dirTNE ] = &DD[dirBSW *size_Mat];
			D.f[dirTSW ] = &DD[dirBNE *size_Mat];
			D.f[dirTSE ] = &DD[dirBNW *size_Mat];
			D.f[dirTNW ] = &DD[dirBSE *size_Mat];
			D.f[dirBNE ] = &DD[dirTSW *size_Mat];
			D.f[dirBSW ] = &DD[dirTNE *size_Mat];
			D.f[dirBSE ] = &DD[dirTNW *size_Mat];
			D.f[dirBNW ] = &DD[dirTSE *size_Mat];
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Test
		//(D.f[dirZERO])[k]=c1o10;
		doubflo rhoDiff = drho - rho[k];
		doubflo VeloX = vx1;
		doubflo VeloY = vx2;
		doubflo VeloZ = vx3;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		q = q_dirE[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+c9over2*( vx1        )*( vx1        )-cu_sq); 
			(D.f[dirW])[kw]=(one-q)/(one+q)*(f_E-f_W+(f_E+f_W-two*feq*om1)/(one-om1))*c1o2+(q*(f_E+f_W)-c2over27*(rhoDiff + six*( VeloX     )))/(one+q);
		}

		q = q_dirW[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+c9over2*(-vx1        )*(-vx1        )-cu_sq); 
			(D.f[dirE])[ke]=(one-q)/(one+q)*(f_W-f_E+(f_W+f_E-two*feq*om1)/(one-om1))*c1o2+(q*(f_W+f_E)-c2over27*(rhoDiff + six*(-VeloX     )))/(one+q);
		}

		q = q_dirN[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+c9over2*(     vx2    )*(     vx2    )-cu_sq); 
			(D.f[dirS])[ks]=(one-q)/(one+q)*(f_N-f_S+(f_N+f_S-two*feq*om1)/(one-om1))*c1o2+(q*(f_N+f_S)-c2over27*(rhoDiff + six*( VeloY     )))/(one+q);
		}

		q = q_dirS[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+c9over2*(    -vx2    )*(    -vx2    )-cu_sq); 
			(D.f[dirN])[kn]=(one-q)/(one+q)*(f_S-f_N+(f_S+f_N-two*feq*om1)/(one-om1))*c1o2+(q*(f_S+f_N)-c2over27*(rhoDiff + six*(-VeloY     )))/(one+q);
		}

		q = q_dirT[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+c9over2*(         vx3)*(         vx3)-cu_sq); 
			(D.f[dirB])[kb]=(one-q)/(one+q)*(f_T-f_B+(f_T+f_B-two*feq*om1)/(one-om1))*c1o2+(q*(f_T+f_B)-c2over27*(rhoDiff + six*( VeloZ     )))/(one+q);
		}

		q = q_dirB[k];
		if (q>=zero && q<=one)
		{
			feq=c2over27* (drho+c9over2*(        -vx3)*(        -vx3)-cu_sq); 
			(D.f[dirT])[kt]=(one-q)/(one+q)*(f_B-f_T+(f_B+f_T-two*feq*om1)/(one-om1))*c1o2+(q*(f_B+f_T)-c2over27*(rhoDiff + six*(-VeloZ     )))/(one+q);
		}

		q = q_dirNE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*( vx1+vx2    )*( vx1+vx2    )-cu_sq); 
			(D.f[dirSW])[ksw]=(one-q)/(one+q)*(f_NE-f_SW+(f_NE+f_SW-two*feq*om1)/(one-om1))*c1o2+(q*(f_NE+f_SW)-c1over54*(rhoDiff + six*(VeloX+VeloY)))/(one+q);
		}

		q = q_dirSW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(-vx1-vx2    )*(-vx1-vx2    )-cu_sq); 
			(D.f[dirNE])[kne]=(one-q)/(one+q)*(f_SW-f_NE+(f_SW+f_NE-two*feq*om1)/(one-om1))*c1o2+(q*(f_SW+f_NE)-c1over54*(rhoDiff + six*(-VeloX-VeloY)))/(one+q);
		}

		q = q_dirSE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*( vx1-vx2    )*( vx1-vx2    )-cu_sq); 
			(D.f[dirNW])[knw]=(one-q)/(one+q)*(f_SE-f_NW+(f_SE+f_NW-two*feq*om1)/(one-om1))*c1o2+(q*(f_SE+f_NW)-c1over54*(rhoDiff + six*( VeloX-VeloY)))/(one+q);
		}

		q = q_dirNW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(-vx1+vx2    )*(-vx1+vx2    )-cu_sq); 
			(D.f[dirSE])[kse]=(one-q)/(one+q)*(f_NW-f_SE+(f_NW+f_SE-two*feq*om1)/(one-om1))*c1o2+(q*(f_NW+f_SE)-c1over54*(rhoDiff + six*(-VeloX+VeloY)))/(one+q);
		}

		q = q_dirTE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*( vx1    +vx3)*( vx1    +vx3)-cu_sq); 
			(D.f[dirBW])[kbw]=(one-q)/(one+q)*(f_TE-f_BW+(f_TE+f_BW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TE+f_BW)-c1over54*(rhoDiff + six*( VeloX+VeloZ)))/(one+q);
		}

		q = q_dirBW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(-vx1    -vx3)*(-vx1    -vx3)-cu_sq); 
			(D.f[dirTE])[kte]=(one-q)/(one+q)*(f_BW-f_TE+(f_BW+f_TE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BW+f_TE)-c1over54*(rhoDiff + six*(-VeloX-VeloZ)))/(one+q);
		}

		q = q_dirBE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*( vx1    -vx3)*( vx1    -vx3)-cu_sq); 
			(D.f[dirTW])[ktw]=(one-q)/(one+q)*(f_BE-f_TW+(f_BE+f_TW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BE+f_TW)-c1over54*(rhoDiff + six*( VeloX-VeloZ)))/(one+q);
		}

		q = q_dirTW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(-vx1    +vx3)*(-vx1    +vx3)-cu_sq); 
			(D.f[dirBE])[kbe]=(one-q)/(one+q)*(f_TW-f_BE+(f_TW+f_BE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TW+f_BE)-c1over54*(rhoDiff + six*(-VeloX+VeloZ)))/(one+q);
		}

		q = q_dirTN[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(     vx2+vx3)*(     vx2+vx3)-cu_sq); 
			(D.f[dirBS])[kbs]=(one-q)/(one+q)*(f_TN-f_BS+(f_TN+f_BS-two*feq*om1)/(one-om1))*c1o2+(q*(f_TN+f_BS)-c1over54*(rhoDiff + six*( VeloY+VeloZ)))/(one+q);
		}

		q = q_dirBS[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(    -vx2-vx3)*(    -vx2-vx3)-cu_sq); 
			(D.f[dirTN])[ktn]=(one-q)/(one+q)*(f_BS-f_TN+(f_BS+f_TN-two*feq*om1)/(one-om1))*c1o2+(q*(f_BS+f_TN)-c1over54*(rhoDiff + six*( -VeloY-VeloZ)))/(one+q);
		}

		q = q_dirBN[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(     vx2-vx3)*(     vx2-vx3)-cu_sq); 
			(D.f[dirTS])[kts]=(one-q)/(one+q)*(f_BN-f_TS+(f_BN+f_TS-two*feq*om1)/(one-om1))*c1o2+(q*(f_BN+f_TS)-c1over54*(rhoDiff + six*( VeloY-VeloZ)))/(one+q);
		}

		q = q_dirTS[k];
		if (q>=zero && q<=one)
		{
			feq=c1over54* (drho+c9over2*(    -vx2+vx3)*(    -vx2+vx3)-cu_sq); 
			(D.f[dirBN])[kbn]=(one-q)/(one+q)*(f_TS-f_BN+(f_TS+f_BN-two*feq*om1)/(one-om1))*c1o2+(q*(f_TS+f_BN)-c1over54*(rhoDiff + six*( -VeloY+VeloZ)))/(one+q);
		}

		q = q_dirTNE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cu_sq); 
			(D.f[dirBSW])[kbsw]=(one-q)/(one+q)*(f_TNE-f_BSW+(f_TNE+f_BSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNE+f_BSW)-c1over216*(rhoDiff + six*( VeloX+VeloY+VeloZ)))/(one+q);
		}

		q = q_dirBSW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cu_sq); 
			(D.f[dirTNE])[ktne]=(one-q)/(one+q)*(f_BSW-f_TNE+(f_BSW+f_TNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSW+f_TNE)-c1over216*(rhoDiff + six*(-VeloX-VeloY-VeloZ)))/(one+q);
		}

		q = q_dirBNE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cu_sq); 
			(D.f[dirTSW])[ktsw]=(one-q)/(one+q)*(f_BNE-f_TSW+(f_BNE+f_TSW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNE+f_TSW)-c1over216*(rhoDiff + six*( VeloX+VeloY-VeloZ)))/(one+q);
		}

		q = q_dirTSW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cu_sq); 
			(D.f[dirBNE])[kbne]=(one-q)/(one+q)*(f_TSW-f_BNE+(f_TSW+f_BNE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSW+f_BNE)-c1over216*(rhoDiff + six*(-VeloX-VeloY+VeloZ)))/(one+q);
		}

		q = q_dirTSE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cu_sq); 
			(D.f[dirBNW])[kbnw]=(one-q)/(one+q)*(f_TSE-f_BNW+(f_TSE+f_BNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_TSE+f_BNW)-c1over216*(rhoDiff + six*( VeloX-VeloY+VeloZ)))/(one+q);
		}

		q = q_dirBNW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cu_sq); 
			(D.f[dirTSE])[ktse]=(one-q)/(one+q)*(f_BNW-f_TSE+(f_BNW+f_TSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_BNW+f_TSE)-c1over216*(rhoDiff + six*(-VeloX+VeloY-VeloZ)))/(one+q);
		}

		q = q_dirBSE[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cu_sq); 
			(D.f[dirTNW])[ktnw]=(one-q)/(one+q)*(f_BSE-f_TNW+(f_BSE+f_TNW-two*feq*om1)/(one-om1))*c1o2+(q*(f_BSE+f_TNW)-c1over216*(rhoDiff + six*( VeloX-VeloY-VeloZ)))/(one+q);
		}

		q = q_dirTNW[k];
		if (q>=zero && q<=one)
		{
			feq=c1over216*(drho+c9over2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cu_sq); 
			(D.f[dirBSE])[kbse]=(one-q)/(one+q)*(f_TNW-f_BSE+(f_TNW+f_BSE-two*feq*om1)/(one-om1))*c1o2+(q*(f_TNW+f_BSE)-c1over216*(rhoDiff + six*(-VeloX+VeloY+VeloZ)))/(one+q);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


