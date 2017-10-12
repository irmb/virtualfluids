/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void QPressDeviceOld27(doubflo* rhoBC,
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

	  drho1 = drho1 - rhoBC[k];
      //////////////////////////////////////////////////////////////////////////

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
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

