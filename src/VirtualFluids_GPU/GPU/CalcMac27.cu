/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcMacCompSP27( doubflo* vxD,
											  doubflo* vyD,
											  doubflo* vzD,
											  doubflo* rhoD,
											  doubflo* pressD,
											  unsigned int* geoD,
											  unsigned int* neighborX,
											  unsigned int* neighborY,
											  unsigned int* neighborZ,
											  unsigned int size_Mat,
											  doubflo* DD,
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

   if(k<size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
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
      //////////////////////////////////////////////////////////////////////////
      pressD[k] = zero;
	  rhoD[k]   = zero;
	  vxD[k]    = zero;
	  vyD[k]    = zero;
	  vzD[k]    = zero;

      if(geoD[k] == GEO_FLUID)
      {
         rhoD[k]    =   (D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirZERO])[kzero]+ 
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw];

         vxD[k]     =  ((D.f[dirE   ])[ke  ]- (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]- (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]- (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]- (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]+ 
						(D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]) / (one + rhoD[k]);

         vyD[k]     =  ((D.f[dirN   ])[kn  ]- (D.f[dirS   ])[ks  ]+
                        (D.f[dirNE  ])[kne ]- (D.f[dirSW  ])[ksw ]-
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]- (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]- (D.f[dirTSW ])[ktsw]- 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw]) / (one + rhoD[k]);

         vzD[k]     =  ((D.f[dirT   ])[kt  ]- (D.f[dirB   ])[kb  ]+
                        (D.f[dirTE  ])[kte ]- (D.f[dirBW  ])[kbw ]-
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]- (D.f[dirBS  ])[kbs ]-
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ]+
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]- 
                        (D.f[dirBNE ])[kbne]- (D.f[dirBSW ])[kbsw]- 
                        (D.f[dirBSE ])[kbse]- (D.f[dirBNW ])[kbnw]) / (one + rhoD[k]);

         pressD[k]  =  ((D.f[dirE   ])[ke  ]+ (D.f[dirW   ])[kw  ]+ 
                        (D.f[dirN   ])[kn  ]+ (D.f[dirS   ])[ks  ]+
                        (D.f[dirT   ])[kt  ]+ (D.f[dirB   ])[kb  ]+
                        two*(
                        (D.f[dirNE  ])[kne ]+ (D.f[dirSW  ])[ksw ]+
                        (D.f[dirSE  ])[kse ]+ (D.f[dirNW  ])[knw ]+
                        (D.f[dirTE  ])[kte ]+ (D.f[dirBW  ])[kbw ]+
                        (D.f[dirBE  ])[kbe ]+ (D.f[dirTW  ])[ktw ]+
                        (D.f[dirTN  ])[ktn ]+ (D.f[dirBS  ])[kbs ]+
                        (D.f[dirBN  ])[kbn ]+ (D.f[dirTS  ])[kts ])+
                        three*(
                        (D.f[dirTNE ])[ktne]+ (D.f[dirTSW ])[ktsw]+ 
                        (D.f[dirTSE ])[ktse]+ (D.f[dirTNW ])[ktnw]+ 
                        (D.f[dirBNE ])[kbne]+ (D.f[dirBSW ])[kbsw]+ 
                        (D.f[dirBSE ])[kbse]+ (D.f[dirBNW ])[kbnw])-
                        rhoD[k]-(vxD[k] * vxD[k] + vyD[k] * vyD[k] + vzD[k] * vzD[k]) * (one+rhoD[k])) * c1o2+rhoD[k]; // times zero for incompressible case   
         //achtung op hart gesetzt Annahme op = 1 ;                                                      ^^^^(1.0/op-0.5)=0.5

      }
   }
}




