/* Device code */
#include "LBM/D3Q27.h"
#include "GPU/constant.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void getSendFsPost27(doubflo* DD,
										   doubflo* bufferFs,
										   int* sendIndex,
                                           int buffmax,
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

   if(k<buffmax)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //set index
      unsigned int kIndex  = sendIndex[k];
      unsigned int kzero   = kIndex;
      unsigned int ke      = kIndex;
      unsigned int kw      = neighborX[kIndex];
      unsigned int kn      = kIndex;
      unsigned int ks      = neighborY[kIndex];
      unsigned int kt      = kIndex;
      unsigned int kb      = neighborZ[kIndex];
      unsigned int ksw     = neighborY[kw];
      unsigned int kne     = kIndex;
      unsigned int kse     = ks;
      unsigned int knw     = kw;
      unsigned int kbw     = neighborZ[kw];
      unsigned int kte     = kIndex;
      unsigned int kbe     = kb;
      unsigned int ktw     = kw;
      unsigned int kbs     = neighborZ[ks];
      unsigned int ktn     = kIndex;
      unsigned int kbn     = kb;
      unsigned int kts     = ks;
      unsigned int ktse    = ks;
      unsigned int kbnw    = kbw;
      unsigned int ktnw    = kw;
      unsigned int kbse    = kbs;
      unsigned int ktsw    = ksw;
      unsigned int kbne    = kb;
      unsigned int ktne    = kIndex;
      unsigned int kbsw    = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Fs
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dirE   ] = &bufferFs[dirE   *buffmax];
      Dbuff.f[dirW   ] = &bufferFs[dirW   *buffmax];
      Dbuff.f[dirN   ] = &bufferFs[dirN   *buffmax];
      Dbuff.f[dirS   ] = &bufferFs[dirS   *buffmax];
      Dbuff.f[dirT   ] = &bufferFs[dirT   *buffmax];
      Dbuff.f[dirB   ] = &bufferFs[dirB   *buffmax];
      Dbuff.f[dirNE  ] = &bufferFs[dirNE  *buffmax];
      Dbuff.f[dirSW  ] = &bufferFs[dirSW  *buffmax];
      Dbuff.f[dirSE  ] = &bufferFs[dirSE  *buffmax];
      Dbuff.f[dirNW  ] = &bufferFs[dirNW  *buffmax];
      Dbuff.f[dirTE  ] = &bufferFs[dirTE  *buffmax];
      Dbuff.f[dirBW  ] = &bufferFs[dirBW  *buffmax];
      Dbuff.f[dirBE  ] = &bufferFs[dirBE  *buffmax];
      Dbuff.f[dirTW  ] = &bufferFs[dirTW  *buffmax];
      Dbuff.f[dirTN  ] = &bufferFs[dirTN  *buffmax];
      Dbuff.f[dirBS  ] = &bufferFs[dirBS  *buffmax];
      Dbuff.f[dirBN  ] = &bufferFs[dirBN  *buffmax];
      Dbuff.f[dirTS  ] = &bufferFs[dirTS  *buffmax];
      Dbuff.f[dirZERO] = &bufferFs[dirZERO*buffmax];
      Dbuff.f[dirTNE ] = &bufferFs[dirTNE *buffmax];
      Dbuff.f[dirTSW ] = &bufferFs[dirTSW *buffmax];
      Dbuff.f[dirTSE ] = &bufferFs[dirTSE *buffmax];
      Dbuff.f[dirTNW ] = &bufferFs[dirTNW *buffmax];
      Dbuff.f[dirBNE ] = &bufferFs[dirBNE *buffmax];
      Dbuff.f[dirBSW ] = &bufferFs[dirBSW *buffmax];
      Dbuff.f[dirBSE ] = &bufferFs[dirBSE *buffmax];
      Dbuff.f[dirBNW ] = &bufferFs[dirBNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      //(Dbuff.f[dirE   ])[k] = (D.f[dirE   ])[ke   ];
      //(Dbuff.f[dirW   ])[k] = (D.f[dirW   ])[kw   ];
      //(Dbuff.f[dirN   ])[k] = (D.f[dirN   ])[kn   ];
      //(Dbuff.f[dirS   ])[k] = (D.f[dirS   ])[ks   ];
      //(Dbuff.f[dirT   ])[k] = (D.f[dirT   ])[kt   ];
      //(Dbuff.f[dirB   ])[k] = (D.f[dirB   ])[kb   ];
      //(Dbuff.f[dirNE  ])[k] = (D.f[dirNE  ])[kne  ];
      //(Dbuff.f[dirSW  ])[k] = (D.f[dirSW  ])[ksw  ];
      //(Dbuff.f[dirSE  ])[k] = (D.f[dirSE  ])[kse  ];
      //(Dbuff.f[dirNW  ])[k] = (D.f[dirNW  ])[knw  ];
      //(Dbuff.f[dirTE  ])[k] = (D.f[dirTE  ])[kte  ];
      //(Dbuff.f[dirBW  ])[k] = (D.f[dirBW  ])[kbw  ];
      //(Dbuff.f[dirBE  ])[k] = (D.f[dirBE  ])[kbe  ];
      //(Dbuff.f[dirTW  ])[k] = (D.f[dirTW  ])[ktw  ];
      //(Dbuff.f[dirTN  ])[k] = (D.f[dirTN  ])[ktn  ];
      //(Dbuff.f[dirBS  ])[k] = (D.f[dirBS  ])[kbs  ];
      //(Dbuff.f[dirBN  ])[k] = (D.f[dirBN  ])[kbn  ];
      //(Dbuff.f[dirTS  ])[k] = (D.f[dirTS  ])[kts  ];
      //(Dbuff.f[dirZERO])[k] = (D.f[dirZERO])[kzero];
      //(Dbuff.f[dirTNE ])[k] = (D.f[dirTNE ])[ktne ];
      //(Dbuff.f[dirTSW ])[k] = (D.f[dirTSW ])[ktsw ];
      //(Dbuff.f[dirTSE ])[k] = (D.f[dirTSE ])[ktse ];
      //(Dbuff.f[dirTNW ])[k] = (D.f[dirTNW ])[ktnw ];
      //(Dbuff.f[dirBNE ])[k] = (D.f[dirBNE ])[kbne ];
      //(Dbuff.f[dirBSW ])[k] = (D.f[dirBSW ])[kbsw ];
      //(Dbuff.f[dirBSE ])[k] = (D.f[dirBSE ])[kbse ];
      //(Dbuff.f[dirBNW ])[k] = (D.f[dirBNW ])[kbnw ];
      (Dbuff.f[dirE   ])[k] = (D.f[dirW   ])[kw   ];
      (Dbuff.f[dirW   ])[k] = (D.f[dirE   ])[ke   ];
      (Dbuff.f[dirN   ])[k] = (D.f[dirS   ])[ks   ];
      (Dbuff.f[dirS   ])[k] = (D.f[dirN   ])[kn   ];
      (Dbuff.f[dirT   ])[k] = (D.f[dirB   ])[kb   ];
      (Dbuff.f[dirB   ])[k] = (D.f[dirT   ])[kt   ];
      (Dbuff.f[dirNE  ])[k] = (D.f[dirSW  ])[ksw  ];
      (Dbuff.f[dirSW  ])[k] = (D.f[dirNE  ])[kne  ];
      (Dbuff.f[dirSE  ])[k] = (D.f[dirNW  ])[knw  ];
      (Dbuff.f[dirNW  ])[k] = (D.f[dirSE  ])[kse  ];
      (Dbuff.f[dirTE  ])[k] = (D.f[dirBW  ])[kbw  ];
      (Dbuff.f[dirBW  ])[k] = (D.f[dirTE  ])[kte  ];
      (Dbuff.f[dirBE  ])[k] = (D.f[dirTW  ])[ktw  ];
      (Dbuff.f[dirTW  ])[k] = (D.f[dirBE  ])[kbe  ];
      (Dbuff.f[dirTN  ])[k] = (D.f[dirBS  ])[kbs  ];
      (Dbuff.f[dirBS  ])[k] = (D.f[dirTN  ])[ktn  ];
      (Dbuff.f[dirBN  ])[k] = (D.f[dirTS  ])[kts  ];
      (Dbuff.f[dirTS  ])[k] = (D.f[dirBN  ])[kbn  ];
      (Dbuff.f[dirZERO])[k] = (D.f[dirZERO])[kzero];
      (Dbuff.f[dirTNE ])[k] = (D.f[dirBSW ])[kbsw ];
      (Dbuff.f[dirTSW ])[k] = (D.f[dirBNE ])[kbne ];
      (Dbuff.f[dirTSE ])[k] = (D.f[dirBNW ])[kbnw ];
      (Dbuff.f[dirTNW ])[k] = (D.f[dirBSE ])[kbse ];
      (Dbuff.f[dirBNE ])[k] = (D.f[dirTSW ])[ktsw ];
      (Dbuff.f[dirBSW ])[k] = (D.f[dirTNE ])[ktne ];
      (Dbuff.f[dirBSE ])[k] = (D.f[dirTNW ])[ktnw ];
      (Dbuff.f[dirBNW ])[k] = (D.f[dirTSE ])[ktse ];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void setRecvFsPost27(doubflo* DD,
										   doubflo* bufferFs,
										   int* recvIndex,
                                           int buffmax,
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

   if(k<buffmax)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //set index
      unsigned int kIndex  = recvIndex[k];
      unsigned int kzero   = kIndex;
      unsigned int ke      = kIndex;
      unsigned int kw      = neighborX[kIndex];
      unsigned int kn      = kIndex;
      unsigned int ks      = neighborY[kIndex];
      unsigned int kt      = kIndex;
      unsigned int kb      = neighborZ[kIndex];
      unsigned int ksw     = neighborY[kw];
      unsigned int kne     = kIndex;
      unsigned int kse     = ks;
      unsigned int knw     = kw;
      unsigned int kbw     = neighborZ[kw];
      unsigned int kte     = kIndex;
      unsigned int kbe     = kb;
      unsigned int ktw     = kw;
      unsigned int kbs     = neighborZ[ks];
      unsigned int ktn     = kIndex;
      unsigned int kbn     = kb;
      unsigned int kts     = ks;
      unsigned int ktse    = ks;
      unsigned int kbnw    = kbw;
      unsigned int ktnw    = kw;
      unsigned int kbse    = kbs;
      unsigned int ktsw    = ksw;
      unsigned int kbne    = kb;
      unsigned int ktne    = kIndex;
      unsigned int kbsw    = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Fs
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dirE   ] = &bufferFs[dirE   *buffmax];
      Dbuff.f[dirW   ] = &bufferFs[dirW   *buffmax];
      Dbuff.f[dirN   ] = &bufferFs[dirN   *buffmax];
      Dbuff.f[dirS   ] = &bufferFs[dirS   *buffmax];
      Dbuff.f[dirT   ] = &bufferFs[dirT   *buffmax];
      Dbuff.f[dirB   ] = &bufferFs[dirB   *buffmax];
      Dbuff.f[dirNE  ] = &bufferFs[dirNE  *buffmax];
      Dbuff.f[dirSW  ] = &bufferFs[dirSW  *buffmax];
      Dbuff.f[dirSE  ] = &bufferFs[dirSE  *buffmax];
      Dbuff.f[dirNW  ] = &bufferFs[dirNW  *buffmax];
      Dbuff.f[dirTE  ] = &bufferFs[dirTE  *buffmax];
      Dbuff.f[dirBW  ] = &bufferFs[dirBW  *buffmax];
      Dbuff.f[dirBE  ] = &bufferFs[dirBE  *buffmax];
      Dbuff.f[dirTW  ] = &bufferFs[dirTW  *buffmax];
      Dbuff.f[dirTN  ] = &bufferFs[dirTN  *buffmax];
      Dbuff.f[dirBS  ] = &bufferFs[dirBS  *buffmax];
      Dbuff.f[dirBN  ] = &bufferFs[dirBN  *buffmax];
      Dbuff.f[dirTS  ] = &bufferFs[dirTS  *buffmax];
      Dbuff.f[dirZERO] = &bufferFs[dirZERO*buffmax];
      Dbuff.f[dirTNE ] = &bufferFs[dirTNE *buffmax];
      Dbuff.f[dirTSW ] = &bufferFs[dirTSW *buffmax];
      Dbuff.f[dirTSE ] = &bufferFs[dirTSE *buffmax];
      Dbuff.f[dirTNW ] = &bufferFs[dirTNW *buffmax];
      Dbuff.f[dirBNE ] = &bufferFs[dirBNE *buffmax];
      Dbuff.f[dirBSW ] = &bufferFs[dirBSW *buffmax];
      Dbuff.f[dirBSE ] = &bufferFs[dirBSE *buffmax];
      Dbuff.f[dirBNW ] = &bufferFs[dirBNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      //(D.f[dirE   ])[ke   ] = (Dbuff.f[dirE   ])[k];
      //(D.f[dirW   ])[kw   ] = (Dbuff.f[dirW   ])[k];
      //(D.f[dirN   ])[kn   ] = (Dbuff.f[dirN   ])[k];
      //(D.f[dirS   ])[ks   ] = (Dbuff.f[dirS   ])[k];
      //(D.f[dirT   ])[kt   ] = (Dbuff.f[dirT   ])[k];
      //(D.f[dirB   ])[kb   ] = (Dbuff.f[dirB   ])[k];
      //(D.f[dirNE  ])[kne  ] = (Dbuff.f[dirNE  ])[k];
      //(D.f[dirSW  ])[ksw  ] = (Dbuff.f[dirSW  ])[k];
      //(D.f[dirSE  ])[kse  ] = (Dbuff.f[dirSE  ])[k];
      //(D.f[dirNW  ])[knw  ] = (Dbuff.f[dirNW  ])[k];
      //(D.f[dirTE  ])[kte  ] = (Dbuff.f[dirTE  ])[k];
      //(D.f[dirBW  ])[kbw  ] = (Dbuff.f[dirBW  ])[k];
      //(D.f[dirBE  ])[kbe  ] = (Dbuff.f[dirBE  ])[k];
      //(D.f[dirTW  ])[ktw  ] = (Dbuff.f[dirTW  ])[k];
      //(D.f[dirTN  ])[ktn  ] = (Dbuff.f[dirTN  ])[k];
      //(D.f[dirBS  ])[kbs  ] = (Dbuff.f[dirBS  ])[k];
      //(D.f[dirBN  ])[kbn  ] = (Dbuff.f[dirBN  ])[k];
      //(D.f[dirTS  ])[kts  ] = (Dbuff.f[dirTS  ])[k];
      //(D.f[dirZERO])[kzero] = (Dbuff.f[dirZERO])[k];
      //(D.f[dirTNE ])[ktne ] = (Dbuff.f[dirTNE ])[k];
      //(D.f[dirTSW ])[ktsw ] = (Dbuff.f[dirTSW ])[k];
      //(D.f[dirTSE ])[ktse ] = (Dbuff.f[dirTSE ])[k];
      //(D.f[dirTNW ])[ktnw ] = (Dbuff.f[dirTNW ])[k];
      //(D.f[dirBNE ])[kbne ] = (Dbuff.f[dirBNE ])[k];
      //(D.f[dirBSW ])[kbsw ] = (Dbuff.f[dirBSW ])[k];
      //(D.f[dirBSE ])[kbse ] = (Dbuff.f[dirBSE ])[k];
      //(D.f[dirBNW ])[kbnw ] = (Dbuff.f[dirBNW ])[k];
      (D.f[dirW   ])[kw   ] = (Dbuff.f[dirE   ])[k];
      (D.f[dirE   ])[ke   ] = (Dbuff.f[dirW   ])[k];
      (D.f[dirS   ])[ks   ] = (Dbuff.f[dirN   ])[k];
      (D.f[dirN   ])[kn   ] = (Dbuff.f[dirS   ])[k];
      (D.f[dirB   ])[kb   ] = (Dbuff.f[dirT   ])[k];
      (D.f[dirT   ])[kt   ] = (Dbuff.f[dirB   ])[k];
      (D.f[dirSW  ])[ksw  ] = (Dbuff.f[dirNE  ])[k];
      (D.f[dirNE  ])[kne  ] = (Dbuff.f[dirSW  ])[k];
      (D.f[dirNW  ])[knw  ] = (Dbuff.f[dirSE  ])[k];
      (D.f[dirSE  ])[kse  ] = (Dbuff.f[dirNW  ])[k];
      (D.f[dirBW  ])[kbw  ] = (Dbuff.f[dirTE  ])[k];
      (D.f[dirTE  ])[kte  ] = (Dbuff.f[dirBW  ])[k];
      (D.f[dirTW  ])[ktw  ] = (Dbuff.f[dirBE  ])[k];
      (D.f[dirBE  ])[kbe  ] = (Dbuff.f[dirTW  ])[k];
      (D.f[dirBS  ])[kbs  ] = (Dbuff.f[dirTN  ])[k];
      (D.f[dirTN  ])[ktn  ] = (Dbuff.f[dirBS  ])[k];
      (D.f[dirTS  ])[kts  ] = (Dbuff.f[dirBN  ])[k];
      (D.f[dirBN  ])[kbn  ] = (Dbuff.f[dirTS  ])[k];
      (D.f[dirZERO])[kzero] = (Dbuff.f[dirZERO])[k];
      (D.f[dirBSW ])[kbsw ] = (Dbuff.f[dirTNE ])[k];
      (D.f[dirBNE ])[kbne ] = (Dbuff.f[dirTSW ])[k];
      (D.f[dirBNW ])[kbnw ] = (Dbuff.f[dirTSE ])[k];
      (D.f[dirBSE ])[kbse ] = (Dbuff.f[dirTNW ])[k];
      (D.f[dirTSW ])[ktsw ] = (Dbuff.f[dirBNE ])[k];
      (D.f[dirTNE ])[ktne ] = (Dbuff.f[dirBSW ])[k];
      (D.f[dirTNW ])[ktnw ] = (Dbuff.f[dirBSE ])[k];
      (D.f[dirTSE ])[ktse ] = (Dbuff.f[dirBNW ])[k];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void getSendFsPre27(doubflo* DD,
										  doubflo* bufferFs,
										  int* sendIndex,
                                          int buffmax,
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

   if(k<buffmax)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //set index
      unsigned int kIndex  = sendIndex[k];
      unsigned int kzero   = kIndex;
      unsigned int ke      = kIndex;
      unsigned int kw      = neighborX[kIndex];
      unsigned int kn      = kIndex;
      unsigned int ks      = neighborY[kIndex];
      unsigned int kt      = kIndex;
      unsigned int kb      = neighborZ[kIndex];
      unsigned int ksw     = neighborY[kw];
      unsigned int kne     = kIndex;
      unsigned int kse     = ks;
      unsigned int knw     = kw;
      unsigned int kbw     = neighborZ[kw];
      unsigned int kte     = kIndex;
      unsigned int kbe     = kb;
      unsigned int ktw     = kw;
      unsigned int kbs     = neighborZ[ks];
      unsigned int ktn     = kIndex;
      unsigned int kbn     = kb;
      unsigned int kts     = ks;
      unsigned int ktse    = ks;
      unsigned int kbnw    = kbw;
      unsigned int ktnw    = kw;
      unsigned int kbse    = kbs;
      unsigned int ktsw    = ksw;
      unsigned int kbne    = kb;
      unsigned int ktne    = kIndex;
      unsigned int kbsw    = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Fs
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dirE   ] = &bufferFs[dirE   *buffmax];
      Dbuff.f[dirW   ] = &bufferFs[dirW   *buffmax];
      Dbuff.f[dirN   ] = &bufferFs[dirN   *buffmax];
      Dbuff.f[dirS   ] = &bufferFs[dirS   *buffmax];
      Dbuff.f[dirT   ] = &bufferFs[dirT   *buffmax];
      Dbuff.f[dirB   ] = &bufferFs[dirB   *buffmax];
      Dbuff.f[dirNE  ] = &bufferFs[dirNE  *buffmax];
      Dbuff.f[dirSW  ] = &bufferFs[dirSW  *buffmax];
      Dbuff.f[dirSE  ] = &bufferFs[dirSE  *buffmax];
      Dbuff.f[dirNW  ] = &bufferFs[dirNW  *buffmax];
      Dbuff.f[dirTE  ] = &bufferFs[dirTE  *buffmax];
      Dbuff.f[dirBW  ] = &bufferFs[dirBW  *buffmax];
      Dbuff.f[dirBE  ] = &bufferFs[dirBE  *buffmax];
      Dbuff.f[dirTW  ] = &bufferFs[dirTW  *buffmax];
      Dbuff.f[dirTN  ] = &bufferFs[dirTN  *buffmax];
      Dbuff.f[dirBS  ] = &bufferFs[dirBS  *buffmax];
      Dbuff.f[dirBN  ] = &bufferFs[dirBN  *buffmax];
      Dbuff.f[dirTS  ] = &bufferFs[dirTS  *buffmax];
      Dbuff.f[dirZERO] = &bufferFs[dirZERO*buffmax];
      Dbuff.f[dirTNE ] = &bufferFs[dirTNE *buffmax];
      Dbuff.f[dirTSW ] = &bufferFs[dirTSW *buffmax];
      Dbuff.f[dirTSE ] = &bufferFs[dirTSE *buffmax];
      Dbuff.f[dirTNW ] = &bufferFs[dirTNW *buffmax];
      Dbuff.f[dirBNE ] = &bufferFs[dirBNE *buffmax];
      Dbuff.f[dirBSW ] = &bufferFs[dirBSW *buffmax];
      Dbuff.f[dirBSE ] = &bufferFs[dirBSE *buffmax];
      Dbuff.f[dirBNW ] = &bufferFs[dirBNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      (Dbuff.f[dirE   ])[k] = (D.f[dirE   ])[ke   ];
      (Dbuff.f[dirW   ])[k] = (D.f[dirW   ])[kw   ];
      (Dbuff.f[dirN   ])[k] = (D.f[dirN   ])[kn   ];
      (Dbuff.f[dirS   ])[k] = (D.f[dirS   ])[ks   ];
      (Dbuff.f[dirT   ])[k] = (D.f[dirT   ])[kt   ];
      (Dbuff.f[dirB   ])[k] = (D.f[dirB   ])[kb   ];
      (Dbuff.f[dirNE  ])[k] = (D.f[dirNE  ])[kne  ];
      (Dbuff.f[dirSW  ])[k] = (D.f[dirSW  ])[ksw  ];
      (Dbuff.f[dirSE  ])[k] = (D.f[dirSE  ])[kse  ];
      (Dbuff.f[dirNW  ])[k] = (D.f[dirNW  ])[knw  ];
      (Dbuff.f[dirTE  ])[k] = (D.f[dirTE  ])[kte  ];
      (Dbuff.f[dirBW  ])[k] = (D.f[dirBW  ])[kbw  ];
      (Dbuff.f[dirBE  ])[k] = (D.f[dirBE  ])[kbe  ];
      (Dbuff.f[dirTW  ])[k] = (D.f[dirTW  ])[ktw  ];
      (Dbuff.f[dirTN  ])[k] = (D.f[dirTN  ])[ktn  ];
      (Dbuff.f[dirBS  ])[k] = (D.f[dirBS  ])[kbs  ];
      (Dbuff.f[dirBN  ])[k] = (D.f[dirBN  ])[kbn  ];
      (Dbuff.f[dirTS  ])[k] = (D.f[dirTS  ])[kts  ];
      (Dbuff.f[dirZERO])[k] = (D.f[dirZERO])[kzero];
      (Dbuff.f[dirTNE ])[k] = (D.f[dirTNE ])[ktne ];
      (Dbuff.f[dirTSW ])[k] = (D.f[dirTSW ])[ktsw ];
      (Dbuff.f[dirTSE ])[k] = (D.f[dirTSE ])[ktse ];
      (Dbuff.f[dirTNW ])[k] = (D.f[dirTNW ])[ktnw ];
      (Dbuff.f[dirBNE ])[k] = (D.f[dirBNE ])[kbne ];
      (Dbuff.f[dirBSW ])[k] = (D.f[dirBSW ])[kbsw ];
      (Dbuff.f[dirBSE ])[k] = (D.f[dirBSE ])[kbse ];
      (Dbuff.f[dirBNW ])[k] = (D.f[dirBNW ])[kbnw ];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void setRecvFsPre27(doubflo* DD,
										  doubflo* bufferFs,
										  int* recvIndex,
                                          int buffmax,
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

   if(k<buffmax)
   {
      ////////////////////////////////////////////////////////////////////////////////
      //set index
      unsigned int kIndex  = recvIndex[k];
      unsigned int kzero   = kIndex;
      unsigned int ke      = kIndex;
      unsigned int kw      = neighborX[kIndex];
      unsigned int kn      = kIndex;
      unsigned int ks      = neighborY[kIndex];
      unsigned int kt      = kIndex;
      unsigned int kb      = neighborZ[kIndex];
      unsigned int ksw     = neighborY[kw];
      unsigned int kne     = kIndex;
      unsigned int kse     = ks;
      unsigned int knw     = kw;
      unsigned int kbw     = neighborZ[kw];
      unsigned int kte     = kIndex;
      unsigned int kbe     = kb;
      unsigned int ktw     = kw;
      unsigned int kbs     = neighborZ[ks];
      unsigned int ktn     = kIndex;
      unsigned int kbn     = kb;
      unsigned int kts     = ks;
      unsigned int ktse    = ks;
      unsigned int kbnw    = kbw;
      unsigned int ktnw    = kw;
      unsigned int kbse    = kbs;
      unsigned int ktsw    = ksw;
      unsigned int kbne    = kb;
      unsigned int ktne    = kIndex;
      unsigned int kbsw    = neighborZ[ksw];
      ////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Fs
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dirE   ] = &bufferFs[dirE   *buffmax];
      Dbuff.f[dirW   ] = &bufferFs[dirW   *buffmax];
      Dbuff.f[dirN   ] = &bufferFs[dirN   *buffmax];
      Dbuff.f[dirS   ] = &bufferFs[dirS   *buffmax];
      Dbuff.f[dirT   ] = &bufferFs[dirT   *buffmax];
      Dbuff.f[dirB   ] = &bufferFs[dirB   *buffmax];
      Dbuff.f[dirNE  ] = &bufferFs[dirNE  *buffmax];
      Dbuff.f[dirSW  ] = &bufferFs[dirSW  *buffmax];
      Dbuff.f[dirSE  ] = &bufferFs[dirSE  *buffmax];
      Dbuff.f[dirNW  ] = &bufferFs[dirNW  *buffmax];
      Dbuff.f[dirTE  ] = &bufferFs[dirTE  *buffmax];
      Dbuff.f[dirBW  ] = &bufferFs[dirBW  *buffmax];
      Dbuff.f[dirBE  ] = &bufferFs[dirBE  *buffmax];
      Dbuff.f[dirTW  ] = &bufferFs[dirTW  *buffmax];
      Dbuff.f[dirTN  ] = &bufferFs[dirTN  *buffmax];
      Dbuff.f[dirBS  ] = &bufferFs[dirBS  *buffmax];
      Dbuff.f[dirBN  ] = &bufferFs[dirBN  *buffmax];
      Dbuff.f[dirTS  ] = &bufferFs[dirTS  *buffmax];
      Dbuff.f[dirZERO] = &bufferFs[dirZERO*buffmax];
      Dbuff.f[dirTNE ] = &bufferFs[dirTNE *buffmax];
      Dbuff.f[dirTSW ] = &bufferFs[dirTSW *buffmax];
      Dbuff.f[dirTSE ] = &bufferFs[dirTSE *buffmax];
      Dbuff.f[dirTNW ] = &bufferFs[dirTNW *buffmax];
      Dbuff.f[dirBNE ] = &bufferFs[dirBNE *buffmax];
      Dbuff.f[dirBSW ] = &bufferFs[dirBSW *buffmax];
      Dbuff.f[dirBSE ] = &bufferFs[dirBSE *buffmax];
      Dbuff.f[dirBNW ] = &bufferFs[dirBNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      (D.f[dirE   ])[ke   ] = (Dbuff.f[dirE   ])[k];
      (D.f[dirW   ])[kw   ] = (Dbuff.f[dirW   ])[k];
      (D.f[dirN   ])[kn   ] = (Dbuff.f[dirN   ])[k];
      (D.f[dirS   ])[ks   ] = (Dbuff.f[dirS   ])[k];
      (D.f[dirT   ])[kt   ] = (Dbuff.f[dirT   ])[k];
      (D.f[dirB   ])[kb   ] = (Dbuff.f[dirB   ])[k];
      (D.f[dirNE  ])[kne  ] = (Dbuff.f[dirNE  ])[k];
      (D.f[dirSW  ])[ksw  ] = (Dbuff.f[dirSW  ])[k];
      (D.f[dirSE  ])[kse  ] = (Dbuff.f[dirSE  ])[k];
      (D.f[dirNW  ])[knw  ] = (Dbuff.f[dirNW  ])[k];
      (D.f[dirTE  ])[kte  ] = (Dbuff.f[dirTE  ])[k];
      (D.f[dirBW  ])[kbw  ] = (Dbuff.f[dirBW  ])[k];
      (D.f[dirBE  ])[kbe  ] = (Dbuff.f[dirBE  ])[k];
      (D.f[dirTW  ])[ktw  ] = (Dbuff.f[dirTW  ])[k];
      (D.f[dirTN  ])[ktn  ] = (Dbuff.f[dirTN  ])[k];
      (D.f[dirBS  ])[kbs  ] = (Dbuff.f[dirBS  ])[k];
      (D.f[dirBN  ])[kbn  ] = (Dbuff.f[dirBN  ])[k];
      (D.f[dirTS  ])[kts  ] = (Dbuff.f[dirTS  ])[k];
      (D.f[dirZERO])[kzero] = (Dbuff.f[dirZERO])[k];
      (D.f[dirTNE ])[ktne ] = (Dbuff.f[dirTNE ])[k];
      (D.f[dirTSW ])[ktsw ] = (Dbuff.f[dirTSW ])[k];
      (D.f[dirTSE ])[ktse ] = (Dbuff.f[dirTSE ])[k];
      (D.f[dirTNW ])[ktnw ] = (Dbuff.f[dirTNW ])[k];
      (D.f[dirBNE ])[kbne ] = (Dbuff.f[dirBNE ])[k];
      (D.f[dirBSW ])[kbsw ] = (Dbuff.f[dirBSW ])[k];
      (D.f[dirBSE ])[kbse ] = (Dbuff.f[dirBSE ])[k];
      (D.f[dirBNW ])[kbnw ] = (Dbuff.f[dirBNW ])[k];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


