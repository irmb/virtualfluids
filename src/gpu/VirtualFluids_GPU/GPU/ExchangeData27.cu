/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void getSendFsPost27(real* DD,
										   real* bufferFs,
										   int* sendIndex,
                                           int buffmax,
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[E   ] = &bufferFs[E   *buffmax];
      Dbuff.f[W   ] = &bufferFs[W   *buffmax];
      Dbuff.f[N   ] = &bufferFs[N   *buffmax];
      Dbuff.f[S   ] = &bufferFs[S   *buffmax];
      Dbuff.f[T   ] = &bufferFs[T   *buffmax];
      Dbuff.f[B   ] = &bufferFs[B   *buffmax];
      Dbuff.f[NE  ] = &bufferFs[NE  *buffmax];
      Dbuff.f[SW  ] = &bufferFs[SW  *buffmax];
      Dbuff.f[SE  ] = &bufferFs[SE  *buffmax];
      Dbuff.f[NW  ] = &bufferFs[NW  *buffmax];
      Dbuff.f[TE  ] = &bufferFs[TE  *buffmax];
      Dbuff.f[BW  ] = &bufferFs[BW  *buffmax];
      Dbuff.f[BE  ] = &bufferFs[BE  *buffmax];
      Dbuff.f[TW  ] = &bufferFs[TW  *buffmax];
      Dbuff.f[TN  ] = &bufferFs[TN  *buffmax];
      Dbuff.f[BS  ] = &bufferFs[BS  *buffmax];
      Dbuff.f[BN  ] = &bufferFs[BN  *buffmax];
      Dbuff.f[TS  ] = &bufferFs[TS  *buffmax];
      Dbuff.f[REST] = &bufferFs[REST*buffmax];
      Dbuff.f[TNE ] = &bufferFs[TNE *buffmax];
      Dbuff.f[TSW ] = &bufferFs[TSW *buffmax];
      Dbuff.f[TSE ] = &bufferFs[TSE *buffmax];
      Dbuff.f[TNW ] = &bufferFs[TNW *buffmax];
      Dbuff.f[BNE ] = &bufferFs[BNE *buffmax];
      Dbuff.f[BSW ] = &bufferFs[BSW *buffmax];
      Dbuff.f[BSE ] = &bufferFs[BSE *buffmax];
      Dbuff.f[BNW ] = &bufferFs[BNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      //(Dbuff.f[E   ])[k] = (D.f[E   ])[ke   ];
      //(Dbuff.f[W   ])[k] = (D.f[W   ])[kw   ];
      //(Dbuff.f[N   ])[k] = (D.f[N   ])[kn   ];
      //(Dbuff.f[S   ])[k] = (D.f[S   ])[ks   ];
      //(Dbuff.f[T   ])[k] = (D.f[T   ])[kt   ];
      //(Dbuff.f[B   ])[k] = (D.f[B   ])[kb   ];
      //(Dbuff.f[NE  ])[k] = (D.f[NE  ])[kne  ];
      //(Dbuff.f[SW  ])[k] = (D.f[SW  ])[ksw  ];
      //(Dbuff.f[SE  ])[k] = (D.f[SE  ])[kse  ];
      //(Dbuff.f[NW  ])[k] = (D.f[NW  ])[knw  ];
      //(Dbuff.f[TE  ])[k] = (D.f[TE  ])[kte  ];
      //(Dbuff.f[BW  ])[k] = (D.f[BW  ])[kbw  ];
      //(Dbuff.f[BE  ])[k] = (D.f[BE  ])[kbe  ];
      //(Dbuff.f[TW  ])[k] = (D.f[TW  ])[ktw  ];
      //(Dbuff.f[TN  ])[k] = (D.f[TN  ])[ktn  ];
      //(Dbuff.f[BS  ])[k] = (D.f[BS  ])[kbs  ];
      //(Dbuff.f[BN  ])[k] = (D.f[BN  ])[kbn  ];
      //(Dbuff.f[TS  ])[k] = (D.f[TS  ])[kts  ];
      //(Dbuff.f[REST])[k] = (D.f[REST])[kzero];
      //(Dbuff.f[TNE ])[k] = (D.f[TNE ])[ktne ];
      //(Dbuff.f[TSW ])[k] = (D.f[TSW ])[ktsw ];
      //(Dbuff.f[TSE ])[k] = (D.f[TSE ])[ktse ];
      //(Dbuff.f[TNW ])[k] = (D.f[TNW ])[ktnw ];
      //(Dbuff.f[BNE ])[k] = (D.f[BNE ])[kbne ];
      //(Dbuff.f[BSW ])[k] = (D.f[BSW ])[kbsw ];
      //(Dbuff.f[BSE ])[k] = (D.f[BSE ])[kbse ];
      //(Dbuff.f[BNW ])[k] = (D.f[BNW ])[kbnw ];
      (Dbuff.f[E   ])[k] = (D.f[W   ])[kw   ];
      (Dbuff.f[W   ])[k] = (D.f[E   ])[ke   ];
      (Dbuff.f[N   ])[k] = (D.f[S   ])[ks   ];
      (Dbuff.f[S   ])[k] = (D.f[N   ])[kn   ];
      (Dbuff.f[T   ])[k] = (D.f[B   ])[kb   ];
      (Dbuff.f[B   ])[k] = (D.f[T   ])[kt   ];
      (Dbuff.f[NE  ])[k] = (D.f[SW  ])[ksw  ];
      (Dbuff.f[SW  ])[k] = (D.f[NE  ])[kne  ];
      (Dbuff.f[SE  ])[k] = (D.f[NW  ])[knw  ];
      (Dbuff.f[NW  ])[k] = (D.f[SE  ])[kse  ];
      (Dbuff.f[TE  ])[k] = (D.f[BW  ])[kbw  ];
      (Dbuff.f[BW  ])[k] = (D.f[TE  ])[kte  ];
      (Dbuff.f[BE  ])[k] = (D.f[TW  ])[ktw  ];
      (Dbuff.f[TW  ])[k] = (D.f[BE  ])[kbe  ];
      (Dbuff.f[TN  ])[k] = (D.f[BS  ])[kbs  ];
      (Dbuff.f[BS  ])[k] = (D.f[TN  ])[ktn  ];
      (Dbuff.f[BN  ])[k] = (D.f[TS  ])[kts  ];
      (Dbuff.f[TS  ])[k] = (D.f[BN  ])[kbn  ];
      (Dbuff.f[REST])[k] = (D.f[REST])[kzero];
      (Dbuff.f[TNE ])[k] = (D.f[BSW ])[kbsw ];
      (Dbuff.f[TSW ])[k] = (D.f[BNE ])[kbne ];
      (Dbuff.f[TSE ])[k] = (D.f[BNW ])[kbnw ];
      (Dbuff.f[TNW ])[k] = (D.f[BSE ])[kbse ];
      (Dbuff.f[BNE ])[k] = (D.f[TSW ])[ktsw ];
      (Dbuff.f[BSW ])[k] = (D.f[TNE ])[ktne ];
      (Dbuff.f[BSE ])[k] = (D.f[TNW ])[ktnw ];
      (Dbuff.f[BNW ])[k] = (D.f[TSE ])[ktse ];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void setRecvFsPost27(real* DD,
										   real* bufferFs,
										   int* recvIndex,
                                           int buffmax,
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[E   ] = &bufferFs[E   *buffmax];
      Dbuff.f[W   ] = &bufferFs[W   *buffmax];
      Dbuff.f[N   ] = &bufferFs[N   *buffmax];
      Dbuff.f[S   ] = &bufferFs[S   *buffmax];
      Dbuff.f[T   ] = &bufferFs[T   *buffmax];
      Dbuff.f[B   ] = &bufferFs[B   *buffmax];
      Dbuff.f[NE  ] = &bufferFs[NE  *buffmax];
      Dbuff.f[SW  ] = &bufferFs[SW  *buffmax];
      Dbuff.f[SE  ] = &bufferFs[SE  *buffmax];
      Dbuff.f[NW  ] = &bufferFs[NW  *buffmax];
      Dbuff.f[TE  ] = &bufferFs[TE  *buffmax];
      Dbuff.f[BW  ] = &bufferFs[BW  *buffmax];
      Dbuff.f[BE  ] = &bufferFs[BE  *buffmax];
      Dbuff.f[TW  ] = &bufferFs[TW  *buffmax];
      Dbuff.f[TN  ] = &bufferFs[TN  *buffmax];
      Dbuff.f[BS  ] = &bufferFs[BS  *buffmax];
      Dbuff.f[BN  ] = &bufferFs[BN  *buffmax];
      Dbuff.f[TS  ] = &bufferFs[TS  *buffmax];
      Dbuff.f[REST] = &bufferFs[REST*buffmax];
      Dbuff.f[TNE ] = &bufferFs[TNE *buffmax];
      Dbuff.f[TSW ] = &bufferFs[TSW *buffmax];
      Dbuff.f[TSE ] = &bufferFs[TSE *buffmax];
      Dbuff.f[TNW ] = &bufferFs[TNW *buffmax];
      Dbuff.f[BNE ] = &bufferFs[BNE *buffmax];
      Dbuff.f[BSW ] = &bufferFs[BSW *buffmax];
      Dbuff.f[BSE ] = &bufferFs[BSE *buffmax];
      Dbuff.f[BNW ] = &bufferFs[BNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      //(D.f[E   ])[ke   ] = (Dbuff.f[E   ])[k];
      //(D.f[W   ])[kw   ] = (Dbuff.f[W   ])[k];
      //(D.f[N   ])[kn   ] = (Dbuff.f[N   ])[k];
      //(D.f[S   ])[ks   ] = (Dbuff.f[S   ])[k];
      //(D.f[T   ])[kt   ] = (Dbuff.f[T   ])[k];
      //(D.f[B   ])[kb   ] = (Dbuff.f[B   ])[k];
      //(D.f[NE  ])[kne  ] = (Dbuff.f[NE  ])[k];
      //(D.f[SW  ])[ksw  ] = (Dbuff.f[SW  ])[k];
      //(D.f[SE  ])[kse  ] = (Dbuff.f[SE  ])[k];
      //(D.f[NW  ])[knw  ] = (Dbuff.f[NW  ])[k];
      //(D.f[TE  ])[kte  ] = (Dbuff.f[TE  ])[k];
      //(D.f[BW  ])[kbw  ] = (Dbuff.f[BW  ])[k];
      //(D.f[BE  ])[kbe  ] = (Dbuff.f[BE  ])[k];
      //(D.f[TW  ])[ktw  ] = (Dbuff.f[TW  ])[k];
      //(D.f[TN  ])[ktn  ] = (Dbuff.f[TN  ])[k];
      //(D.f[BS  ])[kbs  ] = (Dbuff.f[BS  ])[k];
      //(D.f[BN  ])[kbn  ] = (Dbuff.f[BN  ])[k];
      //(D.f[TS  ])[kts  ] = (Dbuff.f[TS  ])[k];
      //(D.f[REST])[kzero] = (Dbuff.f[REST])[k];
      //(D.f[TNE ])[ktne ] = (Dbuff.f[TNE ])[k];
      //(D.f[TSW ])[ktsw ] = (Dbuff.f[TSW ])[k];
      //(D.f[TSE ])[ktse ] = (Dbuff.f[TSE ])[k];
      //(D.f[TNW ])[ktnw ] = (Dbuff.f[TNW ])[k];
      //(D.f[BNE ])[kbne ] = (Dbuff.f[BNE ])[k];
      //(D.f[BSW ])[kbsw ] = (Dbuff.f[BSW ])[k];
      //(D.f[BSE ])[kbse ] = (Dbuff.f[BSE ])[k];
      //(D.f[BNW ])[kbnw ] = (Dbuff.f[BNW ])[k];
      (D.f[W   ])[kw   ] = (Dbuff.f[E   ])[k];
      (D.f[E   ])[ke   ] = (Dbuff.f[W   ])[k];
      (D.f[S   ])[ks   ] = (Dbuff.f[N   ])[k];
      (D.f[N   ])[kn   ] = (Dbuff.f[S   ])[k];
      (D.f[B   ])[kb   ] = (Dbuff.f[T   ])[k];
      (D.f[T   ])[kt   ] = (Dbuff.f[B   ])[k];
      (D.f[SW  ])[ksw  ] = (Dbuff.f[NE  ])[k];
      (D.f[NE  ])[kne  ] = (Dbuff.f[SW  ])[k];
      (D.f[NW  ])[knw  ] = (Dbuff.f[SE  ])[k];
      (D.f[SE  ])[kse  ] = (Dbuff.f[NW  ])[k];
      (D.f[BW  ])[kbw  ] = (Dbuff.f[TE  ])[k];
      (D.f[TE  ])[kte  ] = (Dbuff.f[BW  ])[k];
      (D.f[TW  ])[ktw  ] = (Dbuff.f[BE  ])[k];
      (D.f[BE  ])[kbe  ] = (Dbuff.f[TW  ])[k];
      (D.f[BS  ])[kbs  ] = (Dbuff.f[TN  ])[k];
      (D.f[TN  ])[ktn  ] = (Dbuff.f[BS  ])[k];
      (D.f[TS  ])[kts  ] = (Dbuff.f[BN  ])[k];
      (D.f[BN  ])[kbn  ] = (Dbuff.f[TS  ])[k];
      (D.f[REST])[kzero] = (Dbuff.f[REST])[k];
      (D.f[BSW ])[kbsw ] = (Dbuff.f[TNE ])[k];
      (D.f[BNE ])[kbne ] = (Dbuff.f[TSW ])[k];
      (D.f[BNW ])[kbnw ] = (Dbuff.f[TSE ])[k];
      (D.f[BSE ])[kbse ] = (Dbuff.f[TNW ])[k];
      (D.f[TSW ])[ktsw ] = (Dbuff.f[BNE ])[k];
      (D.f[TNE ])[ktne ] = (Dbuff.f[BSW ])[k];
      (D.f[TNW ])[ktnw ] = (Dbuff.f[BSE ])[k];
      (D.f[TSE ])[ktse ] = (Dbuff.f[BNW ])[k];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void getSendFsPre27(real* DD,
										  real* bufferFs,
										  int* sendIndex,
                                          int buffmax,
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[E   ] = &bufferFs[E   *buffmax];
      Dbuff.f[W   ] = &bufferFs[W   *buffmax];
      Dbuff.f[N   ] = &bufferFs[N   *buffmax];
      Dbuff.f[S   ] = &bufferFs[S   *buffmax];
      Dbuff.f[T   ] = &bufferFs[T   *buffmax];
      Dbuff.f[B   ] = &bufferFs[B   *buffmax];
      Dbuff.f[NE  ] = &bufferFs[NE  *buffmax];
      Dbuff.f[SW  ] = &bufferFs[SW  *buffmax];
      Dbuff.f[SE  ] = &bufferFs[SE  *buffmax];
      Dbuff.f[NW  ] = &bufferFs[NW  *buffmax];
      Dbuff.f[TE  ] = &bufferFs[TE  *buffmax];
      Dbuff.f[BW  ] = &bufferFs[BW  *buffmax];
      Dbuff.f[BE  ] = &bufferFs[BE  *buffmax];
      Dbuff.f[TW  ] = &bufferFs[TW  *buffmax];
      Dbuff.f[TN  ] = &bufferFs[TN  *buffmax];
      Dbuff.f[BS  ] = &bufferFs[BS  *buffmax];
      Dbuff.f[BN  ] = &bufferFs[BN  *buffmax];
      Dbuff.f[TS  ] = &bufferFs[TS  *buffmax];
      Dbuff.f[REST] = &bufferFs[REST*buffmax];
      Dbuff.f[TNE ] = &bufferFs[TNE *buffmax];
      Dbuff.f[TSW ] = &bufferFs[TSW *buffmax];
      Dbuff.f[TSE ] = &bufferFs[TSE *buffmax];
      Dbuff.f[TNW ] = &bufferFs[TNW *buffmax];
      Dbuff.f[BNE ] = &bufferFs[BNE *buffmax];
      Dbuff.f[BSW ] = &bufferFs[BSW *buffmax];
      Dbuff.f[BSE ] = &bufferFs[BSE *buffmax];
      Dbuff.f[BNW ] = &bufferFs[BNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      (Dbuff.f[E   ])[k] = (D.f[E   ])[ke   ];
      (Dbuff.f[W   ])[k] = (D.f[W   ])[kw   ];
      (Dbuff.f[N   ])[k] = (D.f[N   ])[kn   ];
      (Dbuff.f[S   ])[k] = (D.f[S   ])[ks   ];
      (Dbuff.f[T   ])[k] = (D.f[T   ])[kt   ];
      (Dbuff.f[B   ])[k] = (D.f[B   ])[kb   ];
      (Dbuff.f[NE  ])[k] = (D.f[NE  ])[kne  ];
      (Dbuff.f[SW  ])[k] = (D.f[SW  ])[ksw  ];
      (Dbuff.f[SE  ])[k] = (D.f[SE  ])[kse  ];
      (Dbuff.f[NW  ])[k] = (D.f[NW  ])[knw  ];
      (Dbuff.f[TE  ])[k] = (D.f[TE  ])[kte  ];
      (Dbuff.f[BW  ])[k] = (D.f[BW  ])[kbw  ];
      (Dbuff.f[BE  ])[k] = (D.f[BE  ])[kbe  ];
      (Dbuff.f[TW  ])[k] = (D.f[TW  ])[ktw  ];
      (Dbuff.f[TN  ])[k] = (D.f[TN  ])[ktn  ];
      (Dbuff.f[BS  ])[k] = (D.f[BS  ])[kbs  ];
      (Dbuff.f[BN  ])[k] = (D.f[BN  ])[kbn  ];
      (Dbuff.f[TS  ])[k] = (D.f[TS  ])[kts  ];
      (Dbuff.f[REST])[k] = (D.f[REST])[kzero];
      (Dbuff.f[TNE ])[k] = (D.f[TNE ])[ktne ];
      (Dbuff.f[TSW ])[k] = (D.f[TSW ])[ktsw ];
      (Dbuff.f[TSE ])[k] = (D.f[TSE ])[ktse ];
      (Dbuff.f[TNW ])[k] = (D.f[TNW ])[ktnw ];
      (Dbuff.f[BNE ])[k] = (D.f[BNE ])[kbne ];
      (Dbuff.f[BSW ])[k] = (D.f[BSW ])[kbsw ];
      (Dbuff.f[BSE ])[k] = (D.f[BSE ])[kbse ];
      (Dbuff.f[BNW ])[k] = (D.f[BNW ])[kbnw ];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void setRecvFsPre27(real* DD,
										  real* bufferFs,
										  int* recvIndex,
                                          int buffmax,
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
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[E   ] = &bufferFs[E   *buffmax];
      Dbuff.f[W   ] = &bufferFs[W   *buffmax];
      Dbuff.f[N   ] = &bufferFs[N   *buffmax];
      Dbuff.f[S   ] = &bufferFs[S   *buffmax];
      Dbuff.f[T   ] = &bufferFs[T   *buffmax];
      Dbuff.f[B   ] = &bufferFs[B   *buffmax];
      Dbuff.f[NE  ] = &bufferFs[NE  *buffmax];
      Dbuff.f[SW  ] = &bufferFs[SW  *buffmax];
      Dbuff.f[SE  ] = &bufferFs[SE  *buffmax];
      Dbuff.f[NW  ] = &bufferFs[NW  *buffmax];
      Dbuff.f[TE  ] = &bufferFs[TE  *buffmax];
      Dbuff.f[BW  ] = &bufferFs[BW  *buffmax];
      Dbuff.f[BE  ] = &bufferFs[BE  *buffmax];
      Dbuff.f[TW  ] = &bufferFs[TW  *buffmax];
      Dbuff.f[TN  ] = &bufferFs[TN  *buffmax];
      Dbuff.f[BS  ] = &bufferFs[BS  *buffmax];
      Dbuff.f[BN  ] = &bufferFs[BN  *buffmax];
      Dbuff.f[TS  ] = &bufferFs[TS  *buffmax];
      Dbuff.f[REST] = &bufferFs[REST*buffmax];
      Dbuff.f[TNE ] = &bufferFs[TNE *buffmax];
      Dbuff.f[TSW ] = &bufferFs[TSW *buffmax];
      Dbuff.f[TSE ] = &bufferFs[TSE *buffmax];
      Dbuff.f[TNW ] = &bufferFs[TNW *buffmax];
      Dbuff.f[BNE ] = &bufferFs[BNE *buffmax];
      Dbuff.f[BSW ] = &bufferFs[BSW *buffmax];
      Dbuff.f[BSE ] = &bufferFs[BSE *buffmax];
      Dbuff.f[BNW ] = &bufferFs[BNW *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      (D.f[E   ])[ke   ] = (Dbuff.f[E   ])[k];
      (D.f[W   ])[kw   ] = (Dbuff.f[W   ])[k];
      (D.f[N   ])[kn   ] = (Dbuff.f[N   ])[k];
      (D.f[S   ])[ks   ] = (Dbuff.f[S   ])[k];
      (D.f[T   ])[kt   ] = (Dbuff.f[T   ])[k];
      (D.f[B   ])[kb   ] = (Dbuff.f[B   ])[k];
      (D.f[NE  ])[kne  ] = (Dbuff.f[NE  ])[k];
      (D.f[SW  ])[ksw  ] = (Dbuff.f[SW  ])[k];
      (D.f[SE  ])[kse  ] = (Dbuff.f[SE  ])[k];
      (D.f[NW  ])[knw  ] = (Dbuff.f[NW  ])[k];
      (D.f[TE  ])[kte  ] = (Dbuff.f[TE  ])[k];
      (D.f[BW  ])[kbw  ] = (Dbuff.f[BW  ])[k];
      (D.f[BE  ])[kbe  ] = (Dbuff.f[BE  ])[k];
      (D.f[TW  ])[ktw  ] = (Dbuff.f[TW  ])[k];
      (D.f[TN  ])[ktn  ] = (Dbuff.f[TN  ])[k];
      (D.f[BS  ])[kbs  ] = (Dbuff.f[BS  ])[k];
      (D.f[BN  ])[kbn  ] = (Dbuff.f[BN  ])[k];
      (D.f[TS  ])[kts  ] = (Dbuff.f[TS  ])[k];
      (D.f[REST])[kzero] = (Dbuff.f[REST])[k];
      (D.f[TNE ])[ktne ] = (Dbuff.f[TNE ])[k];
      (D.f[TSW ])[ktsw ] = (Dbuff.f[TSW ])[k];
      (D.f[TSE ])[ktse ] = (Dbuff.f[TSE ])[k];
      (D.f[TNW ])[ktnw ] = (Dbuff.f[TNW ])[k];
      (D.f[BNE ])[kbne ] = (Dbuff.f[BNE ])[k];
      (D.f[BSW ])[kbsw ] = (Dbuff.f[BSW ])[k];
      (D.f[BSE ])[kbse ] = (Dbuff.f[BSE ])[k];
      (D.f[BNW ])[kbnw ] = (Dbuff.f[BNW ])[k];
   }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void getSendGsF3(
	real* G6,
	real* bufferGs,
	int* sendIndex,
	int buffmax,
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

	if (k < buffmax)
	{
		////////////////////////////////////////////////////////////////////////////////
		//set index
		unsigned int kIndex = sendIndex[k];
		unsigned int kr = kIndex;
		unsigned int kw = neighborX[kIndex];
		unsigned int ks = neighborY[kIndex];
		unsigned int kb = neighborZ[kIndex];
		////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Gs
		Distributions6 G;
		if (isEvenTimestep)
		{
			G.g[E] = &G6[E   *size_Mat];
			G.g[W] = &G6[W   *size_Mat];
			G.g[N] = &G6[N   *size_Mat];
			G.g[S] = &G6[S   *size_Mat];
			G.g[T] = &G6[T   *size_Mat];
			G.g[B] = &G6[B   *size_Mat];
		}
		else
		{
			G.g[W] = &G6[E   *size_Mat];
			G.g[E] = &G6[W   *size_Mat];
			G.g[S] = &G6[N   *size_Mat];
			G.g[N] = &G6[S   *size_Mat];
			G.g[B] = &G6[T   *size_Mat];
			G.g[T] = &G6[B   *size_Mat];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Buffer Gs
		Distributions6 Dbuff;
		Dbuff.g[E] = &bufferGs[E   *buffmax];
		Dbuff.g[W] = &bufferGs[W   *buffmax];
		Dbuff.g[N] = &bufferGs[N   *buffmax];
		Dbuff.g[S] = &bufferGs[S   *buffmax];
		Dbuff.g[T] = &bufferGs[T   *buffmax];
		Dbuff.g[B] = &bufferGs[B   *buffmax];
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//write Gs to buffer
		(Dbuff.g[E])[k] = (G.g[W])[kw];
		(Dbuff.g[W])[k] = (G.g[E])[kr];
		(Dbuff.g[N])[k] = (G.g[S])[ks];
		(Dbuff.g[S])[k] = (G.g[N])[kr];
		(Dbuff.g[T])[k] = (G.g[B])[kb];
		(Dbuff.g[B])[k] = (G.g[T])[kr];
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void setRecvGsF3(
	real* G6,
	real* bufferGs,
	int* recvIndex,
	int buffmax,
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

	if (k < buffmax)
	{
		////////////////////////////////////////////////////////////////////////////////
		//set index
		unsigned int kIndex = recvIndex[k];
		unsigned int kr = kIndex;
		unsigned int kw = neighborX[kIndex];
		unsigned int ks = neighborY[kIndex];
		unsigned int kb = neighborZ[kIndex];
		////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Gs
		Distributions6 G;
		if (isEvenTimestep)
		{
			G.g[E] = &G6[E   *size_Mat];
			G.g[W] = &G6[W   *size_Mat];
			G.g[N] = &G6[N   *size_Mat];
			G.g[S] = &G6[S   *size_Mat];
			G.g[T] = &G6[T   *size_Mat];
			G.g[B] = &G6[B   *size_Mat];
		}
		else
		{
			G.g[W] = &G6[E   *size_Mat];
			G.g[E] = &G6[W   *size_Mat];
			G.g[S] = &G6[N   *size_Mat];
			G.g[N] = &G6[S   *size_Mat];
			G.g[B] = &G6[T   *size_Mat];
			G.g[T] = &G6[B   *size_Mat];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Buffer Gs
		Distributions6 Dbuff;
		Dbuff.g[E] = &bufferGs[E   *buffmax];
		Dbuff.g[W] = &bufferGs[W   *buffmax];
		Dbuff.g[N] = &bufferGs[N   *buffmax];
		Dbuff.g[S] = &bufferGs[S   *buffmax];
		Dbuff.g[T] = &bufferGs[T   *buffmax];
		Dbuff.g[B] = &bufferGs[B   *buffmax];
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//write buffer to Gs
		(G.g[W])[kw] = (Dbuff.g[E])[k];
		(G.g[E])[kr] = (Dbuff.g[W])[k];
		(G.g[S])[ks] = (Dbuff.g[N])[k];
		(G.g[N])[kr] = (Dbuff.g[S])[k];
		(G.g[B])[kb] = (Dbuff.g[T])[k];
		(G.g[T])[kr] = (Dbuff.g[B])[k];
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
