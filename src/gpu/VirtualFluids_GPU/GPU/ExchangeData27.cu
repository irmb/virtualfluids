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
         D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
      } 
      else
      {
         D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[DIR_P00   ] = &bufferFs[DIR_P00   *buffmax];
      Dbuff.f[DIR_M00   ] = &bufferFs[DIR_M00   *buffmax];
      Dbuff.f[DIR_0P0   ] = &bufferFs[DIR_0P0   *buffmax];
      Dbuff.f[DIR_0M0   ] = &bufferFs[DIR_0M0   *buffmax];
      Dbuff.f[DIR_00P   ] = &bufferFs[DIR_00P   *buffmax];
      Dbuff.f[DIR_00M   ] = &bufferFs[DIR_00M   *buffmax];
      Dbuff.f[DIR_PP0  ] = &bufferFs[DIR_PP0  *buffmax];
      Dbuff.f[DIR_MM0  ] = &bufferFs[DIR_MM0  *buffmax];
      Dbuff.f[DIR_PM0  ] = &bufferFs[DIR_PM0  *buffmax];
      Dbuff.f[DIR_MP0  ] = &bufferFs[DIR_MP0  *buffmax];
      Dbuff.f[DIR_P0P  ] = &bufferFs[DIR_P0P  *buffmax];
      Dbuff.f[DIR_M0M  ] = &bufferFs[DIR_M0M  *buffmax];
      Dbuff.f[DIR_P0M  ] = &bufferFs[DIR_P0M  *buffmax];
      Dbuff.f[DIR_M0P  ] = &bufferFs[DIR_M0P  *buffmax];
      Dbuff.f[DIR_0PP  ] = &bufferFs[DIR_0PP  *buffmax];
      Dbuff.f[DIR_0MM  ] = &bufferFs[DIR_0MM  *buffmax];
      Dbuff.f[DIR_0PM  ] = &bufferFs[DIR_0PM  *buffmax];
      Dbuff.f[DIR_0MP  ] = &bufferFs[DIR_0MP  *buffmax];
      Dbuff.f[DIR_000] = &bufferFs[DIR_000*buffmax];
      Dbuff.f[DIR_PPP ] = &bufferFs[DIR_PPP *buffmax];
      Dbuff.f[DIR_MMP ] = &bufferFs[DIR_MMP *buffmax];
      Dbuff.f[DIR_PMP ] = &bufferFs[DIR_PMP *buffmax];
      Dbuff.f[DIR_MPP ] = &bufferFs[DIR_MPP *buffmax];
      Dbuff.f[DIR_PPM ] = &bufferFs[DIR_PPM *buffmax];
      Dbuff.f[DIR_MMM ] = &bufferFs[DIR_MMM *buffmax];
      Dbuff.f[DIR_PMM ] = &bufferFs[DIR_PMM *buffmax];
      Dbuff.f[DIR_MPM ] = &bufferFs[DIR_MPM *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      //(Dbuff.f[DIR_P00   ])[k] = (D.f[DIR_P00   ])[ke   ];
      //(Dbuff.f[DIR_M00   ])[k] = (D.f[DIR_M00   ])[kw   ];
      //(Dbuff.f[DIR_0P0   ])[k] = (D.f[DIR_0P0   ])[kn   ];
      //(Dbuff.f[DIR_0M0   ])[k] = (D.f[DIR_0M0   ])[ks   ];
      //(Dbuff.f[DIR_00P   ])[k] = (D.f[DIR_00P   ])[kt   ];
      //(Dbuff.f[DIR_00M   ])[k] = (D.f[DIR_00M   ])[kb   ];
      //(Dbuff.f[DIR_PP0  ])[k] = (D.f[DIR_PP0  ])[kne  ];
      //(Dbuff.f[DIR_MM0  ])[k] = (D.f[DIR_MM0  ])[ksw  ];
      //(Dbuff.f[DIR_PM0  ])[k] = (D.f[DIR_PM0  ])[kse  ];
      //(Dbuff.f[DIR_MP0  ])[k] = (D.f[DIR_MP0  ])[knw  ];
      //(Dbuff.f[DIR_P0P  ])[k] = (D.f[DIR_P0P  ])[kte  ];
      //(Dbuff.f[DIR_M0M  ])[k] = (D.f[DIR_M0M  ])[kbw  ];
      //(Dbuff.f[DIR_P0M  ])[k] = (D.f[DIR_P0M  ])[kbe  ];
      //(Dbuff.f[DIR_M0P  ])[k] = (D.f[DIR_M0P  ])[ktw  ];
      //(Dbuff.f[DIR_0PP  ])[k] = (D.f[DIR_0PP  ])[ktn  ];
      //(Dbuff.f[DIR_0MM  ])[k] = (D.f[DIR_0MM  ])[kbs  ];
      //(Dbuff.f[DIR_0PM  ])[k] = (D.f[DIR_0PM  ])[kbn  ];
      //(Dbuff.f[DIR_0MP  ])[k] = (D.f[DIR_0MP  ])[kts  ];
      //(Dbuff.f[DIR_000])[k] = (D.f[DIR_000])[kzero];
      //(Dbuff.f[DIR_PPP ])[k] = (D.f[DIR_PPP ])[ktne ];
      //(Dbuff.f[DIR_MMP ])[k] = (D.f[DIR_MMP ])[ktsw ];
      //(Dbuff.f[DIR_PMP ])[k] = (D.f[DIR_PMP ])[ktse ];
      //(Dbuff.f[DIR_MPP ])[k] = (D.f[DIR_MPP ])[ktnw ];
      //(Dbuff.f[DIR_PPM ])[k] = (D.f[DIR_PPM ])[kbne ];
      //(Dbuff.f[DIR_MMM ])[k] = (D.f[DIR_MMM ])[kbsw ];
      //(Dbuff.f[DIR_PMM ])[k] = (D.f[DIR_PMM ])[kbse ];
      //(Dbuff.f[DIR_MPM ])[k] = (D.f[DIR_MPM ])[kbnw ];
      (Dbuff.f[DIR_P00   ])[k] = (D.f[DIR_M00   ])[kw   ];
      (Dbuff.f[DIR_M00   ])[k] = (D.f[DIR_P00   ])[ke   ];
      (Dbuff.f[DIR_0P0   ])[k] = (D.f[DIR_0M0   ])[ks   ];
      (Dbuff.f[DIR_0M0   ])[k] = (D.f[DIR_0P0   ])[kn   ];
      (Dbuff.f[DIR_00P   ])[k] = (D.f[DIR_00M   ])[kb   ];
      (Dbuff.f[DIR_00M   ])[k] = (D.f[DIR_00P   ])[kt   ];
      (Dbuff.f[DIR_PP0  ])[k] = (D.f[DIR_MM0  ])[ksw  ];
      (Dbuff.f[DIR_MM0  ])[k] = (D.f[DIR_PP0  ])[kne  ];
      (Dbuff.f[DIR_PM0  ])[k] = (D.f[DIR_MP0  ])[knw  ];
      (Dbuff.f[DIR_MP0  ])[k] = (D.f[DIR_PM0  ])[kse  ];
      (Dbuff.f[DIR_P0P  ])[k] = (D.f[DIR_M0M  ])[kbw  ];
      (Dbuff.f[DIR_M0M  ])[k] = (D.f[DIR_P0P  ])[kte  ];
      (Dbuff.f[DIR_P0M  ])[k] = (D.f[DIR_M0P  ])[ktw  ];
      (Dbuff.f[DIR_M0P  ])[k] = (D.f[DIR_P0M  ])[kbe  ];
      (Dbuff.f[DIR_0PP  ])[k] = (D.f[DIR_0MM  ])[kbs  ];
      (Dbuff.f[DIR_0MM  ])[k] = (D.f[DIR_0PP  ])[ktn  ];
      (Dbuff.f[DIR_0PM  ])[k] = (D.f[DIR_0MP  ])[kts  ];
      (Dbuff.f[DIR_0MP  ])[k] = (D.f[DIR_0PM  ])[kbn  ];
      (Dbuff.f[DIR_000])[k] = (D.f[DIR_000])[kzero];
      (Dbuff.f[DIR_PPP ])[k] = (D.f[DIR_MMM ])[kbsw ];
      (Dbuff.f[DIR_MMP ])[k] = (D.f[DIR_PPM ])[kbne ];
      (Dbuff.f[DIR_PMP ])[k] = (D.f[DIR_MPM ])[kbnw ];
      (Dbuff.f[DIR_MPP ])[k] = (D.f[DIR_PMM ])[kbse ];
      (Dbuff.f[DIR_PPM ])[k] = (D.f[DIR_MMP ])[ktsw ];
      (Dbuff.f[DIR_MMM ])[k] = (D.f[DIR_PPP ])[ktne ];
      (Dbuff.f[DIR_PMM ])[k] = (D.f[DIR_MPP ])[ktnw ];
      (Dbuff.f[DIR_MPM ])[k] = (D.f[DIR_PMP ])[ktse ];
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
         D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
      } 
      else
      {
         D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[DIR_P00   ] = &bufferFs[DIR_P00   *buffmax];
      Dbuff.f[DIR_M00   ] = &bufferFs[DIR_M00   *buffmax];
      Dbuff.f[DIR_0P0   ] = &bufferFs[DIR_0P0   *buffmax];
      Dbuff.f[DIR_0M0   ] = &bufferFs[DIR_0M0   *buffmax];
      Dbuff.f[DIR_00P   ] = &bufferFs[DIR_00P   *buffmax];
      Dbuff.f[DIR_00M   ] = &bufferFs[DIR_00M   *buffmax];
      Dbuff.f[DIR_PP0  ] = &bufferFs[DIR_PP0  *buffmax];
      Dbuff.f[DIR_MM0  ] = &bufferFs[DIR_MM0  *buffmax];
      Dbuff.f[DIR_PM0  ] = &bufferFs[DIR_PM0  *buffmax];
      Dbuff.f[DIR_MP0  ] = &bufferFs[DIR_MP0  *buffmax];
      Dbuff.f[DIR_P0P  ] = &bufferFs[DIR_P0P  *buffmax];
      Dbuff.f[DIR_M0M  ] = &bufferFs[DIR_M0M  *buffmax];
      Dbuff.f[DIR_P0M  ] = &bufferFs[DIR_P0M  *buffmax];
      Dbuff.f[DIR_M0P  ] = &bufferFs[DIR_M0P  *buffmax];
      Dbuff.f[DIR_0PP  ] = &bufferFs[DIR_0PP  *buffmax];
      Dbuff.f[DIR_0MM  ] = &bufferFs[DIR_0MM  *buffmax];
      Dbuff.f[DIR_0PM  ] = &bufferFs[DIR_0PM  *buffmax];
      Dbuff.f[DIR_0MP  ] = &bufferFs[DIR_0MP  *buffmax];
      Dbuff.f[DIR_000] = &bufferFs[DIR_000*buffmax];
      Dbuff.f[DIR_PPP ] = &bufferFs[DIR_PPP *buffmax];
      Dbuff.f[DIR_MMP ] = &bufferFs[DIR_MMP *buffmax];
      Dbuff.f[DIR_PMP ] = &bufferFs[DIR_PMP *buffmax];
      Dbuff.f[DIR_MPP ] = &bufferFs[DIR_MPP *buffmax];
      Dbuff.f[DIR_PPM ] = &bufferFs[DIR_PPM *buffmax];
      Dbuff.f[DIR_MMM ] = &bufferFs[DIR_MMM *buffmax];
      Dbuff.f[DIR_PMM ] = &bufferFs[DIR_PMM *buffmax];
      Dbuff.f[DIR_MPM ] = &bufferFs[DIR_MPM *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      //(D.f[DIR_P00   ])[ke   ] = (Dbuff.f[DIR_P00   ])[k];
      //(D.f[DIR_M00   ])[kw   ] = (Dbuff.f[DIR_M00   ])[k];
      //(D.f[DIR_0P0   ])[kn   ] = (Dbuff.f[DIR_0P0   ])[k];
      //(D.f[DIR_0M0   ])[ks   ] = (Dbuff.f[DIR_0M0   ])[k];
      //(D.f[DIR_00P   ])[kt   ] = (Dbuff.f[DIR_00P   ])[k];
      //(D.f[DIR_00M   ])[kb   ] = (Dbuff.f[DIR_00M   ])[k];
      //(D.f[DIR_PP0  ])[kne  ] = (Dbuff.f[DIR_PP0  ])[k];
      //(D.f[DIR_MM0  ])[ksw  ] = (Dbuff.f[DIR_MM0  ])[k];
      //(D.f[DIR_PM0  ])[kse  ] = (Dbuff.f[DIR_PM0  ])[k];
      //(D.f[DIR_MP0  ])[knw  ] = (Dbuff.f[DIR_MP0  ])[k];
      //(D.f[DIR_P0P  ])[kte  ] = (Dbuff.f[DIR_P0P  ])[k];
      //(D.f[DIR_M0M  ])[kbw  ] = (Dbuff.f[DIR_M0M  ])[k];
      //(D.f[DIR_P0M  ])[kbe  ] = (Dbuff.f[DIR_P0M  ])[k];
      //(D.f[DIR_M0P  ])[ktw  ] = (Dbuff.f[DIR_M0P  ])[k];
      //(D.f[DIR_0PP  ])[ktn  ] = (Dbuff.f[DIR_0PP  ])[k];
      //(D.f[DIR_0MM  ])[kbs  ] = (Dbuff.f[DIR_0MM  ])[k];
      //(D.f[DIR_0PM  ])[kbn  ] = (Dbuff.f[DIR_0PM  ])[k];
      //(D.f[DIR_0MP  ])[kts  ] = (Dbuff.f[DIR_0MP  ])[k];
      //(D.f[DIR_000])[kzero] = (Dbuff.f[DIR_000])[k];
      //(D.f[DIR_PPP ])[ktne ] = (Dbuff.f[DIR_PPP ])[k];
      //(D.f[DIR_MMP ])[ktsw ] = (Dbuff.f[DIR_MMP ])[k];
      //(D.f[DIR_PMP ])[ktse ] = (Dbuff.f[DIR_PMP ])[k];
      //(D.f[DIR_MPP ])[ktnw ] = (Dbuff.f[DIR_MPP ])[k];
      //(D.f[DIR_PPM ])[kbne ] = (Dbuff.f[DIR_PPM ])[k];
      //(D.f[DIR_MMM ])[kbsw ] = (Dbuff.f[DIR_MMM ])[k];
      //(D.f[DIR_PMM ])[kbse ] = (Dbuff.f[DIR_PMM ])[k];
      //(D.f[DIR_MPM ])[kbnw ] = (Dbuff.f[DIR_MPM ])[k];
      (D.f[DIR_M00   ])[kw   ] = (Dbuff.f[DIR_P00   ])[k];
      (D.f[DIR_P00   ])[ke   ] = (Dbuff.f[DIR_M00   ])[k];
      (D.f[DIR_0M0   ])[ks   ] = (Dbuff.f[DIR_0P0   ])[k];
      (D.f[DIR_0P0   ])[kn   ] = (Dbuff.f[DIR_0M0   ])[k];
      (D.f[DIR_00M   ])[kb   ] = (Dbuff.f[DIR_00P   ])[k];
      (D.f[DIR_00P   ])[kt   ] = (Dbuff.f[DIR_00M   ])[k];
      (D.f[DIR_MM0  ])[ksw  ] = (Dbuff.f[DIR_PP0  ])[k];
      (D.f[DIR_PP0  ])[kne  ] = (Dbuff.f[DIR_MM0  ])[k];
      (D.f[DIR_MP0  ])[knw  ] = (Dbuff.f[DIR_PM0  ])[k];
      (D.f[DIR_PM0  ])[kse  ] = (Dbuff.f[DIR_MP0  ])[k];
      (D.f[DIR_M0M  ])[kbw  ] = (Dbuff.f[DIR_P0P  ])[k];
      (D.f[DIR_P0P  ])[kte  ] = (Dbuff.f[DIR_M0M  ])[k];
      (D.f[DIR_M0P  ])[ktw  ] = (Dbuff.f[DIR_P0M  ])[k];
      (D.f[DIR_P0M  ])[kbe  ] = (Dbuff.f[DIR_M0P  ])[k];
      (D.f[DIR_0MM  ])[kbs  ] = (Dbuff.f[DIR_0PP  ])[k];
      (D.f[DIR_0PP  ])[ktn  ] = (Dbuff.f[DIR_0MM  ])[k];
      (D.f[DIR_0MP  ])[kts  ] = (Dbuff.f[DIR_0PM  ])[k];
      (D.f[DIR_0PM  ])[kbn  ] = (Dbuff.f[DIR_0MP  ])[k];
      (D.f[DIR_000])[kzero] = (Dbuff.f[DIR_000])[k];
      (D.f[DIR_MMM ])[kbsw ] = (Dbuff.f[DIR_PPP ])[k];
      (D.f[DIR_PPM ])[kbne ] = (Dbuff.f[DIR_MMP ])[k];
      (D.f[DIR_MPM ])[kbnw ] = (Dbuff.f[DIR_PMP ])[k];
      (D.f[DIR_PMM ])[kbse ] = (Dbuff.f[DIR_MPP ])[k];
      (D.f[DIR_MMP ])[ktsw ] = (Dbuff.f[DIR_PPM ])[k];
      (D.f[DIR_PPP ])[ktne ] = (Dbuff.f[DIR_MMM ])[k];
      (D.f[DIR_MPP ])[ktnw ] = (Dbuff.f[DIR_PMM ])[k];
      (D.f[DIR_PMP ])[ktse ] = (Dbuff.f[DIR_MPM ])[k];
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
         D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
      } 
      else
      {
         D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[DIR_P00   ] = &bufferFs[DIR_P00   *buffmax];
      Dbuff.f[DIR_M00   ] = &bufferFs[DIR_M00   *buffmax];
      Dbuff.f[DIR_0P0   ] = &bufferFs[DIR_0P0   *buffmax];
      Dbuff.f[DIR_0M0   ] = &bufferFs[DIR_0M0   *buffmax];
      Dbuff.f[DIR_00P   ] = &bufferFs[DIR_00P   *buffmax];
      Dbuff.f[DIR_00M   ] = &bufferFs[DIR_00M   *buffmax];
      Dbuff.f[DIR_PP0  ] = &bufferFs[DIR_PP0  *buffmax];
      Dbuff.f[DIR_MM0  ] = &bufferFs[DIR_MM0  *buffmax];
      Dbuff.f[DIR_PM0  ] = &bufferFs[DIR_PM0  *buffmax];
      Dbuff.f[DIR_MP0  ] = &bufferFs[DIR_MP0  *buffmax];
      Dbuff.f[DIR_P0P  ] = &bufferFs[DIR_P0P  *buffmax];
      Dbuff.f[DIR_M0M  ] = &bufferFs[DIR_M0M  *buffmax];
      Dbuff.f[DIR_P0M  ] = &bufferFs[DIR_P0M  *buffmax];
      Dbuff.f[DIR_M0P  ] = &bufferFs[DIR_M0P  *buffmax];
      Dbuff.f[DIR_0PP  ] = &bufferFs[DIR_0PP  *buffmax];
      Dbuff.f[DIR_0MM  ] = &bufferFs[DIR_0MM  *buffmax];
      Dbuff.f[DIR_0PM  ] = &bufferFs[DIR_0PM  *buffmax];
      Dbuff.f[DIR_0MP  ] = &bufferFs[DIR_0MP  *buffmax];
      Dbuff.f[DIR_000] = &bufferFs[DIR_000*buffmax];
      Dbuff.f[DIR_PPP ] = &bufferFs[DIR_PPP *buffmax];
      Dbuff.f[DIR_MMP ] = &bufferFs[DIR_MMP *buffmax];
      Dbuff.f[DIR_PMP ] = &bufferFs[DIR_PMP *buffmax];
      Dbuff.f[DIR_MPP ] = &bufferFs[DIR_MPP *buffmax];
      Dbuff.f[DIR_PPM ] = &bufferFs[DIR_PPM *buffmax];
      Dbuff.f[DIR_MMM ] = &bufferFs[DIR_MMM *buffmax];
      Dbuff.f[DIR_PMM ] = &bufferFs[DIR_PMM *buffmax];
      Dbuff.f[DIR_MPM ] = &bufferFs[DIR_MPM *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      (Dbuff.f[DIR_P00   ])[k] = (D.f[DIR_P00   ])[ke   ];
      (Dbuff.f[DIR_M00   ])[k] = (D.f[DIR_M00   ])[kw   ];
      (Dbuff.f[DIR_0P0   ])[k] = (D.f[DIR_0P0   ])[kn   ];
      (Dbuff.f[DIR_0M0   ])[k] = (D.f[DIR_0M0   ])[ks   ];
      (Dbuff.f[DIR_00P   ])[k] = (D.f[DIR_00P   ])[kt   ];
      (Dbuff.f[DIR_00M   ])[k] = (D.f[DIR_00M   ])[kb   ];
      (Dbuff.f[DIR_PP0  ])[k] = (D.f[DIR_PP0  ])[kne  ];
      (Dbuff.f[DIR_MM0  ])[k] = (D.f[DIR_MM0  ])[ksw  ];
      (Dbuff.f[DIR_PM0  ])[k] = (D.f[DIR_PM0  ])[kse  ];
      (Dbuff.f[DIR_MP0  ])[k] = (D.f[DIR_MP0  ])[knw  ];
      (Dbuff.f[DIR_P0P  ])[k] = (D.f[DIR_P0P  ])[kte  ];
      (Dbuff.f[DIR_M0M  ])[k] = (D.f[DIR_M0M  ])[kbw  ];
      (Dbuff.f[DIR_P0M  ])[k] = (D.f[DIR_P0M  ])[kbe  ];
      (Dbuff.f[DIR_M0P  ])[k] = (D.f[DIR_M0P  ])[ktw  ];
      (Dbuff.f[DIR_0PP  ])[k] = (D.f[DIR_0PP  ])[ktn  ];
      (Dbuff.f[DIR_0MM  ])[k] = (D.f[DIR_0MM  ])[kbs  ];
      (Dbuff.f[DIR_0PM  ])[k] = (D.f[DIR_0PM  ])[kbn  ];
      (Dbuff.f[DIR_0MP  ])[k] = (D.f[DIR_0MP  ])[kts  ];
      (Dbuff.f[DIR_000])[k] = (D.f[DIR_000])[kzero];
      (Dbuff.f[DIR_PPP ])[k] = (D.f[DIR_PPP ])[ktne ];
      (Dbuff.f[DIR_MMP ])[k] = (D.f[DIR_MMP ])[ktsw ];
      (Dbuff.f[DIR_PMP ])[k] = (D.f[DIR_PMP ])[ktse ];
      (Dbuff.f[DIR_MPP ])[k] = (D.f[DIR_MPP ])[ktnw ];
      (Dbuff.f[DIR_PPM ])[k] = (D.f[DIR_PPM ])[kbne ];
      (Dbuff.f[DIR_MMM ])[k] = (D.f[DIR_MMM ])[kbsw ];
      (Dbuff.f[DIR_PMM ])[k] = (D.f[DIR_PMM ])[kbse ];
      (Dbuff.f[DIR_MPM ])[k] = (D.f[DIR_MPM ])[kbnw ];
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
         D.f[DIR_P00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_M00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_PMP *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_MPM *size_Mat];
      } 
      else
      {
         D.f[DIR_M00   ] = &DD[DIR_P00   *size_Mat];
         D.f[DIR_P00   ] = &DD[DIR_M00   *size_Mat];
         D.f[DIR_0M0   ] = &DD[DIR_0P0   *size_Mat];
         D.f[DIR_0P0   ] = &DD[DIR_0M0   *size_Mat];
         D.f[DIR_00M   ] = &DD[DIR_00P   *size_Mat];
         D.f[DIR_00P   ] = &DD[DIR_00M   *size_Mat];
         D.f[DIR_MM0  ] = &DD[DIR_PP0  *size_Mat];
         D.f[DIR_PP0  ] = &DD[DIR_MM0  *size_Mat];
         D.f[DIR_MP0  ] = &DD[DIR_PM0  *size_Mat];
         D.f[DIR_PM0  ] = &DD[DIR_MP0  *size_Mat];
         D.f[DIR_M0M  ] = &DD[DIR_P0P  *size_Mat];
         D.f[DIR_P0P  ] = &DD[DIR_M0M  *size_Mat];
         D.f[DIR_M0P  ] = &DD[DIR_P0M  *size_Mat];
         D.f[DIR_P0M  ] = &DD[DIR_M0P  *size_Mat];
         D.f[DIR_0MM  ] = &DD[DIR_0PP  *size_Mat];
         D.f[DIR_0PP  ] = &DD[DIR_0MM  *size_Mat];
         D.f[DIR_0MP  ] = &DD[DIR_0PM  *size_Mat];
         D.f[DIR_0PM  ] = &DD[DIR_0MP  *size_Mat];
         D.f[DIR_000] = &DD[DIR_000*size_Mat];
         D.f[DIR_PPP ] = &DD[DIR_MMM *size_Mat];
         D.f[DIR_MMP ] = &DD[DIR_PPM *size_Mat];
         D.f[DIR_PMP ] = &DD[DIR_MPM *size_Mat];
         D.f[DIR_MPP ] = &DD[DIR_PMM *size_Mat];
         D.f[DIR_PPM ] = &DD[DIR_MMP *size_Mat];
         D.f[DIR_MMM ] = &DD[DIR_PPP *size_Mat];
         D.f[DIR_PMM ] = &DD[DIR_MPP *size_Mat];
         D.f[DIR_MPM ] = &DD[DIR_PMP *size_Mat];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[DIR_P00   ] = &bufferFs[DIR_P00   *buffmax];
      Dbuff.f[DIR_M00   ] = &bufferFs[DIR_M00   *buffmax];
      Dbuff.f[DIR_0P0   ] = &bufferFs[DIR_0P0   *buffmax];
      Dbuff.f[DIR_0M0   ] = &bufferFs[DIR_0M0   *buffmax];
      Dbuff.f[DIR_00P   ] = &bufferFs[DIR_00P   *buffmax];
      Dbuff.f[DIR_00M   ] = &bufferFs[DIR_00M   *buffmax];
      Dbuff.f[DIR_PP0  ] = &bufferFs[DIR_PP0  *buffmax];
      Dbuff.f[DIR_MM0  ] = &bufferFs[DIR_MM0  *buffmax];
      Dbuff.f[DIR_PM0  ] = &bufferFs[DIR_PM0  *buffmax];
      Dbuff.f[DIR_MP0  ] = &bufferFs[DIR_MP0  *buffmax];
      Dbuff.f[DIR_P0P  ] = &bufferFs[DIR_P0P  *buffmax];
      Dbuff.f[DIR_M0M  ] = &bufferFs[DIR_M0M  *buffmax];
      Dbuff.f[DIR_P0M  ] = &bufferFs[DIR_P0M  *buffmax];
      Dbuff.f[DIR_M0P  ] = &bufferFs[DIR_M0P  *buffmax];
      Dbuff.f[DIR_0PP  ] = &bufferFs[DIR_0PP  *buffmax];
      Dbuff.f[DIR_0MM  ] = &bufferFs[DIR_0MM  *buffmax];
      Dbuff.f[DIR_0PM  ] = &bufferFs[DIR_0PM  *buffmax];
      Dbuff.f[DIR_0MP  ] = &bufferFs[DIR_0MP  *buffmax];
      Dbuff.f[DIR_000] = &bufferFs[DIR_000*buffmax];
      Dbuff.f[DIR_PPP ] = &bufferFs[DIR_PPP *buffmax];
      Dbuff.f[DIR_MMP ] = &bufferFs[DIR_MMP *buffmax];
      Dbuff.f[DIR_PMP ] = &bufferFs[DIR_PMP *buffmax];
      Dbuff.f[DIR_MPP ] = &bufferFs[DIR_MPP *buffmax];
      Dbuff.f[DIR_PPM ] = &bufferFs[DIR_PPM *buffmax];
      Dbuff.f[DIR_MMM ] = &bufferFs[DIR_MMM *buffmax];
      Dbuff.f[DIR_PMM ] = &bufferFs[DIR_PMM *buffmax];
      Dbuff.f[DIR_MPM ] = &bufferFs[DIR_MPM *buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      (D.f[DIR_P00   ])[ke   ] = (Dbuff.f[DIR_P00   ])[k];
      (D.f[DIR_M00   ])[kw   ] = (Dbuff.f[DIR_M00   ])[k];
      (D.f[DIR_0P0   ])[kn   ] = (Dbuff.f[DIR_0P0   ])[k];
      (D.f[DIR_0M0   ])[ks   ] = (Dbuff.f[DIR_0M0   ])[k];
      (D.f[DIR_00P   ])[kt   ] = (Dbuff.f[DIR_00P   ])[k];
      (D.f[DIR_00M   ])[kb   ] = (Dbuff.f[DIR_00M   ])[k];
      (D.f[DIR_PP0  ])[kne  ] = (Dbuff.f[DIR_PP0  ])[k];
      (D.f[DIR_MM0  ])[ksw  ] = (Dbuff.f[DIR_MM0  ])[k];
      (D.f[DIR_PM0  ])[kse  ] = (Dbuff.f[DIR_PM0  ])[k];
      (D.f[DIR_MP0  ])[knw  ] = (Dbuff.f[DIR_MP0  ])[k];
      (D.f[DIR_P0P  ])[kte  ] = (Dbuff.f[DIR_P0P  ])[k];
      (D.f[DIR_M0M  ])[kbw  ] = (Dbuff.f[DIR_M0M  ])[k];
      (D.f[DIR_P0M  ])[kbe  ] = (Dbuff.f[DIR_P0M  ])[k];
      (D.f[DIR_M0P  ])[ktw  ] = (Dbuff.f[DIR_M0P  ])[k];
      (D.f[DIR_0PP  ])[ktn  ] = (Dbuff.f[DIR_0PP  ])[k];
      (D.f[DIR_0MM  ])[kbs  ] = (Dbuff.f[DIR_0MM  ])[k];
      (D.f[DIR_0PM  ])[kbn  ] = (Dbuff.f[DIR_0PM  ])[k];
      (D.f[DIR_0MP  ])[kts  ] = (Dbuff.f[DIR_0MP  ])[k];
      (D.f[DIR_000])[kzero] = (Dbuff.f[DIR_000])[k];
      (D.f[DIR_PPP ])[ktne ] = (Dbuff.f[DIR_PPP ])[k];
      (D.f[DIR_MMP ])[ktsw ] = (Dbuff.f[DIR_MMP ])[k];
      (D.f[DIR_PMP ])[ktse ] = (Dbuff.f[DIR_PMP ])[k];
      (D.f[DIR_MPP ])[ktnw ] = (Dbuff.f[DIR_MPP ])[k];
      (D.f[DIR_PPM ])[kbne ] = (Dbuff.f[DIR_PPM ])[k];
      (D.f[DIR_MMM ])[kbsw ] = (Dbuff.f[DIR_MMM ])[k];
      (D.f[DIR_PMM ])[kbse ] = (Dbuff.f[DIR_PMM ])[k];
      (D.f[DIR_MPM ])[kbnw ] = (Dbuff.f[DIR_MPM ])[k];
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
			G.g[DIR_P00] = &G6[DIR_P00   *size_Mat];
			G.g[DIR_M00] = &G6[DIR_M00   *size_Mat];
			G.g[DIR_0P0] = &G6[DIR_0P0   *size_Mat];
			G.g[DIR_0M0] = &G6[DIR_0M0   *size_Mat];
			G.g[DIR_00P] = &G6[DIR_00P   *size_Mat];
			G.g[DIR_00M] = &G6[DIR_00M   *size_Mat];
		}
		else
		{
			G.g[DIR_M00] = &G6[DIR_P00   *size_Mat];
			G.g[DIR_P00] = &G6[DIR_M00   *size_Mat];
			G.g[DIR_0M0] = &G6[DIR_0P0   *size_Mat];
			G.g[DIR_0P0] = &G6[DIR_0M0   *size_Mat];
			G.g[DIR_00M] = &G6[DIR_00P   *size_Mat];
			G.g[DIR_00P] = &G6[DIR_00M   *size_Mat];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Buffer Gs
		Distributions6 Dbuff;
		Dbuff.g[DIR_P00] = &bufferGs[DIR_P00   *buffmax];
		Dbuff.g[DIR_M00] = &bufferGs[DIR_M00   *buffmax];
		Dbuff.g[DIR_0P0] = &bufferGs[DIR_0P0   *buffmax];
		Dbuff.g[DIR_0M0] = &bufferGs[DIR_0M0   *buffmax];
		Dbuff.g[DIR_00P] = &bufferGs[DIR_00P   *buffmax];
		Dbuff.g[DIR_00M] = &bufferGs[DIR_00M   *buffmax];
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//write Gs to buffer
		(Dbuff.g[DIR_P00])[k] = (G.g[DIR_M00])[kw];
		(Dbuff.g[DIR_M00])[k] = (G.g[DIR_P00])[kr];
		(Dbuff.g[DIR_0P0])[k] = (G.g[DIR_0M0])[ks];
		(Dbuff.g[DIR_0M0])[k] = (G.g[DIR_0P0])[kr];
		(Dbuff.g[DIR_00P])[k] = (G.g[DIR_00M])[kb];
		(Dbuff.g[DIR_00M])[k] = (G.g[DIR_00P])[kr];
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
			G.g[DIR_P00] = &G6[DIR_P00   *size_Mat];
			G.g[DIR_M00] = &G6[DIR_M00   *size_Mat];
			G.g[DIR_0P0] = &G6[DIR_0P0   *size_Mat];
			G.g[DIR_0M0] = &G6[DIR_0M0   *size_Mat];
			G.g[DIR_00P] = &G6[DIR_00P   *size_Mat];
			G.g[DIR_00M] = &G6[DIR_00M   *size_Mat];
		}
		else
		{
			G.g[DIR_M00] = &G6[DIR_P00   *size_Mat];
			G.g[DIR_P00] = &G6[DIR_M00   *size_Mat];
			G.g[DIR_0M0] = &G6[DIR_0P0   *size_Mat];
			G.g[DIR_0P0] = &G6[DIR_0M0   *size_Mat];
			G.g[DIR_00M] = &G6[DIR_00P   *size_Mat];
			G.g[DIR_00P] = &G6[DIR_00M   *size_Mat];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Buffer Gs
		Distributions6 Dbuff;
		Dbuff.g[DIR_P00] = &bufferGs[DIR_P00   *buffmax];
		Dbuff.g[DIR_M00] = &bufferGs[DIR_M00   *buffmax];
		Dbuff.g[DIR_0P0] = &bufferGs[DIR_0P0   *buffmax];
		Dbuff.g[DIR_0M0] = &bufferGs[DIR_0M0   *buffmax];
		Dbuff.g[DIR_00P] = &bufferGs[DIR_00P   *buffmax];
		Dbuff.g[DIR_00M] = &bufferGs[DIR_00M   *buffmax];
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//write buffer to Gs
		(G.g[DIR_M00])[kw] = (Dbuff.g[DIR_P00])[k];
		(G.g[DIR_P00])[kr] = (Dbuff.g[DIR_M00])[k];
		(G.g[DIR_0M0])[ks] = (Dbuff.g[DIR_0P0])[k];
		(G.g[DIR_0P0])[kr] = (Dbuff.g[DIR_0M0])[k];
		(G.g[DIR_00M])[kb] = (Dbuff.g[DIR_00P])[k];
		(G.g[DIR_00P])[kr] = (Dbuff.g[DIR_00M])[k];
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
