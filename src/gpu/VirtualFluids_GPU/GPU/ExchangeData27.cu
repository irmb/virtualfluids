/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void getSendFsPost27(real* DD,
										   real* bufferFs,
										   int* sendIndex,
                                           int buffmax,
                                           unsigned int* neighborX,
                                           unsigned int* neighborY,
                                           unsigned int* neighborZ,
                                           unsigned long long numberOfLBnodes, 
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
         D.f[dP00] = &DD[dP00 * numberOfLBnodes];
         D.f[dM00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D.f[dM00] = &DD[dP00 * numberOfLBnodes];
         D.f[dP00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dP00] = &bufferFs[dP00 * buffmax];
      Dbuff.f[dM00] = &bufferFs[dM00 * buffmax];
      Dbuff.f[DIR_0P0] = &bufferFs[DIR_0P0 * buffmax];
      Dbuff.f[DIR_0M0] = &bufferFs[DIR_0M0 * buffmax];
      Dbuff.f[DIR_00P] = &bufferFs[DIR_00P * buffmax];
      Dbuff.f[DIR_00M] = &bufferFs[DIR_00M * buffmax];
      Dbuff.f[DIR_PP0] = &bufferFs[DIR_PP0 * buffmax];
      Dbuff.f[DIR_MM0] = &bufferFs[DIR_MM0 * buffmax];
      Dbuff.f[DIR_PM0] = &bufferFs[DIR_PM0 * buffmax];
      Dbuff.f[DIR_MP0] = &bufferFs[DIR_MP0 * buffmax];
      Dbuff.f[DIR_P0P] = &bufferFs[DIR_P0P * buffmax];
      Dbuff.f[DIR_M0M] = &bufferFs[DIR_M0M * buffmax];
      Dbuff.f[DIR_P0M] = &bufferFs[DIR_P0M * buffmax];
      Dbuff.f[DIR_M0P] = &bufferFs[DIR_M0P * buffmax];
      Dbuff.f[DIR_0PP] = &bufferFs[DIR_0PP * buffmax];
      Dbuff.f[DIR_0MM] = &bufferFs[DIR_0MM * buffmax];
      Dbuff.f[DIR_0PM] = &bufferFs[DIR_0PM * buffmax];
      Dbuff.f[DIR_0MP] = &bufferFs[DIR_0MP * buffmax];
      Dbuff.f[d000] = &bufferFs[d000 * buffmax];
      Dbuff.f[DIR_PPP] = &bufferFs[DIR_PPP * buffmax];
      Dbuff.f[DIR_MMP] = &bufferFs[DIR_MMP * buffmax];
      Dbuff.f[DIR_PMP] = &bufferFs[DIR_PMP * buffmax];
      Dbuff.f[DIR_MPP] = &bufferFs[DIR_MPP * buffmax];
      Dbuff.f[DIR_PPM] = &bufferFs[DIR_PPM * buffmax];
      Dbuff.f[DIR_MMM] = &bufferFs[DIR_MMM * buffmax];
      Dbuff.f[DIR_PMM] = &bufferFs[DIR_PMM * buffmax];
      Dbuff.f[DIR_MPM] = &bufferFs[DIR_MPM * buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      //(Dbuff.f[dP00])[k] = (D.f[dP00])[ke   ];
      //(Dbuff.f[dM00])[k] = (D.f[dM00])[kw   ];
      //(Dbuff.f[DIR_0P0])[k] = (D.f[DIR_0P0])[kn   ];
      //(Dbuff.f[DIR_0M0])[k] = (D.f[DIR_0M0])[ks   ];
      //(Dbuff.f[DIR_00P])[k] = (D.f[DIR_00P])[kt   ];
      //(Dbuff.f[DIR_00M])[k] = (D.f[DIR_00M])[kb   ];
      //(Dbuff.f[DIR_PP0])[k] = (D.f[DIR_PP0])[kne  ];
      //(Dbuff.f[DIR_MM0])[k] = (D.f[DIR_MM0])[ksw  ];
      //(Dbuff.f[DIR_PM0])[k] = (D.f[DIR_PM0])[kse  ];
      //(Dbuff.f[DIR_MP0])[k] = (D.f[DIR_MP0])[knw  ];
      //(Dbuff.f[DIR_P0P])[k] = (D.f[DIR_P0P])[kte  ];
      //(Dbuff.f[DIR_M0M])[k] = (D.f[DIR_M0M])[kbw  ];
      //(Dbuff.f[DIR_P0M])[k] = (D.f[DIR_P0M])[kbe  ];
      //(Dbuff.f[DIR_M0P])[k] = (D.f[DIR_M0P])[ktw  ];
      //(Dbuff.f[DIR_0PP])[k] = (D.f[DIR_0PP])[ktn  ];
      //(Dbuff.f[DIR_0MM])[k] = (D.f[DIR_0MM])[kbs  ];
      //(Dbuff.f[DIR_0PM])[k] = (D.f[DIR_0PM])[kbn  ];
      //(Dbuff.f[DIR_0MP])[k] = (D.f[DIR_0MP])[kts  ];
      //(Dbuff.f[d000])[k] = (D.f[d000])[kzero];
      //(Dbuff.f[DIR_PPP])[k] = (D.f[DIR_PPP])[ktne ];
      //(Dbuff.f[DIR_MMP])[k] = (D.f[DIR_MMP])[ktsw ];
      //(Dbuff.f[DIR_PMP])[k] = (D.f[DIR_PMP])[ktse ];
      //(Dbuff.f[DIR_MPP])[k] = (D.f[DIR_MPP])[ktnw ];
      //(Dbuff.f[DIR_PPM])[k] = (D.f[DIR_PPM])[kbne ];
      //(Dbuff.f[DIR_MMM])[k] = (D.f[DIR_MMM])[kbsw ];
      //(Dbuff.f[DIR_PMM])[k] = (D.f[DIR_PMM])[kbse ];
      //(Dbuff.f[DIR_MPM])[k] = (D.f[DIR_MPM])[kbnw ];
      (Dbuff.f[dP00])[k] = (D.f[dM00])[kw   ];
      (Dbuff.f[dM00])[k] = (D.f[dP00])[ke   ];
      (Dbuff.f[DIR_0P0])[k] = (D.f[DIR_0M0])[ks   ];
      (Dbuff.f[DIR_0M0])[k] = (D.f[DIR_0P0])[kn   ];
      (Dbuff.f[DIR_00P])[k] = (D.f[DIR_00M])[kb   ];
      (Dbuff.f[DIR_00M])[k] = (D.f[DIR_00P])[kt   ];
      (Dbuff.f[DIR_PP0])[k] = (D.f[DIR_MM0])[ksw  ];
      (Dbuff.f[DIR_MM0])[k] = (D.f[DIR_PP0])[kne  ];
      (Dbuff.f[DIR_PM0])[k] = (D.f[DIR_MP0])[knw  ];
      (Dbuff.f[DIR_MP0])[k] = (D.f[DIR_PM0])[kse  ];
      (Dbuff.f[DIR_P0P])[k] = (D.f[DIR_M0M])[kbw  ];
      (Dbuff.f[DIR_M0M])[k] = (D.f[DIR_P0P])[kte  ];
      (Dbuff.f[DIR_P0M])[k] = (D.f[DIR_M0P])[ktw  ];
      (Dbuff.f[DIR_M0P])[k] = (D.f[DIR_P0M])[kbe  ];
      (Dbuff.f[DIR_0PP])[k] = (D.f[DIR_0MM])[kbs  ];
      (Dbuff.f[DIR_0MM])[k] = (D.f[DIR_0PP])[ktn  ];
      (Dbuff.f[DIR_0PM])[k] = (D.f[DIR_0MP])[kts  ];
      (Dbuff.f[DIR_0MP])[k] = (D.f[DIR_0PM])[kbn  ];
      (Dbuff.f[d000])[k] = (D.f[d000])[kzero];
      (Dbuff.f[DIR_PPP])[k] = (D.f[DIR_MMM])[kbsw ];
      (Dbuff.f[DIR_MMP])[k] = (D.f[DIR_PPM])[kbne ];
      (Dbuff.f[DIR_PMP])[k] = (D.f[DIR_MPM])[kbnw ];
      (Dbuff.f[DIR_MPP])[k] = (D.f[DIR_PMM])[kbse ];
      (Dbuff.f[DIR_PPM])[k] = (D.f[DIR_MMP])[ktsw ];
      (Dbuff.f[DIR_MMM])[k] = (D.f[DIR_PPP])[ktne ];
      (Dbuff.f[DIR_PMM])[k] = (D.f[DIR_MPP])[ktnw ];
      (Dbuff.f[DIR_MPM])[k] = (D.f[DIR_PMP])[ktse ];
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
                                           unsigned long long numberOfLBnodes, 
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
         D.f[dP00] = &DD[dP00 * numberOfLBnodes];
         D.f[dM00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D.f[dM00] = &DD[dP00 * numberOfLBnodes];
         D.f[dP00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dP00] = &bufferFs[dP00 * buffmax];
      Dbuff.f[dM00] = &bufferFs[dM00 * buffmax];
      Dbuff.f[DIR_0P0] = &bufferFs[DIR_0P0 * buffmax];
      Dbuff.f[DIR_0M0] = &bufferFs[DIR_0M0 * buffmax];
      Dbuff.f[DIR_00P] = &bufferFs[DIR_00P * buffmax];
      Dbuff.f[DIR_00M] = &bufferFs[DIR_00M * buffmax];
      Dbuff.f[DIR_PP0] = &bufferFs[DIR_PP0 * buffmax];
      Dbuff.f[DIR_MM0] = &bufferFs[DIR_MM0 * buffmax];
      Dbuff.f[DIR_PM0] = &bufferFs[DIR_PM0 * buffmax];
      Dbuff.f[DIR_MP0] = &bufferFs[DIR_MP0 * buffmax];
      Dbuff.f[DIR_P0P] = &bufferFs[DIR_P0P * buffmax];
      Dbuff.f[DIR_M0M] = &bufferFs[DIR_M0M * buffmax];
      Dbuff.f[DIR_P0M] = &bufferFs[DIR_P0M * buffmax];
      Dbuff.f[DIR_M0P] = &bufferFs[DIR_M0P * buffmax];
      Dbuff.f[DIR_0PP] = &bufferFs[DIR_0PP * buffmax];
      Dbuff.f[DIR_0MM] = &bufferFs[DIR_0MM * buffmax];
      Dbuff.f[DIR_0PM] = &bufferFs[DIR_0PM * buffmax];
      Dbuff.f[DIR_0MP] = &bufferFs[DIR_0MP * buffmax];
      Dbuff.f[d000] = &bufferFs[d000 * buffmax];
      Dbuff.f[DIR_PPP] = &bufferFs[DIR_PPP * buffmax];
      Dbuff.f[DIR_MMP] = &bufferFs[DIR_MMP * buffmax];
      Dbuff.f[DIR_PMP] = &bufferFs[DIR_PMP * buffmax];
      Dbuff.f[DIR_MPP] = &bufferFs[DIR_MPP * buffmax];
      Dbuff.f[DIR_PPM] = &bufferFs[DIR_PPM * buffmax];
      Dbuff.f[DIR_MMM] = &bufferFs[DIR_MMM * buffmax];
      Dbuff.f[DIR_PMM] = &bufferFs[DIR_PMM * buffmax];
      Dbuff.f[DIR_MPM] = &bufferFs[DIR_MPM * buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      //(D.f[dP00])[ke   ] = (Dbuff.f[dP00])[k];
      //(D.f[dM00])[kw   ] = (Dbuff.f[dM00])[k];
      //(D.f[DIR_0P0])[kn   ] = (Dbuff.f[DIR_0P0])[k];
      //(D.f[DIR_0M0])[ks   ] = (Dbuff.f[DIR_0M0])[k];
      //(D.f[DIR_00P])[kt   ] = (Dbuff.f[DIR_00P])[k];
      //(D.f[DIR_00M])[kb   ] = (Dbuff.f[DIR_00M])[k];
      //(D.f[DIR_PP0])[kne  ] = (Dbuff.f[DIR_PP0])[k];
      //(D.f[DIR_MM0])[ksw  ] = (Dbuff.f[DIR_MM0])[k];
      //(D.f[DIR_PM0])[kse  ] = (Dbuff.f[DIR_PM0])[k];
      //(D.f[DIR_MP0])[knw  ] = (Dbuff.f[DIR_MP0])[k];
      //(D.f[DIR_P0P])[kte  ] = (Dbuff.f[DIR_P0P])[k];
      //(D.f[DIR_M0M])[kbw  ] = (Dbuff.f[DIR_M0M])[k];
      //(D.f[DIR_P0M])[kbe  ] = (Dbuff.f[DIR_P0M])[k];
      //(D.f[DIR_M0P])[ktw  ] = (Dbuff.f[DIR_M0P])[k];
      //(D.f[DIR_0PP])[ktn  ] = (Dbuff.f[DIR_0PP])[k];
      //(D.f[DIR_0MM])[kbs  ] = (Dbuff.f[DIR_0MM])[k];
      //(D.f[DIR_0PM])[kbn  ] = (Dbuff.f[DIR_0PM])[k];
      //(D.f[DIR_0MP])[kts  ] = (Dbuff.f[DIR_0MP])[k];
      //(D.f[d000])[kzero] = (Dbuff.f[d000])[k];
      //(D.f[DIR_PPP])[ktne ] = (Dbuff.f[DIR_PPP])[k];
      //(D.f[DIR_MMP])[ktsw ] = (Dbuff.f[DIR_MMP])[k];
      //(D.f[DIR_PMP])[ktse ] = (Dbuff.f[DIR_PMP])[k];
      //(D.f[DIR_MPP])[ktnw ] = (Dbuff.f[DIR_MPP])[k];
      //(D.f[DIR_PPM])[kbne ] = (Dbuff.f[DIR_PPM])[k];
      //(D.f[DIR_MMM])[kbsw ] = (Dbuff.f[DIR_MMM])[k];
      //(D.f[DIR_PMM])[kbse ] = (Dbuff.f[DIR_PMM])[k];
      //(D.f[DIR_MPM])[kbnw ] = (Dbuff.f[DIR_MPM])[k];
      (D.f[dM00])[kw   ] = (Dbuff.f[dP00])[k];
      (D.f[dP00])[ke   ] = (Dbuff.f[dM00])[k];
      (D.f[DIR_0M0])[ks   ] = (Dbuff.f[DIR_0P0])[k];
      (D.f[DIR_0P0])[kn   ] = (Dbuff.f[DIR_0M0])[k];
      (D.f[DIR_00M])[kb   ] = (Dbuff.f[DIR_00P])[k];
      (D.f[DIR_00P])[kt   ] = (Dbuff.f[DIR_00M])[k];
      (D.f[DIR_MM0])[ksw  ] = (Dbuff.f[DIR_PP0])[k];
      (D.f[DIR_PP0])[kne  ] = (Dbuff.f[DIR_MM0])[k];
      (D.f[DIR_MP0])[knw  ] = (Dbuff.f[DIR_PM0])[k];
      (D.f[DIR_PM0])[kse  ] = (Dbuff.f[DIR_MP0])[k];
      (D.f[DIR_M0M])[kbw  ] = (Dbuff.f[DIR_P0P])[k];
      (D.f[DIR_P0P])[kte  ] = (Dbuff.f[DIR_M0M])[k];
      (D.f[DIR_M0P])[ktw  ] = (Dbuff.f[DIR_P0M])[k];
      (D.f[DIR_P0M])[kbe  ] = (Dbuff.f[DIR_M0P])[k];
      (D.f[DIR_0MM])[kbs  ] = (Dbuff.f[DIR_0PP])[k];
      (D.f[DIR_0PP])[ktn  ] = (Dbuff.f[DIR_0MM])[k];
      (D.f[DIR_0MP])[kts  ] = (Dbuff.f[DIR_0PM])[k];
      (D.f[DIR_0PM])[kbn  ] = (Dbuff.f[DIR_0MP])[k];
      (D.f[d000])[kzero] = (Dbuff.f[d000])[k];
      (D.f[DIR_MMM])[kbsw ] = (Dbuff.f[DIR_PPP])[k];
      (D.f[DIR_PPM])[kbne ] = (Dbuff.f[DIR_MMP])[k];
      (D.f[DIR_MPM])[kbnw ] = (Dbuff.f[DIR_PMP])[k];
      (D.f[DIR_PMM])[kbse ] = (Dbuff.f[DIR_MPP])[k];
      (D.f[DIR_MMP])[ktsw ] = (Dbuff.f[DIR_PPM])[k];
      (D.f[DIR_PPP])[ktne ] = (Dbuff.f[DIR_MMM])[k];
      (D.f[DIR_MPP])[ktnw ] = (Dbuff.f[DIR_PMM])[k];
      (D.f[DIR_PMP])[ktse ] = (Dbuff.f[DIR_MPM])[k];
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
                                          unsigned long long numberOfLBnodes, 
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
         D.f[dP00] = &DD[dP00 * numberOfLBnodes];
         D.f[dM00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D.f[dM00] = &DD[dP00 * numberOfLBnodes];
         D.f[dP00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dP00] = &bufferFs[dP00 * buffmax];
      Dbuff.f[dM00] = &bufferFs[dM00 * buffmax];
      Dbuff.f[DIR_0P0] = &bufferFs[DIR_0P0 * buffmax];
      Dbuff.f[DIR_0M0] = &bufferFs[DIR_0M0 * buffmax];
      Dbuff.f[DIR_00P] = &bufferFs[DIR_00P * buffmax];
      Dbuff.f[DIR_00M] = &bufferFs[DIR_00M * buffmax];
      Dbuff.f[DIR_PP0] = &bufferFs[DIR_PP0 * buffmax];
      Dbuff.f[DIR_MM0] = &bufferFs[DIR_MM0 * buffmax];
      Dbuff.f[DIR_PM0] = &bufferFs[DIR_PM0 * buffmax];
      Dbuff.f[DIR_MP0] = &bufferFs[DIR_MP0 * buffmax];
      Dbuff.f[DIR_P0P] = &bufferFs[DIR_P0P * buffmax];
      Dbuff.f[DIR_M0M] = &bufferFs[DIR_M0M * buffmax];
      Dbuff.f[DIR_P0M] = &bufferFs[DIR_P0M * buffmax];
      Dbuff.f[DIR_M0P] = &bufferFs[DIR_M0P * buffmax];
      Dbuff.f[DIR_0PP] = &bufferFs[DIR_0PP * buffmax];
      Dbuff.f[DIR_0MM] = &bufferFs[DIR_0MM * buffmax];
      Dbuff.f[DIR_0PM] = &bufferFs[DIR_0PM * buffmax];
      Dbuff.f[DIR_0MP] = &bufferFs[DIR_0MP * buffmax];
      Dbuff.f[d000] = &bufferFs[d000 * buffmax];
      Dbuff.f[DIR_PPP] = &bufferFs[DIR_PPP * buffmax];
      Dbuff.f[DIR_MMP] = &bufferFs[DIR_MMP * buffmax];
      Dbuff.f[DIR_PMP] = &bufferFs[DIR_PMP * buffmax];
      Dbuff.f[DIR_MPP] = &bufferFs[DIR_MPP * buffmax];
      Dbuff.f[DIR_PPM] = &bufferFs[DIR_PPM * buffmax];
      Dbuff.f[DIR_MMM] = &bufferFs[DIR_MMM * buffmax];
      Dbuff.f[DIR_PMM] = &bufferFs[DIR_PMM * buffmax];
      Dbuff.f[DIR_MPM] = &bufferFs[DIR_MPM * buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy to buffer
      (Dbuff.f[dP00])[k] = (D.f[dP00])[ke   ];
      (Dbuff.f[dM00])[k] = (D.f[dM00])[kw   ];
      (Dbuff.f[DIR_0P0])[k] = (D.f[DIR_0P0])[kn   ];
      (Dbuff.f[DIR_0M0])[k] = (D.f[DIR_0M0])[ks   ];
      (Dbuff.f[DIR_00P])[k] = (D.f[DIR_00P])[kt   ];
      (Dbuff.f[DIR_00M])[k] = (D.f[DIR_00M])[kb   ];
      (Dbuff.f[DIR_PP0])[k] = (D.f[DIR_PP0])[kne  ];
      (Dbuff.f[DIR_MM0])[k] = (D.f[DIR_MM0])[ksw  ];
      (Dbuff.f[DIR_PM0])[k] = (D.f[DIR_PM0])[kse  ];
      (Dbuff.f[DIR_MP0])[k] = (D.f[DIR_MP0])[knw  ];
      (Dbuff.f[DIR_P0P])[k] = (D.f[DIR_P0P])[kte  ];
      (Dbuff.f[DIR_M0M])[k] = (D.f[DIR_M0M])[kbw  ];
      (Dbuff.f[DIR_P0M])[k] = (D.f[DIR_P0M])[kbe  ];
      (Dbuff.f[DIR_M0P])[k] = (D.f[DIR_M0P])[ktw  ];
      (Dbuff.f[DIR_0PP])[k] = (D.f[DIR_0PP])[ktn  ];
      (Dbuff.f[DIR_0MM])[k] = (D.f[DIR_0MM])[kbs  ];
      (Dbuff.f[DIR_0PM])[k] = (D.f[DIR_0PM])[kbn  ];
      (Dbuff.f[DIR_0MP])[k] = (D.f[DIR_0MP])[kts  ];
      (Dbuff.f[d000])[k] = (D.f[d000])[kzero];
      (Dbuff.f[DIR_PPP])[k] = (D.f[DIR_PPP])[ktne ];
      (Dbuff.f[DIR_MMP])[k] = (D.f[DIR_MMP])[ktsw ];
      (Dbuff.f[DIR_PMP])[k] = (D.f[DIR_PMP])[ktse ];
      (Dbuff.f[DIR_MPP])[k] = (D.f[DIR_MPP])[ktnw ];
      (Dbuff.f[DIR_PPM])[k] = (D.f[DIR_PPM])[kbne ];
      (Dbuff.f[DIR_MMM])[k] = (D.f[DIR_MMM])[kbsw ];
      (Dbuff.f[DIR_PMM])[k] = (D.f[DIR_PMM])[kbse ];
      (Dbuff.f[DIR_MPM])[k] = (D.f[DIR_MPM])[kbnw ];
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
                                          unsigned long long numberOfLBnodes, 
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
         D.f[dP00] = &DD[dP00 * numberOfLBnodes];
         D.f[dM00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_PMP * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_MPM * numberOfLBnodes];
      } 
      else
      {
         D.f[dM00] = &DD[dP00 * numberOfLBnodes];
         D.f[dP00] = &DD[dM00 * numberOfLBnodes];
         D.f[DIR_0M0] = &DD[DIR_0P0 * numberOfLBnodes];
         D.f[DIR_0P0] = &DD[DIR_0M0 * numberOfLBnodes];
         D.f[DIR_00M] = &DD[DIR_00P * numberOfLBnodes];
         D.f[DIR_00P] = &DD[DIR_00M * numberOfLBnodes];
         D.f[DIR_MM0] = &DD[DIR_PP0 * numberOfLBnodes];
         D.f[DIR_PP0] = &DD[DIR_MM0 * numberOfLBnodes];
         D.f[DIR_MP0] = &DD[DIR_PM0 * numberOfLBnodes];
         D.f[DIR_PM0] = &DD[DIR_MP0 * numberOfLBnodes];
         D.f[DIR_M0M] = &DD[DIR_P0P * numberOfLBnodes];
         D.f[DIR_P0P] = &DD[DIR_M0M * numberOfLBnodes];
         D.f[DIR_M0P] = &DD[DIR_P0M * numberOfLBnodes];
         D.f[DIR_P0M] = &DD[DIR_M0P * numberOfLBnodes];
         D.f[DIR_0MM] = &DD[DIR_0PP * numberOfLBnodes];
         D.f[DIR_0PP] = &DD[DIR_0MM * numberOfLBnodes];
         D.f[DIR_0MP] = &DD[DIR_0PM * numberOfLBnodes];
         D.f[DIR_0PM] = &DD[DIR_0MP * numberOfLBnodes];
         D.f[d000] = &DD[d000 * numberOfLBnodes];
         D.f[DIR_PPP] = &DD[DIR_MMM * numberOfLBnodes];
         D.f[DIR_MMP] = &DD[DIR_PPM * numberOfLBnodes];
         D.f[DIR_PMP] = &DD[DIR_MPM * numberOfLBnodes];
         D.f[DIR_MPP] = &DD[DIR_PMM * numberOfLBnodes];
         D.f[DIR_PPM] = &DD[DIR_MMP * numberOfLBnodes];
         D.f[DIR_MMM] = &DD[DIR_PPP * numberOfLBnodes];
         D.f[DIR_PMM] = &DD[DIR_MPP * numberOfLBnodes];
         D.f[DIR_MPM] = &DD[DIR_PMP * numberOfLBnodes];
      }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //set Pointer for Buffer Fs
      Distributions27 Dbuff;
      Dbuff.f[dP00] = &bufferFs[dP00 * buffmax];
      Dbuff.f[dM00] = &bufferFs[dM00 * buffmax];
      Dbuff.f[DIR_0P0] = &bufferFs[DIR_0P0 * buffmax];
      Dbuff.f[DIR_0M0] = &bufferFs[DIR_0M0 * buffmax];
      Dbuff.f[DIR_00P] = &bufferFs[DIR_00P * buffmax];
      Dbuff.f[DIR_00M] = &bufferFs[DIR_00M * buffmax];
      Dbuff.f[DIR_PP0] = &bufferFs[DIR_PP0 * buffmax];
      Dbuff.f[DIR_MM0] = &bufferFs[DIR_MM0 * buffmax];
      Dbuff.f[DIR_PM0] = &bufferFs[DIR_PM0 * buffmax];
      Dbuff.f[DIR_MP0] = &bufferFs[DIR_MP0 * buffmax];
      Dbuff.f[DIR_P0P] = &bufferFs[DIR_P0P * buffmax];
      Dbuff.f[DIR_M0M] = &bufferFs[DIR_M0M * buffmax];
      Dbuff.f[DIR_P0M] = &bufferFs[DIR_P0M * buffmax];
      Dbuff.f[DIR_M0P] = &bufferFs[DIR_M0P * buffmax];
      Dbuff.f[DIR_0PP] = &bufferFs[DIR_0PP * buffmax];
      Dbuff.f[DIR_0MM] = &bufferFs[DIR_0MM * buffmax];
      Dbuff.f[DIR_0PM] = &bufferFs[DIR_0PM * buffmax];
      Dbuff.f[DIR_0MP] = &bufferFs[DIR_0MP * buffmax];
      Dbuff.f[d000] = &bufferFs[d000 * buffmax];
      Dbuff.f[DIR_PPP] = &bufferFs[DIR_PPP * buffmax];
      Dbuff.f[DIR_MMP] = &bufferFs[DIR_MMP * buffmax];
      Dbuff.f[DIR_PMP] = &bufferFs[DIR_PMP * buffmax];
      Dbuff.f[DIR_MPP] = &bufferFs[DIR_MPP * buffmax];
      Dbuff.f[DIR_PPM] = &bufferFs[DIR_PPM * buffmax];
      Dbuff.f[DIR_MMM] = &bufferFs[DIR_MMM * buffmax];
      Dbuff.f[DIR_PMM] = &bufferFs[DIR_PMM * buffmax];
      Dbuff.f[DIR_MPM] = &bufferFs[DIR_MPM * buffmax];
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //copy from buffer
      (D.f[dP00])[ke   ] = (Dbuff.f[dP00])[k];
      (D.f[dM00])[kw   ] = (Dbuff.f[dM00])[k];
      (D.f[DIR_0P0])[kn   ] = (Dbuff.f[DIR_0P0])[k];
      (D.f[DIR_0M0])[ks   ] = (Dbuff.f[DIR_0M0])[k];
      (D.f[DIR_00P])[kt   ] = (Dbuff.f[DIR_00P])[k];
      (D.f[DIR_00M])[kb   ] = (Dbuff.f[DIR_00M])[k];
      (D.f[DIR_PP0])[kne  ] = (Dbuff.f[DIR_PP0])[k];
      (D.f[DIR_MM0])[ksw  ] = (Dbuff.f[DIR_MM0])[k];
      (D.f[DIR_PM0])[kse  ] = (Dbuff.f[DIR_PM0])[k];
      (D.f[DIR_MP0])[knw  ] = (Dbuff.f[DIR_MP0])[k];
      (D.f[DIR_P0P])[kte  ] = (Dbuff.f[DIR_P0P])[k];
      (D.f[DIR_M0M])[kbw  ] = (Dbuff.f[DIR_M0M])[k];
      (D.f[DIR_P0M])[kbe  ] = (Dbuff.f[DIR_P0M])[k];
      (D.f[DIR_M0P])[ktw  ] = (Dbuff.f[DIR_M0P])[k];
      (D.f[DIR_0PP])[ktn  ] = (Dbuff.f[DIR_0PP])[k];
      (D.f[DIR_0MM])[kbs  ] = (Dbuff.f[DIR_0MM])[k];
      (D.f[DIR_0PM])[kbn  ] = (Dbuff.f[DIR_0PM])[k];
      (D.f[DIR_0MP])[kts  ] = (Dbuff.f[DIR_0MP])[k];
      (D.f[d000])[kzero] = (Dbuff.f[d000])[k];
      (D.f[DIR_PPP])[ktne ] = (Dbuff.f[DIR_PPP])[k];
      (D.f[DIR_MMP])[ktsw ] = (Dbuff.f[DIR_MMP])[k];
      (D.f[DIR_PMP])[ktse ] = (Dbuff.f[DIR_PMP])[k];
      (D.f[DIR_MPP])[ktnw ] = (Dbuff.f[DIR_MPP])[k];
      (D.f[DIR_PPM])[kbne ] = (Dbuff.f[DIR_PPM])[k];
      (D.f[DIR_MMM])[kbsw ] = (Dbuff.f[DIR_MMM])[k];
      (D.f[DIR_PMM])[kbse ] = (Dbuff.f[DIR_PMM])[k];
      (D.f[DIR_MPM])[kbnw ] = (Dbuff.f[DIR_MPM])[k];
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
	unsigned long long numberOfLBnodes,
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
			G.g[dP00] = &G6[dP00 * numberOfLBnodes];
			G.g[dM00] = &G6[dM00 * numberOfLBnodes];
			G.g[DIR_0P0] = &G6[DIR_0P0 * numberOfLBnodes];
			G.g[DIR_0M0] = &G6[DIR_0M0 * numberOfLBnodes];
			G.g[DIR_00P] = &G6[DIR_00P * numberOfLBnodes];
			G.g[DIR_00M] = &G6[DIR_00M * numberOfLBnodes];
		}
		else
		{
			G.g[dM00] = &G6[dP00 * numberOfLBnodes];
			G.g[dP00] = &G6[dM00 * numberOfLBnodes];
			G.g[DIR_0M0] = &G6[DIR_0P0 * numberOfLBnodes];
			G.g[DIR_0P0] = &G6[DIR_0M0 * numberOfLBnodes];
			G.g[DIR_00M] = &G6[DIR_00P * numberOfLBnodes];
			G.g[DIR_00P] = &G6[DIR_00M * numberOfLBnodes];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Buffer Gs
		Distributions6 Dbuff;
		Dbuff.g[dP00] = &bufferGs[dP00 * buffmax];
		Dbuff.g[dM00] = &bufferGs[dM00 * buffmax];
		Dbuff.g[DIR_0P0] = &bufferGs[DIR_0P0 * buffmax];
		Dbuff.g[DIR_0M0] = &bufferGs[DIR_0M0 * buffmax];
		Dbuff.g[DIR_00P] = &bufferGs[DIR_00P * buffmax];
		Dbuff.g[DIR_00M] = &bufferGs[DIR_00M * buffmax];
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//write Gs to buffer
		(Dbuff.g[dP00])[k] = (G.g[dM00])[kw];
		(Dbuff.g[dM00])[k] = (G.g[dP00])[kr];
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
	unsigned long long numberOfLBnodes,
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
			G.g[dP00] = &G6[dP00 * numberOfLBnodes];
			G.g[dM00] = &G6[dM00 * numberOfLBnodes];
			G.g[DIR_0P0] = &G6[DIR_0P0 * numberOfLBnodes];
			G.g[DIR_0M0] = &G6[DIR_0M0 * numberOfLBnodes];
			G.g[DIR_00P] = &G6[DIR_00P * numberOfLBnodes];
			G.g[DIR_00M] = &G6[DIR_00M * numberOfLBnodes];
		}
		else
		{
			G.g[dM00] = &G6[dP00 * numberOfLBnodes];
			G.g[dP00] = &G6[dM00 * numberOfLBnodes];
			G.g[DIR_0M0] = &G6[DIR_0P0 * numberOfLBnodes];
			G.g[DIR_0P0] = &G6[DIR_0M0 * numberOfLBnodes];
			G.g[DIR_00M] = &G6[DIR_00P * numberOfLBnodes];
			G.g[DIR_00P] = &G6[DIR_00M * numberOfLBnodes];
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//set Pointer for Buffer Gs
		Distributions6 Dbuff;
		Dbuff.g[dP00] = &bufferGs[dP00 * buffmax];
		Dbuff.g[dM00] = &bufferGs[dM00 * buffmax];
		Dbuff.g[DIR_0P0] = &bufferGs[DIR_0P0 * buffmax];
		Dbuff.g[DIR_0M0] = &bufferGs[DIR_0M0 * buffmax];
		Dbuff.g[DIR_00P] = &bufferGs[DIR_00P * buffmax];
		Dbuff.g[DIR_00M] = &bufferGs[DIR_00M * buffmax];
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//write buffer to Gs
		(G.g[dM00])[kw] = (Dbuff.g[dP00])[k];
		(G.g[dP00])[kr] = (Dbuff.g[dM00])[k];
		(G.g[DIR_0M0])[ks] = (Dbuff.g[DIR_0P0])[k];
		(G.g[DIR_0P0])[kr] = (Dbuff.g[DIR_0M0])[k];
		(G.g[DIR_00M])[kb] = (Dbuff.g[DIR_00P])[k];
		(G.g[DIR_00P])[kr] = (Dbuff.g[DIR_00M])[k];
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
