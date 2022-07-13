/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc2ndMomentsIncompSP27(  real* kxyFromfcNEQ,
														real* kyzFromfcNEQ,
														real* kxzFromfcNEQ,
														real* kxxMyyFromfcNEQ,
														real* kxxMzzFromfcNEQ,
														unsigned int* geoD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														unsigned int size_Mat,
														real* DD,
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

   if(k < size_Mat)
   {
      //////////////////////////////////////////////////////////////////////////
      //index
      //unsigned int kzero= k;
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
      real        f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,/*f_ZERO,*/f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;
	  f_E    = (D.f[E   ])[ke   ];
	  f_W    = (D.f[W   ])[kw   ];
	  f_N    = (D.f[N   ])[kn   ];
	  f_S    = (D.f[S   ])[ks   ];
	  f_T    = (D.f[T   ])[kt   ];
	  f_B    = (D.f[B   ])[kb   ];
	  f_NE   = (D.f[NE  ])[kne  ];
	  f_SW   = (D.f[SW  ])[ksw  ];
	  f_SE   = (D.f[SE  ])[kse  ];
	  f_NW   = (D.f[NW  ])[knw  ];
	  f_TE   = (D.f[TE  ])[kte  ];
	  f_BW   = (D.f[BW  ])[kbw  ];
	  f_BE   = (D.f[BE  ])[kbe  ];
	  f_TW   = (D.f[TW  ])[ktw  ];
	  f_TN   = (D.f[TN  ])[ktn  ];
	  f_BS   = (D.f[BS  ])[kbs  ];
	  f_BN   = (D.f[BN  ])[kbn  ];
	  f_TS   = (D.f[TS  ])[kts  ];
	  //f_ZERO = (D.f[REST])[kzero];
	  f_TNE  = (D.f[TNE ])[ktne ];
	  f_TSW  = (D.f[TSW ])[ktsw ];
	  f_TSE  = (D.f[TSE ])[ktse ];
	  f_TNW  = (D.f[TNW ])[ktnw ];
	  f_BNE  = (D.f[BNE ])[kbne ];
	  f_BSW  = (D.f[BSW ])[kbsw ];
	  f_BSE  = (D.f[BSE ])[kbse ];
	  f_BNW  = (D.f[BNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
	  real vx1, vx2, vx3;
      kxyFromfcNEQ[k]       = c0o1;
	  kyzFromfcNEQ[k]       = c0o1;
	  kxzFromfcNEQ[k]       = c0o1;
	  kxxMyyFromfcNEQ[k]    = c0o1;
	  kxxMzzFromfcNEQ[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
		  vx1                = ((f_TNE-f_BSW)+(f_BSE-f_TNW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)) + (((f_NE-f_SW)+(f_TE-f_BW))+((f_SE-f_NW)+(f_BE-f_TW))) + (f_E-f_W);
		  vx2                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_BNE-f_TSW)+(f_BNW-f_TSE)) + (((f_NE-f_SW)+(f_TN-f_BS))+((f_BN-f_TS)+(f_NW-f_SE))) + (f_N-f_S);
		  vx3                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_TSW-f_BNE)+(f_TSE-f_BNW)) + (((f_TE-f_BW)+(f_TN-f_BS))+((f_TW-f_BE)+(f_TS-f_BN))) + (f_T-f_B);
		  kxyFromfcNEQ[k]    = -c3o1 *(f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE-(vx1*vx2));
		  kyzFromfcNEQ[k]    = -c3o1 *(f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW-(vx2*vx3));
		  kxzFromfcNEQ[k]    = -c3o1 *(f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE-(vx1*vx3));
		  kxxMyyFromfcNEQ[k] = -c3o2 * (f_BW+f_W+f_TW-f_BS-f_S-f_TS-f_BN-f_N-f_TN+f_BE+f_E+f_TE-(vx1*vx1-vx2*vx2));		//all E+W minus all N+S (no combinations of xy left)
		  kxxMzzFromfcNEQ[k] = -c3o2 * (f_SW+f_W+f_NW-f_BS-f_TS-f_B-f_T-f_BN-f_TN+f_SE+f_E+f_NE-(vx1*vx1-vx3*vx3));		//all E+W minus all T+B (no combinations of xz left)
      }
   }
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc2ndMomentsCompSP27(real* kxyFromfcNEQ,
													real* kyzFromfcNEQ,
													real* kxzFromfcNEQ,
													real* kxxMyyFromfcNEQ,
													real* kxxMzzFromfcNEQ,
													unsigned int* geoD,
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													unsigned int size_Mat,
													real* DD,
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

   if(k < size_Mat)
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
      real f_ZERO;
      real        f_E,f_W,f_N,f_S,f_T,f_B,f_NE,f_SW,f_SE,f_NW,f_TE,f_BW,f_BE,f_TW,f_TN,f_BS,f_BN,f_TS,f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;
	  f_E    = (D.f[E   ])[ke   ];
	  f_W    = (D.f[W   ])[kw   ];
	  f_N    = (D.f[N   ])[kn   ];
	  f_S    = (D.f[S   ])[ks   ];
	  f_T    = (D.f[T   ])[kt   ];
	  f_B    = (D.f[B   ])[kb   ];
	  f_NE   = (D.f[NE  ])[kne  ];
	  f_SW   = (D.f[SW  ])[ksw  ];
	  f_SE   = (D.f[SE  ])[kse  ];
	  f_NW   = (D.f[NW  ])[knw  ];
	  f_TE   = (D.f[TE  ])[kte  ];
	  f_BW   = (D.f[BW  ])[kbw  ];
	  f_BE   = (D.f[BE  ])[kbe  ];
	  f_TW   = (D.f[TW  ])[ktw  ];
	  f_TN   = (D.f[TN  ])[ktn  ];
	  f_BS   = (D.f[BS  ])[kbs  ];
	  f_BN   = (D.f[BN  ])[kbn  ];
	  f_TS   = (D.f[TS  ])[kts  ];
	  f_ZERO = (D.f[REST])[kzero];
	  f_TNE  = (D.f[TNE ])[ktne ];
	  f_TSW  = (D.f[TSW ])[ktsw ];
	  f_TSE  = (D.f[TSE ])[ktse ];
	  f_TNW  = (D.f[TNW ])[ktnw ];
	  f_BNE  = (D.f[BNE ])[kbne ];
	  f_BSW  = (D.f[BSW ])[kbsw ];
	  f_BSE  = (D.f[BSE ])[kbse ];
	  f_BNW  = (D.f[BNW ])[kbnw ];
      //////////////////////////////////////////////////////////////////////////
	  real drho;
	  real vx1, vx2, vx3, rho;
      kxyFromfcNEQ[k]       = c0o1;
	  kyzFromfcNEQ[k]       = c0o1;
	  kxzFromfcNEQ[k]       = c0o1;
	  kxxMyyFromfcNEQ[k]    = c0o1;
	  kxxMzzFromfcNEQ[k]    = c0o1;

      if(geoD[k] == GEO_FLUID)
      {
          drho               = ((f_TNE+f_BSW)+(f_BSE+f_TNW)+(f_BNE+f_TSW)+(f_TSE+f_BNW)) +
		 					   ((f_NE+f_SW)+(f_TE+f_BW)+(f_SE+f_NW)+(f_BE+f_TW)+(f_BN+f_TS)+(f_TN+f_BS)) +
		 					   ((f_E-f_W) + (f_N-f_S) + (f_T-f_B)) + f_ZERO;
		  rho                = drho + c1o1;
		  vx1                = ((f_TNE-f_BSW)+(f_BSE-f_TNW)+(f_BNE-f_TSW)+(f_TSE-f_BNW)) + (((f_NE-f_SW)+(f_TE-f_BW))+((f_SE-f_NW)+(f_BE-f_TW))) + (f_E-f_W) / rho;
		  vx2                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_BNE-f_TSW)+(f_BNW-f_TSE)) + (((f_NE-f_SW)+(f_TN-f_BS))+((f_BN-f_TS)+(f_NW-f_SE))) + (f_N-f_S) / rho;
		  vx3                = ((f_TNE-f_BSW)+(f_TNW-f_BSE)+(f_TSW-f_BNE)+(f_TSE-f_BNW)) + (((f_TE-f_BW)+(f_TN-f_BS))+((f_TW-f_BE)+(f_TS-f_BN))) + (f_T-f_B) / rho;
		  kxyFromfcNEQ[k]    = -c3o1 *(f_SW+f_BSW+f_TSW-f_NW-f_BNW-f_TNW-f_SE-f_BSE-f_TSE+f_NE+f_BNE+f_TNE-(vx1*vx2));
		  kyzFromfcNEQ[k]    = -c3o1 *(f_BS+f_BSE+f_BSW-f_TS-f_TSE-f_TSW-f_BN-f_BNE-f_BNW+f_TN+f_TNE+f_TNW-(vx2*vx3));
		  kxzFromfcNEQ[k]    = -c3o1 *(f_BW+f_BSW+f_BNW-f_TW-f_TSW-f_TNW-f_BE-f_BSE-f_BNE+f_TE+f_TSE+f_TNE-(vx1*vx3));
		  kxxMyyFromfcNEQ[k] = -c3o2 * (f_BW+f_W+f_TW-f_BS-f_S-f_TS-f_BN-f_N-f_TN+f_BE+f_E+f_TE-(vx1*vx1-vx2*vx2));		//all E+W minus all N+S (no combinations of xy left)
		  kxxMzzFromfcNEQ[k] = -c3o2 * (f_SW+f_W+f_NW-f_BS-f_TS-f_B-f_T-f_BN-f_TN+f_SE+f_E+f_NE-(vx1*vx1-vx3*vx3));		//all E+W minus all T+B (no combinations of xz left)
      }
   }
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc3rdMomentsIncompSP27(  real* CUMbbb,
														real* CUMabc,
														real* CUMbac,
														real* CUMbca,
														real* CUMcba,
														real* CUMacb,
														real* CUMcab,
														unsigned int* bcMatD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* DDStart,
														int size_Mat,
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
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
			{
				D.f[E   ] = &DDStart[E   *size_Mat];
				D.f[W   ] = &DDStart[W   *size_Mat];
				D.f[N   ] = &DDStart[N   *size_Mat];
				D.f[S   ] = &DDStart[S   *size_Mat];
				D.f[T   ] = &DDStart[T   *size_Mat];
				D.f[B   ] = &DDStart[B   *size_Mat];
				D.f[NE  ] = &DDStart[NE  *size_Mat];
				D.f[SW  ] = &DDStart[SW  *size_Mat];
				D.f[SE  ] = &DDStart[SE  *size_Mat];
				D.f[NW  ] = &DDStart[NW  *size_Mat];
				D.f[TE  ] = &DDStart[TE  *size_Mat];
				D.f[BW  ] = &DDStart[BW  *size_Mat];
				D.f[BE  ] = &DDStart[BE  *size_Mat];
				D.f[TW  ] = &DDStart[TW  *size_Mat];
				D.f[TN  ] = &DDStart[TN  *size_Mat];
				D.f[BS  ] = &DDStart[BS  *size_Mat];
				D.f[BN  ] = &DDStart[BN  *size_Mat];
				D.f[TS  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[TNE ] = &DDStart[TNE *size_Mat];
				D.f[TSW ] = &DDStart[TSW *size_Mat];
				D.f[TSE ] = &DDStart[TSE *size_Mat];
				D.f[TNW ] = &DDStart[TNW *size_Mat];
				D.f[BNE ] = &DDStart[BNE *size_Mat];
				D.f[BSW ] = &DDStart[BSW *size_Mat];
				D.f[BSE ] = &DDStart[BSE *size_Mat];
				D.f[BNW ] = &DDStart[BNW *size_Mat];
			}
			else
			{
				D.f[W   ] = &DDStart[E   *size_Mat];
				D.f[E   ] = &DDStart[W   *size_Mat];
				D.f[S   ] = &DDStart[N   *size_Mat];
				D.f[N   ] = &DDStart[S   *size_Mat];
				D.f[B   ] = &DDStart[T   *size_Mat];
				D.f[T   ] = &DDStart[B   *size_Mat];
				D.f[SW  ] = &DDStart[NE  *size_Mat];
				D.f[NE  ] = &DDStart[SW  *size_Mat];
				D.f[NW  ] = &DDStart[SE  *size_Mat];
				D.f[SE  ] = &DDStart[NW  *size_Mat];
				D.f[BW  ] = &DDStart[TE  *size_Mat];
				D.f[TE  ] = &DDStart[BW  *size_Mat];
				D.f[TW  ] = &DDStart[BE  *size_Mat];
				D.f[BE  ] = &DDStart[TW  *size_Mat];
				D.f[BS  ] = &DDStart[TN  *size_Mat];
				D.f[TN  ] = &DDStart[BS  *size_Mat];
				D.f[TS  ] = &DDStart[BN  *size_Mat];
				D.f[BN  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[BSW ] = &DDStart[TNE *size_Mat];
				D.f[BNE ] = &DDStart[TSW *size_Mat];
				D.f[BNW ] = &DDStart[TSE *size_Mat];
				D.f[BSE ] = &DDStart[TNW *size_Mat];
				D.f[TSW ] = &DDStart[BNE *size_Mat];
				D.f[TNE ] = &DDStart[BSW *size_Mat];
				D.f[TNW ] = &DDStart[BSE *size_Mat];
				D.f[TSE ] = &DDStart[BNW *size_Mat];
			}

			////////////////////////////////////////////////////////////////////////////////
			//index
			unsigned int kw   = neighborX[k];
			unsigned int ks   = neighborY[k];
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			unsigned int kbw  = neighborZ[kw];
			unsigned int kbs  = neighborZ[ks];
			unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[E   ])[k  ];
			real mfabb = (D.f[W   ])[kw ];
			real mfbcb = (D.f[N   ])[k  ];
			real mfbab = (D.f[S   ])[ks ];
			real mfbbc = (D.f[T   ])[k  ];
			real mfbba = (D.f[B   ])[kb ];
			real mfccb = (D.f[NE  ])[k  ];
			real mfaab = (D.f[SW  ])[ksw];
			real mfcab = (D.f[SE  ])[ks ];
			real mfacb = (D.f[NW  ])[kw ];
			real mfcbc = (D.f[TE  ])[k  ];
			real mfaba = (D.f[BW  ])[kbw];
			real mfcba = (D.f[BE  ])[kb ];
			real mfabc = (D.f[TW  ])[kw ];
			real mfbcc = (D.f[TN  ])[k  ];
			real mfbaa = (D.f[BS  ])[kbs];
			real mfbca = (D.f[BN  ])[kb ];
			real mfbac = (D.f[TS  ])[ks ];
			real mfbbb = (D.f[REST])[k  ];
			real mfccc = (D.f[TNE ])[k  ];
			real mfaac = (D.f[TSW ])[ksw];
			real mfcac = (D.f[TSE ])[ks ];
			real mfacc = (D.f[TNW ])[kw ];
			real mfcca = (D.f[BNE ])[kb ];
			real mfaaa = (D.f[BSW ])[kbsw];
			real mfcaa = (D.f[BSE ])[kbs];
			real mfaca = (D.f[BNW ])[kbw];
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb));
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab));
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba));
			////////////////////////////////////////////////////////////////////////////////////
			real oMdrho = c1o1 - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
								   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
								   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);
			////////////////////////////////////////////////////////////////////////////////////
			real m0, m1, m2;	
			real vx2;
			real vy2;
			real vz2;
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			//3.
			CUMbbb[k] = mfbbb;
			CUMabc[k] = mfabc;
			CUMbac[k] = mfbac;
			CUMbca[k] = mfbca;
			CUMcba[k] = mfcba;
			CUMacb[k] = mfacb;
			CUMcab[k] = mfcab;
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalc3rdMomentsCompSP27(real* CUMbbb,
													real* CUMabc,
													real* CUMbac,
													real* CUMbca,
													real* CUMcba,
													real* CUMacb,
													real* CUMcab,
													unsigned int* bcMatD,
													unsigned int* neighborX,
													unsigned int* neighborY,
													unsigned int* neighborZ,
													real* DDStart,
													int size_Mat,
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
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
			{
				D.f[E   ] = &DDStart[E   *size_Mat];
				D.f[W   ] = &DDStart[W   *size_Mat];
				D.f[N   ] = &DDStart[N   *size_Mat];
				D.f[S   ] = &DDStart[S   *size_Mat];
				D.f[T   ] = &DDStart[T   *size_Mat];
				D.f[B   ] = &DDStart[B   *size_Mat];
				D.f[NE  ] = &DDStart[NE  *size_Mat];
				D.f[SW  ] = &DDStart[SW  *size_Mat];
				D.f[SE  ] = &DDStart[SE  *size_Mat];
				D.f[NW  ] = &DDStart[NW  *size_Mat];
				D.f[TE  ] = &DDStart[TE  *size_Mat];
				D.f[BW  ] = &DDStart[BW  *size_Mat];
				D.f[BE  ] = &DDStart[BE  *size_Mat];
				D.f[TW  ] = &DDStart[TW  *size_Mat];
				D.f[TN  ] = &DDStart[TN  *size_Mat];
				D.f[BS  ] = &DDStart[BS  *size_Mat];
				D.f[BN  ] = &DDStart[BN  *size_Mat];
				D.f[TS  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[TNE ] = &DDStart[TNE *size_Mat];
				D.f[TSW ] = &DDStart[TSW *size_Mat];
				D.f[TSE ] = &DDStart[TSE *size_Mat];
				D.f[TNW ] = &DDStart[TNW *size_Mat];
				D.f[BNE ] = &DDStart[BNE *size_Mat];
				D.f[BSW ] = &DDStart[BSW *size_Mat];
				D.f[BSE ] = &DDStart[BSE *size_Mat];
				D.f[BNW ] = &DDStart[BNW *size_Mat];
			}
			else
			{
				D.f[W   ] = &DDStart[E   *size_Mat];
				D.f[E   ] = &DDStart[W   *size_Mat];
				D.f[S   ] = &DDStart[N   *size_Mat];
				D.f[N   ] = &DDStart[S   *size_Mat];
				D.f[B   ] = &DDStart[T   *size_Mat];
				D.f[T   ] = &DDStart[B   *size_Mat];
				D.f[SW  ] = &DDStart[NE  *size_Mat];
				D.f[NE  ] = &DDStart[SW  *size_Mat];
				D.f[NW  ] = &DDStart[SE  *size_Mat];
				D.f[SE  ] = &DDStart[NW  *size_Mat];
				D.f[BW  ] = &DDStart[TE  *size_Mat];
				D.f[TE  ] = &DDStart[BW  *size_Mat];
				D.f[TW  ] = &DDStart[BE  *size_Mat];
				D.f[BE  ] = &DDStart[TW  *size_Mat];
				D.f[BS  ] = &DDStart[TN  *size_Mat];
				D.f[TN  ] = &DDStart[BS  *size_Mat];
				D.f[TS  ] = &DDStart[BN  *size_Mat];
				D.f[BN  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[BSW ] = &DDStart[TNE *size_Mat];
				D.f[BNE ] = &DDStart[TSW *size_Mat];
				D.f[BNW ] = &DDStart[TSE *size_Mat];
				D.f[BSE ] = &DDStart[TNW *size_Mat];
				D.f[TSW ] = &DDStart[BNE *size_Mat];
				D.f[TNE ] = &DDStart[BSW *size_Mat];
				D.f[TNW ] = &DDStart[BSE *size_Mat];
				D.f[TSE ] = &DDStart[BNW *size_Mat];
			}

			////////////////////////////////////////////////////////////////////////////////
			//index
			unsigned int kw   = neighborX[k];
			unsigned int ks   = neighborY[k];
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			unsigned int kbw  = neighborZ[kw];
			unsigned int kbs  = neighborZ[ks];
			unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[E   ])[k  ];
			real mfabb = (D.f[W   ])[kw ];
			real mfbcb = (D.f[N   ])[k  ];
			real mfbab = (D.f[S   ])[ks ];
			real mfbbc = (D.f[T   ])[k  ];
			real mfbba = (D.f[B   ])[kb ];
			real mfccb = (D.f[NE  ])[k  ];
			real mfaab = (D.f[SW  ])[ksw];
			real mfcab = (D.f[SE  ])[ks ];
			real mfacb = (D.f[NW  ])[kw ];
			real mfcbc = (D.f[TE  ])[k  ];
			real mfaba = (D.f[BW  ])[kbw];
			real mfcba = (D.f[BE  ])[kb ];
			real mfabc = (D.f[TW  ])[kw ];
			real mfbcc = (D.f[TN  ])[k  ];
			real mfbaa = (D.f[BS  ])[kbs];
			real mfbca = (D.f[BN  ])[kb ];
			real mfbac = (D.f[TS  ])[ks ];
			real mfbbb = (D.f[REST])[k  ];
			real mfccc = (D.f[TNE ])[k  ];
			real mfaac = (D.f[TSW ])[ksw];
			real mfcac = (D.f[TSE ])[ks ];
			real mfacc = (D.f[TNW ])[kw ];
			real mfcca = (D.f[BNE ])[kb ];
			real mfaaa = (D.f[BSW ])[kbsw];
			real mfcaa = (D.f[BSE ])[kbs];
			real mfaca = (D.f[BNW ])[kbw];
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb)) + (mfbba+mfbbc)) + mfbbb;

			real rho = c1o1+drho;
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			real oMdrho = c1o1; // comp special
			////////////////////////////////////////////////////////////////////////////////////
			real m0, m1, m2;	
			real vx2;
			real vy2;
			real vz2;
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			//3.
			CUMbbb[k] = mfbbb;
			CUMabc[k] = mfabc;
			CUMbac[k] = mfbac;
			CUMbca[k] = mfbca;
			CUMcba[k] = mfcba;
			CUMacb[k] = mfacb;
			CUMcab[k] = mfcab;
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcHigherMomentsIncompSP27(   real* CUMcbb,
															real* CUMbcb,
															real* CUMbbc,
															real* CUMcca,
															real* CUMcac,
															real* CUMacc,
															real* CUMbcc,
															real* CUMcbc,
															real* CUMccb,
															real* CUMccc,
															unsigned int* bcMatD,
															unsigned int* neighborX,
															unsigned int* neighborY,
															unsigned int* neighborZ,
															real* DDStart,
															int size_Mat,
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
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
			{
				D.f[E   ] = &DDStart[E   *size_Mat];
				D.f[W   ] = &DDStart[W   *size_Mat];
				D.f[N   ] = &DDStart[N   *size_Mat];
				D.f[S   ] = &DDStart[S   *size_Mat];
				D.f[T   ] = &DDStart[T   *size_Mat];
				D.f[B   ] = &DDStart[B   *size_Mat];
				D.f[NE  ] = &DDStart[NE  *size_Mat];
				D.f[SW  ] = &DDStart[SW  *size_Mat];
				D.f[SE  ] = &DDStart[SE  *size_Mat];
				D.f[NW  ] = &DDStart[NW  *size_Mat];
				D.f[TE  ] = &DDStart[TE  *size_Mat];
				D.f[BW  ] = &DDStart[BW  *size_Mat];
				D.f[BE  ] = &DDStart[BE  *size_Mat];
				D.f[TW  ] = &DDStart[TW  *size_Mat];
				D.f[TN  ] = &DDStart[TN  *size_Mat];
				D.f[BS  ] = &DDStart[BS  *size_Mat];
				D.f[BN  ] = &DDStart[BN  *size_Mat];
				D.f[TS  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[TNE ] = &DDStart[TNE *size_Mat];
				D.f[TSW ] = &DDStart[TSW *size_Mat];
				D.f[TSE ] = &DDStart[TSE *size_Mat];
				D.f[TNW ] = &DDStart[TNW *size_Mat];
				D.f[BNE ] = &DDStart[BNE *size_Mat];
				D.f[BSW ] = &DDStart[BSW *size_Mat];
				D.f[BSE ] = &DDStart[BSE *size_Mat];
				D.f[BNW ] = &DDStart[BNW *size_Mat];
			}
			else
			{
				D.f[W   ] = &DDStart[E   *size_Mat];
				D.f[E   ] = &DDStart[W   *size_Mat];
				D.f[S   ] = &DDStart[N   *size_Mat];
				D.f[N   ] = &DDStart[S   *size_Mat];
				D.f[B   ] = &DDStart[T   *size_Mat];
				D.f[T   ] = &DDStart[B   *size_Mat];
				D.f[SW  ] = &DDStart[NE  *size_Mat];
				D.f[NE  ] = &DDStart[SW  *size_Mat];
				D.f[NW  ] = &DDStart[SE  *size_Mat];
				D.f[SE  ] = &DDStart[NW  *size_Mat];
				D.f[BW  ] = &DDStart[TE  *size_Mat];
				D.f[TE  ] = &DDStart[BW  *size_Mat];
				D.f[TW  ] = &DDStart[BE  *size_Mat];
				D.f[BE  ] = &DDStart[TW  *size_Mat];
				D.f[BS  ] = &DDStart[TN  *size_Mat];
				D.f[TN  ] = &DDStart[BS  *size_Mat];
				D.f[TS  ] = &DDStart[BN  *size_Mat];
				D.f[BN  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[BSW ] = &DDStart[TNE *size_Mat];
				D.f[BNE ] = &DDStart[TSW *size_Mat];
				D.f[BNW ] = &DDStart[TSE *size_Mat];
				D.f[BSE ] = &DDStart[TNW *size_Mat];
				D.f[TSW ] = &DDStart[BNE *size_Mat];
				D.f[TNE ] = &DDStart[BSW *size_Mat];
				D.f[TNW ] = &DDStart[BSE *size_Mat];
				D.f[TSE ] = &DDStart[BNW *size_Mat];
			}

			////////////////////////////////////////////////////////////////////////////////
			//index
			unsigned int kw   = neighborX[k];
			unsigned int ks   = neighborY[k];
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			unsigned int kbw  = neighborZ[kw];
			unsigned int kbs  = neighborZ[ks];
			unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[E   ])[k  ];
			real mfabb = (D.f[W   ])[kw ];
			real mfbcb = (D.f[N   ])[k  ];
			real mfbab = (D.f[S   ])[ks ];
			real mfbbc = (D.f[T   ])[k  ];
			real mfbba = (D.f[B   ])[kb ];
			real mfccb = (D.f[NE  ])[k  ];
			real mfaab = (D.f[SW  ])[ksw];
			real mfcab = (D.f[SE  ])[ks ];
			real mfacb = (D.f[NW  ])[kw ];
			real mfcbc = (D.f[TE  ])[k  ];
			real mfaba = (D.f[BW  ])[kbw];
			real mfcba = (D.f[BE  ])[kb ];
			real mfabc = (D.f[TW  ])[kw ];
			real mfbcc = (D.f[TN  ])[k  ];
			real mfbaa = (D.f[BS  ])[kbs];
			real mfbca = (D.f[BN  ])[kb ];
			real mfbac = (D.f[TS  ])[ks ];
			real mfbbb = (D.f[REST])[k  ];
			real mfccc = (D.f[TNE ])[k  ];
			real mfaac = (D.f[TSW ])[ksw];
			real mfcac = (D.f[TSE ])[ks ];
			real mfacc = (D.f[TNW ])[kw ];
			real mfcca = (D.f[BNE ])[kb ];
			real mfaaa = (D.f[BSW ])[kbsw];
			real mfcaa = (D.f[BSE ])[kbs];
			real mfaca = (D.f[BNW ])[kbw];
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb));
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab));
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba));
			////////////////////////////////////////////////////////////////////////////////////
			real oMdrho = c1o1 - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
								   mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
								   mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);
			////////////////////////////////////////////////////////////////////////////////////
			real m0, m1, m2;	
			real vx2;
			real vy2;
			real vz2;
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			//Cum 4.
			CUMcbb[k]      = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + c2o1 * mfbba * mfbab); 
			CUMbcb[k]      = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + c2o1 * mfbba * mfabb); 
			CUMbbc[k]      = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + c2o1 * mfbab * mfabb); 

			CUMcca[k]      = mfcca - ((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
			CUMcac[k]      = mfcac - ((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);
			CUMacc[k]      = mfacc - ((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-c1o1)*oMdrho);

			//Cum 5.
			CUMbcc[k]      = mfbcc - (mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) - c1o3 * (mfbca + mfbac) * oMdrho;
			CUMcbc[k]      = mfcbc - (mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) - c1o3 * (mfcba + mfabc) * oMdrho;
			CUMccb[k]      = mfccb - (mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) - c1o3 * (mfacb + mfcab) * oMdrho;

			//Cum 6.
			CUMccc[k]      = mfccc  +((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb))
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb)
							-    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
							-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(c1o1-c2o1* oMdrho)- c1o27* oMdrho * oMdrho*(-c2o1* oMdrho)
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) * c2o3*oMdrho) +c1o27*oMdrho;
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}































////////////////////////////////////////////////////////////////////////////////
__global__ void LBCalcHigherMomentsCompSP27( real* CUMcbb,
														real* CUMbcb,
														real* CUMbbc,
														real* CUMcca,
														real* CUMcac,
														real* CUMacc,
														real* CUMbcc,
														real* CUMcbc,
														real* CUMccb,
														real* CUMccc,
														unsigned int* bcMatD,
														unsigned int* neighborX,
														unsigned int* neighborY,
														unsigned int* neighborZ,
														real* DDStart,
														int size_Mat,
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
		BC = bcMatD[k];

		if( (BC != GEO_SOLID) && (BC != GEO_VOID) )
		{
			Distributions27 D;
			if (EvenOrOdd==true)
			{
				D.f[E   ] = &DDStart[E   *size_Mat];
				D.f[W   ] = &DDStart[W   *size_Mat];
				D.f[N   ] = &DDStart[N   *size_Mat];
				D.f[S   ] = &DDStart[S   *size_Mat];
				D.f[T   ] = &DDStart[T   *size_Mat];
				D.f[B   ] = &DDStart[B   *size_Mat];
				D.f[NE  ] = &DDStart[NE  *size_Mat];
				D.f[SW  ] = &DDStart[SW  *size_Mat];
				D.f[SE  ] = &DDStart[SE  *size_Mat];
				D.f[NW  ] = &DDStart[NW  *size_Mat];
				D.f[TE  ] = &DDStart[TE  *size_Mat];
				D.f[BW  ] = &DDStart[BW  *size_Mat];
				D.f[BE  ] = &DDStart[BE  *size_Mat];
				D.f[TW  ] = &DDStart[TW  *size_Mat];
				D.f[TN  ] = &DDStart[TN  *size_Mat];
				D.f[BS  ] = &DDStart[BS  *size_Mat];
				D.f[BN  ] = &DDStart[BN  *size_Mat];
				D.f[TS  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[TNE ] = &DDStart[TNE *size_Mat];
				D.f[TSW ] = &DDStart[TSW *size_Mat];
				D.f[TSE ] = &DDStart[TSE *size_Mat];
				D.f[TNW ] = &DDStart[TNW *size_Mat];
				D.f[BNE ] = &DDStart[BNE *size_Mat];
				D.f[BSW ] = &DDStart[BSW *size_Mat];
				D.f[BSE ] = &DDStart[BSE *size_Mat];
				D.f[BNW ] = &DDStart[BNW *size_Mat];
			}
			else
			{
				D.f[W   ] = &DDStart[E   *size_Mat];
				D.f[E   ] = &DDStart[W   *size_Mat];
				D.f[S   ] = &DDStart[N   *size_Mat];
				D.f[N   ] = &DDStart[S   *size_Mat];
				D.f[B   ] = &DDStart[T   *size_Mat];
				D.f[T   ] = &DDStart[B   *size_Mat];
				D.f[SW  ] = &DDStart[NE  *size_Mat];
				D.f[NE  ] = &DDStart[SW  *size_Mat];
				D.f[NW  ] = &DDStart[SE  *size_Mat];
				D.f[SE  ] = &DDStart[NW  *size_Mat];
				D.f[BW  ] = &DDStart[TE  *size_Mat];
				D.f[TE  ] = &DDStart[BW  *size_Mat];
				D.f[TW  ] = &DDStart[BE  *size_Mat];
				D.f[BE  ] = &DDStart[TW  *size_Mat];
				D.f[BS  ] = &DDStart[TN  *size_Mat];
				D.f[TN  ] = &DDStart[BS  *size_Mat];
				D.f[TS  ] = &DDStart[BN  *size_Mat];
				D.f[BN  ] = &DDStart[TS  *size_Mat];
				D.f[REST] = &DDStart[REST*size_Mat];
				D.f[BSW ] = &DDStart[TNE *size_Mat];
				D.f[BNE ] = &DDStart[TSW *size_Mat];
				D.f[BNW ] = &DDStart[TSE *size_Mat];
				D.f[BSE ] = &DDStart[TNW *size_Mat];
				D.f[TSW ] = &DDStart[BNE *size_Mat];
				D.f[TNE ] = &DDStart[BSW *size_Mat];
				D.f[TNW ] = &DDStart[BSE *size_Mat];
				D.f[TSE ] = &DDStart[BNW *size_Mat];
			}

			////////////////////////////////////////////////////////////////////////////////
			//index
			unsigned int kw   = neighborX[k];
			unsigned int ks   = neighborY[k];
			unsigned int kb   = neighborZ[k];
			unsigned int ksw  = neighborY[kw];
			unsigned int kbw  = neighborZ[kw];
			unsigned int kbs  = neighborZ[ks];
			unsigned int kbsw = neighborZ[ksw];
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			real mfcbb = (D.f[E   ])[k  ];
			real mfabb = (D.f[W   ])[kw ];
			real mfbcb = (D.f[N   ])[k  ];
			real mfbab = (D.f[S   ])[ks ];
			real mfbbc = (D.f[T   ])[k  ];
			real mfbba = (D.f[B   ])[kb ];
			real mfccb = (D.f[NE  ])[k  ];
			real mfaab = (D.f[SW  ])[ksw];
			real mfcab = (D.f[SE  ])[ks ];
			real mfacb = (D.f[NW  ])[kw ];
			real mfcbc = (D.f[TE  ])[k  ];
			real mfaba = (D.f[BW  ])[kbw];
			real mfcba = (D.f[BE  ])[kb ];
			real mfabc = (D.f[TW  ])[kw ];
			real mfbcc = (D.f[TN  ])[k  ];
			real mfbaa = (D.f[BS  ])[kbs];
			real mfbca = (D.f[BN  ])[kb ];
			real mfbac = (D.f[TS  ])[ks ];
			real mfbbb = (D.f[REST])[k  ];
			real mfccc = (D.f[TNE ])[k  ];
			real mfaac = (D.f[TSW ])[ksw];
			real mfcac = (D.f[TSE ])[ks ];
			real mfacc = (D.f[TNW ])[kw ];
			real mfcca = (D.f[BNE ])[kb ];
			real mfaaa = (D.f[BSW ])[kbsw];
			real mfcaa = (D.f[BSE ])[kbs];
			real mfaca = (D.f[BNW ])[kbw];
			////////////////////////////////////////////////////////////////////////////////////
			real drho = ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
							(((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
							((mfabb+mfcbb) + (mfbab+mfbcb)) + (mfbba+mfbbc)) + mfbbb;

			real rho = c1o1+drho;
			////////////////////////////////////////////////////////////////////////////////////
			real vvx    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
						     (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
						       (mfcbb-mfabb)) / rho;
			real vvy    =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
				             (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
				               (mfbcb-mfbab)) / rho;
			real vvz    =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
				             (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
				               (mfbbc-mfbba)) / rho;
			////////////////////////////////////////////////////////////////////////////////////
			real oMdrho = c1o1; // comp special
			////////////////////////////////////////////////////////////////////////////////////
			real m0, m1, m2;	
			real vx2;
			real vy2;
			real vz2;
			vx2=vvx*vvx;
			vy2=vvy*vvy;
			vz2=vvz*vvz;
			////////////////////////////////////////////////////////////////////////////////////
			//Hin
			////////////////////////////////////////////////////////////////////////////////////
			// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Z - Dir
			m2    = mfaaa	+ mfaac;
			m1    = mfaac	- mfaaa;
			m0    = m2		+ mfaab;
			mfaaa = m0;
			m0   += c1o36 * oMdrho;	
			mfaab = m1 -		m0 * vvz;
			mfaac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfabc;
			m1    = mfabc  - mfaba;
			m0    = m2		+ mfabb;
			mfaba = m0;
			m0   += c1o9 * oMdrho;
			mfabb = m1 -		m0 * vvz;
			mfabc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfacc;
			m1    = mfacc  - mfaca;
			m0    = m2		+ mfacb;
			mfaca = m0;
			m0   += c1o36 * oMdrho;
			mfacb = m1 -		m0 * vvz;
			mfacc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbac;
			m1    = mfbac	- mfbaa;
			m0    = m2		+ mfbab;
			mfbaa = m0;
			m0   += c1o9 * oMdrho;
			mfbab = m1 -		m0 * vvz;
			mfbac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbba  + mfbbc;
			m1    = mfbbc  - mfbba;
			m0    = m2		+ mfbbb;
			mfbba = m0;
			m0   += c4o9 * oMdrho;
			mfbbb = m1 -		m0 * vvz;
			mfbbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbca  + mfbcc;
			m1    = mfbcc  - mfbca;
			m0    = m2		+ mfbcb;
			mfbca = m0;
			m0   += c1o9 * oMdrho;
			mfbcb = m1 -		m0 * vvz;
			mfbcc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcac;
			m1    = mfcac	- mfcaa;
			m0    = m2		+ mfcab;
			mfcaa = m0;
			m0   += c1o36 * oMdrho;
			mfcab = m1 -		m0 * vvz;
			mfcac = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcba  + mfcbc;
			m1    = mfcbc  - mfcba;
			m0    = m2		+ mfcbb;
			mfcba = m0;
			m0   += c1o9 * oMdrho;
			mfcbb = m1 -		m0 * vvz;
			mfcbc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcca  + mfccc;
			m1    = mfccc  - mfcca;
			m0    = m2		+ mfccb;
			mfcca = m0;
			m0   += c1o36 * oMdrho;
			mfccb = m1 -		m0 * vvz;
			mfccc = m2 - c2o1*	m1 * vvz + vz2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// Y - Dir
			m2    = mfaaa	+ mfaca;
			m1    = mfaca	- mfaaa;
			m0    = m2		+ mfaba;
			mfaaa = m0;
			m0   += c1o6 * oMdrho;
			mfaba = m1 -		m0 * vvy;
			mfaca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab  + mfacb;
			m1    = mfacb  - mfaab;
			m0    = m2		+ mfabb;
			mfaab = m0;
			mfabb = m1 -		m0 * vvy;
			mfacb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac  + mfacc;
			m1    = mfacc  - mfaac;
			m0    = m2		+ mfabc;
			mfaac = m0;
			m0   += c1o18 * oMdrho;
			mfabc = m1 -		m0 * vvy;
			mfacc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbaa	+ mfbca;
			m1    = mfbca	- mfbaa;
			m0    = m2		+ mfbba;
			mfbaa = m0;
			m0   += c2o3 * oMdrho;
			mfbba = m1 -		m0 * vvy;
			mfbca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbab  + mfbcb;
			m1    = mfbcb  - mfbab;
			m0    = m2		+ mfbbb;
			mfbab = m0;
			mfbbb = m1 -		m0 * vvy;
			mfbcb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfbac  + mfbcc;
			m1    = mfbcc  - mfbac;
			m0    = m2		+ mfbbc;
			mfbac = m0;
			m0   += c2o9 * oMdrho;
			mfbbc = m1 -		m0 * vvy;
			mfbcc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcaa	+ mfcca;
			m1    = mfcca	- mfcaa;
			m0    = m2		+ mfcba;
			mfcaa = m0;
			m0   += c1o6 * oMdrho;
			mfcba = m1 -		m0 * vvy;
			mfcca = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcab  + mfccb;
			m1    = mfccb  - mfcab;
			m0    = m2		+ mfcbb;
			mfcab = m0;
			mfcbb = m1 -		m0 * vvy;
			mfccb = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfcac  + mfccc;
			m1    = mfccc  - mfcac;
			m0    = m2		+ mfcbc;
			mfcac = m0;
			m0   += c1o18 * oMdrho;
			mfcbc = m1 -		m0 * vvy;
			mfccc = m2 - c2o1*	m1 * vvy + vy2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
			////////////////////////////////////////////////////////////////////////////////////
			// X - Dir
			m2    = mfaaa	+ mfcaa;
			m1    = mfcaa	- mfaaa;
			m0    = m2		+ mfbaa;
			mfaaa = m0;
			m0   += c1o1* oMdrho;
			mfbaa = m1 -		m0 * vvx;
			mfcaa = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaba  + mfcba;
			m1    = mfcba  - mfaba;
			m0    = m2		+ mfbba;
			mfaba = m0;
			mfbba = m1 -		m0 * vvx;
			mfcba = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaca  + mfcca;
			m1    = mfcca  - mfaca;
			m0    = m2		+ mfbca;
			mfaca = m0;
			m0   += c1o3 * oMdrho;
			mfbca = m1 -		m0 * vvx;
			mfcca = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaab	+ mfcab;
			m1    = mfcab	- mfaab;
			m0    = m2		+ mfbab;
			mfaab = m0;
			mfbab = m1 -		m0 * vvx;
			mfcab = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabb  + mfcbb;
			m1    = mfcbb  - mfabb;
			m0    = m2		+ mfbbb;
			mfabb = m0;
			mfbbb = m1 -		m0 * vvx;
			mfcbb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacb  + mfccb;
			m1    = mfccb  - mfacb;
			m0    = m2		+ mfbcb;
			mfacb = m0;
			mfbcb = m1 -		m0 * vvx;
			mfccb = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfaac	+ mfcac;
			m1    = mfcac	- mfaac;
			m0    = m2		+ mfbac;
			mfaac = m0;
			m0   += c1o3 * oMdrho;
			mfbac = m1 -		m0 * vvx;
			mfcac = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfabc  + mfcbc;
			m1    = mfcbc  - mfabc;
			m0    = m2		+ mfbbc;
			mfabc = m0;
			mfbbc = m1 -		m0 * vvx;
			mfcbc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			m2    = mfacc  + mfccc;
			m1    = mfccc  - mfacc;
			m0    = m2		+ mfbcc;
			mfacc = m0;
			m0   += c1o9 * oMdrho;
			mfbcc = m1 -		m0 * vvx;
			mfccc = m2 - c2o1*	m1 * vvx + vx2 * m0;
			////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////

			real OxxPyyPzz = c1o1;
			real omega = c1o1 / (c3o1*0.001 + c1o2);
			real B = (c4o1 * omega * OxxPyyPzz * (c9o1 * omega - c16o1) - c4o1 * omega * omega - c2o1 * OxxPyyPzz * OxxPyyPzz * (c2o1 + c9o1 * omega * (omega - c2o1))) /
				(c3o1 * (omega - OxxPyyPzz) * (OxxPyyPzz * (c2o1 + c3o1 * omega) - c8o1 * omega));

			CUMbcc[k] = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)*(c1o1 + rho*c6o1*B / (c2o1 + c3o1 * B))) / rho;
			CUMcbc[k] = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)*(c1o1 + rho*c6o1*B / (c2o1 + c3o1 * B))) / rho;
			CUMccb[k] = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)*(c1o1 + rho*c6o1*B / (c2o1 + c3o1 * B))) / rho;

			////////////////////////////////////////////////////////////////////////////////////
			// Cumulants
			////////////////////////////////////////////////////////////////////////////////////
			//central moments to cumulants
			//4.
			CUMcbb[k]      = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;	
			CUMbcb[k]      = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho; 
			CUMbbc[k]      = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho; 
			 		
			CUMcca[k]      = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho  - c1o9*(drho/rho));
			CUMcac[k]      = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho  - c1o9*(drho/rho));
			CUMacc[k]      = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho  - c1o9*(drho/rho));

			//5.
			//CUMbcc[k]      = mfbcc - ((mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac) ) / rho ;
			//CUMcbc[k]      = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc) ) / rho ;
			//CUMccb[k]      = mfccb - ((mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab) ) / rho ;
			
			//6.
			CUMccc[k]      = mfccc + ((-c4o1 *  mfbbb * mfbbb  
							-           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
							-    c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
							-     c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
							+(   c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
							+     c2o1 * (mfcaa * mfaca * mfaac)
							+ c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
							-    c1o3 * (mfacc + mfcac + mfcca) /rho 
							-    c1o9 * (mfcaa + mfaca + mfaac) /rho 
							+(    c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba) 
							+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3 
							+ c1o27*((drho * drho - drho)/(rho*rho)));
			////////////////////////////////////////////////////////////////////////////////////
		}                                                                                                                    
	}
}
