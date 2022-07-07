/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
extern "C" __global__ void GetVeloforForcing27( real* DD, 
												int* bcIndex, 
												int nonAtBC, 
												real* Vx,
												real* Vy,
												real* Vz,
												unsigned int* neighborX,
												unsigned int* neighborY,
												unsigned int* neighborZ,
												unsigned int size_Mat, 
												bool isEvenTimestep)
{
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
	////////////////////////////////////////////////////////////////////////////////
	const unsigned  x = threadIdx.x;  // Globaler x-Index 
	const unsigned  y = blockIdx.x;   // Globaler y-Index 
	const unsigned  z = blockIdx.y;   // Globaler z-Index 

	const unsigned nx = blockDim.x;
	const unsigned ny = gridDim.x;

	const unsigned k = nx*(ny*z + y) + x;
	//////////////////////////////////////////////////////////////////////////
	if(k < nonAtBC)
	{
		////////////////////////////////////////////////////////////////////////////////
		//index
		unsigned int KQK  = bcIndex[k];
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
		real mfcbb = (D.f[E   ])[ke   ];
		real mfabb = (D.f[W   ])[kw   ];
		real mfbcb = (D.f[N   ])[kn   ];
		real mfbab = (D.f[S   ])[ks   ];
		real mfbbc = (D.f[T   ])[kt   ];
		real mfbba = (D.f[B   ])[kb   ];
		real mfccb = (D.f[NE  ])[kne  ];
		real mfaab = (D.f[SW  ])[ksw  ];
		real mfcab = (D.f[SE  ])[kse  ];
		real mfacb = (D.f[NW  ])[knw  ];
		real mfcbc = (D.f[TE  ])[kte  ];
		real mfaba = (D.f[BW  ])[kbw  ];
		real mfcba = (D.f[BE  ])[kbe  ];
		real mfabc = (D.f[TW  ])[ktw  ];
		real mfbcc = (D.f[TN  ])[ktn  ];
		real mfbaa = (D.f[BS  ])[kbs  ];
		real mfbca = (D.f[BN  ])[kbn  ];
		real mfbac = (D.f[TS  ])[kts  ];
		real mfbbb = (D.f[REST])[kzero];
		real mfccc = (D.f[TNE ])[ktne ];
		real mfaac = (D.f[TSW ])[ktsw ];
		real mfcac = (D.f[TSE ])[ktse ];
		real mfacc = (D.f[TNW ])[ktnw ];
		real mfcca = (D.f[BNE ])[kbne ];
		real mfaaa = (D.f[BSW ])[kbsw ];
		real mfcaa = (D.f[BSE ])[kbse ];
		real mfaca = (D.f[BNW ])[kbnw ];
		////////////////////////////////////////////////////////////////////////////////////
		real rho   = (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
					 	 mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
						 mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb + c1o1);//!!!!Achtung + one
		////////////////////////////////////////////////////////////////////////////////////
		real vx =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfcaa-mfacc) + (mfcca-mfaac))) + 
			         (((mfcba-mfabc) + (mfcbc-mfaba)) + ((mfcab-mfacb) + (mfccb-mfaab))) +
				       (mfcbb-mfabb))/ rho;
		real vy =((((mfccc-mfaaa) + (mfaca-mfcac)) + ((mfacc-mfcaa) + (mfcca-mfaac))) + 
			         (((mfbca-mfbac) + (mfbcc-mfbaa)) + ((mfacb-mfcab) + (mfccb-mfaab))) +
			           (mfbcb-mfbab)) / rho;
		real vz =((((mfccc-mfaaa) + (mfcac-mfaca)) + ((mfacc-mfcaa) + (mfaac-mfcca))) + 
			         (((mfbac-mfbca) + (mfbcc-mfbaa)) + ((mfabc-mfcba) + (mfcbc-mfaba))) +
			           (mfbbc-mfbba)) / rho;
		////////////////////////////////////////////////////////////////////////////////////
		Vx[k] = vx;
		Vy[k] = vy;
		Vz[k] = vz;
		////////////////////////////////////////////////////////////////////////////////////
	}
}

