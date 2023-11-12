/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
__global__ void GetVeloforForcing27( real* DD, 
                                                int* bcIndex, 
                                                int nonAtBC, 
                                                real* Vx,
                                                real* Vy,
                                                real* Vz,
                                                unsigned int* neighborX,
                                                unsigned int* neighborY,
                                                unsigned int* neighborZ,
                                                unsigned long long numberOfLBnodes, 
                                                bool isEvenTimestep)
{
    Distributions27 D;
    if (isEvenTimestep==false)
    {
        D.f[dP00] = &DD[dP00 * numberOfLBnodes];
        D.f[dM00] = &DD[dM00 * numberOfLBnodes];
        D.f[d0P0] = &DD[d0P0 * numberOfLBnodes];
        D.f[d0M0] = &DD[d0M0 * numberOfLBnodes];
        D.f[d00P] = &DD[d00P * numberOfLBnodes];
        D.f[d00M] = &DD[d00M * numberOfLBnodes];
        D.f[dPP0] = &DD[dPP0 * numberOfLBnodes];
        D.f[dMM0] = &DD[dMM0 * numberOfLBnodes];
        D.f[dPM0] = &DD[dPM0 * numberOfLBnodes];
        D.f[dMP0] = &DD[dMP0 * numberOfLBnodes];
        D.f[dP0P] = &DD[dP0P * numberOfLBnodes];
        D.f[dM0M] = &DD[dM0M * numberOfLBnodes];
        D.f[dP0M] = &DD[dP0M * numberOfLBnodes];
        D.f[dM0P] = &DD[dM0P * numberOfLBnodes];
        D.f[d0PP] = &DD[d0PP * numberOfLBnodes];
        D.f[d0MM] = &DD[d0MM * numberOfLBnodes];
        D.f[d0PM] = &DD[d0PM * numberOfLBnodes];
        D.f[d0MP] = &DD[d0MP * numberOfLBnodes];
        D.f[d000] = &DD[d000 * numberOfLBnodes];
        D.f[dPPP] = &DD[dPPP * numberOfLBnodes];
        D.f[dMMP] = &DD[dMMP * numberOfLBnodes];
        D.f[dPMP] = &DD[dPMP * numberOfLBnodes];
        D.f[dMPP] = &DD[dMPP * numberOfLBnodes];
        D.f[dPPM] = &DD[dPPM * numberOfLBnodes];
        D.f[dMMM] = &DD[dMMM * numberOfLBnodes];
        D.f[dPMM] = &DD[dPMM * numberOfLBnodes];
        D.f[dMPM] = &DD[dMPM * numberOfLBnodes];
    } 
    else
    {
        D.f[dM00] = &DD[dP00 * numberOfLBnodes];
        D.f[dP00] = &DD[dM00 * numberOfLBnodes];
        D.f[d0M0] = &DD[d0P0 * numberOfLBnodes];
        D.f[d0P0] = &DD[d0M0 * numberOfLBnodes];
        D.f[d00M] = &DD[d00P * numberOfLBnodes];
        D.f[d00P] = &DD[d00M * numberOfLBnodes];
        D.f[dMM0] = &DD[dPP0 * numberOfLBnodes];
        D.f[dPP0] = &DD[dMM0 * numberOfLBnodes];
        D.f[dMP0] = &DD[dPM0 * numberOfLBnodes];
        D.f[dPM0] = &DD[dMP0 * numberOfLBnodes];
        D.f[dM0M] = &DD[dP0P * numberOfLBnodes];
        D.f[dP0P] = &DD[dM0M * numberOfLBnodes];
        D.f[dM0P] = &DD[dP0M * numberOfLBnodes];
        D.f[dP0M] = &DD[dM0P * numberOfLBnodes];
        D.f[d0MM] = &DD[d0PP * numberOfLBnodes];
        D.f[d0PP] = &DD[d0MM * numberOfLBnodes];
        D.f[d0MP] = &DD[d0PM * numberOfLBnodes];
        D.f[d0PM] = &DD[d0MP * numberOfLBnodes];
        D.f[d000] = &DD[d000 * numberOfLBnodes];
        D.f[dPPP] = &DD[dMMM * numberOfLBnodes];
        D.f[dMMP] = &DD[dPPM * numberOfLBnodes];
        D.f[dPMP] = &DD[dMPM * numberOfLBnodes];
        D.f[dMPP] = &DD[dPMM * numberOfLBnodes];
        D.f[dPPM] = &DD[dMMP * numberOfLBnodes];
        D.f[dMMM] = &DD[dPPP * numberOfLBnodes];
        D.f[dPMM] = &DD[dMPP * numberOfLBnodes];
        D.f[dMPM] = &DD[dPMP * numberOfLBnodes];
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
        real mfcbb = (D.f[dP00])[ke   ];
        real mfabb = (D.f[dM00])[kw   ];
        real mfbcb = (D.f[d0P0])[kn   ];
        real mfbab = (D.f[d0M0])[ks   ];
        real mfbbc = (D.f[d00P])[kt   ];
        real mfbba = (D.f[d00M])[kb   ];
        real mfccb = (D.f[dPP0])[kne  ];
        real mfaab = (D.f[dMM0])[ksw  ];
        real mfcab = (D.f[dPM0])[kse  ];
        real mfacb = (D.f[dMP0])[knw  ];
        real mfcbc = (D.f[dP0P])[kte  ];
        real mfaba = (D.f[dM0M])[kbw  ];
        real mfcba = (D.f[dP0M])[kbe  ];
        real mfabc = (D.f[dM0P])[ktw  ];
        real mfbcc = (D.f[d0PP])[ktn  ];
        real mfbaa = (D.f[d0MM])[kbs  ];
        real mfbca = (D.f[d0PM])[kbn  ];
        real mfbac = (D.f[d0MP])[kts  ];
        real mfbbb = (D.f[d000])[kzero];
        real mfccc = (D.f[dPPP])[ktne ];
        real mfaac = (D.f[dMMP])[ktsw ];
        real mfcac = (D.f[dPMP])[ktse ];
        real mfacc = (D.f[dMPP])[ktnw ];
        real mfcca = (D.f[dPPM])[kbne ];
        real mfaaa = (D.f[dMMM])[kbsw ];
        real mfcaa = (D.f[dPMM])[kbse ];
        real mfaca = (D.f[dMPM])[kbnw ];
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

