#include "Calculation/Calculation.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void B92IncompressibleNavierStokes_Device(real omega,
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

    if (k < size_Mat)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if ((BC != GEO_SOLID) && (BC != GEO_VOID))
        {
            Distributions27 D;
            if (EvenOrOdd == true)
            {
                D.f[dP00] = &DDStart[dP00 * size_Mat];
                D.f[dM00] = &DDStart[dM00 * size_Mat];
                D.f[d0P0] = &DDStart[d0P0 * size_Mat];
                D.f[d0M0] = &DDStart[d0M0 * size_Mat];
                D.f[d00P] = &DDStart[d00P * size_Mat];
                D.f[d00M] = &DDStart[d00M * size_Mat];
                D.f[dPP0] = &DDStart[dPP0 * size_Mat];
                D.f[dMM0] = &DDStart[dMM0 * size_Mat];
                D.f[dPM0] = &DDStart[dPM0 * size_Mat];
                D.f[dMP0] = &DDStart[dMP0 * size_Mat];
                D.f[dP0P] = &DDStart[dP0P * size_Mat];
                D.f[dM0M] = &DDStart[dM0M * size_Mat];
                D.f[dP0M] = &DDStart[dP0M * size_Mat];
                D.f[dM0P] = &DDStart[dM0P * size_Mat];
                D.f[d0PP] = &DDStart[d0PP * size_Mat];
                D.f[d0MM] = &DDStart[d0MM * size_Mat];
                D.f[d0PM] = &DDStart[d0PM * size_Mat];
                D.f[d0MP] = &DDStart[d0MP * size_Mat];
                D.f[d000] = &DDStart[d000 * size_Mat];
                D.f[dPPP] = &DDStart[dPPP * size_Mat];
                D.f[dMMP] = &DDStart[dMMP * size_Mat];
                D.f[dPMP] = &DDStart[dPMP * size_Mat];
                D.f[dMPP] = &DDStart[dMPP * size_Mat];
                D.f[dPPM] = &DDStart[dPPM * size_Mat];
                D.f[dMMM] = &DDStart[dMMM * size_Mat];
                D.f[dPMM] = &DDStart[dPMM * size_Mat];
                D.f[dMPM] = &DDStart[dMPM * size_Mat];
            }
            else
            {
                D.f[dM00] = &DDStart[dP00 * size_Mat];
                D.f[dP00] = &DDStart[dM00 * size_Mat];
                D.f[d0M0] = &DDStart[d0P0 * size_Mat];
                D.f[d0P0] = &DDStart[d0M0 * size_Mat];
                D.f[d00M] = &DDStart[d00P * size_Mat];
                D.f[d00P] = &DDStart[d00M * size_Mat];
                D.f[dMM0] = &DDStart[dPP0 * size_Mat];
                D.f[dPP0] = &DDStart[dMM0 * size_Mat];
                D.f[dMP0] = &DDStart[dPM0 * size_Mat];
                D.f[dPM0] = &DDStart[dMP0 * size_Mat];
                D.f[dM0M] = &DDStart[dP0P * size_Mat];
                D.f[dP0P] = &DDStart[dM0M * size_Mat];
                D.f[dM0P] = &DDStart[dP0M * size_Mat];
                D.f[dP0M] = &DDStart[dM0P * size_Mat];
                D.f[d0MM] = &DDStart[d0PP * size_Mat];
                D.f[d0PP] = &DDStart[d0MM * size_Mat];
                D.f[d0MP] = &DDStart[d0PM * size_Mat];
                D.f[d0PM] = &DDStart[d0MP * size_Mat];
                D.f[d000] = &DDStart[d000 * size_Mat];
                D.f[dMMM] = &DDStart[dPPP * size_Mat];
                D.f[dPPM] = &DDStart[dMMP * size_Mat];
                D.f[dMPM] = &DDStart[dPMP * size_Mat];
                D.f[dPMM] = &DDStart[dMPP * size_Mat];
                D.f[dMMP] = &DDStart[dPPM * size_Mat];
                D.f[dPPP] = &DDStart[dMMM * size_Mat];
                D.f[dMPP] = &DDStart[dPMM * size_Mat];
                D.f[dPMP] = &DDStart[dMPM * size_Mat];
            }

            ////////////////////////////////////////////////////////////////////////////////
            //index
            //unsigned int kzero= k;
            //unsigned int ke   = k;
            unsigned int kw = neighborX[k];
            //unsigned int kn   = k;
            unsigned int ks = neighborY[k];
            //unsigned int kt   = k;
            unsigned int kb = neighborZ[k];
            unsigned int ksw = neighborY[kw];
            //unsigned int kne  = k;
            //unsigned int kse  = ks;
            //unsigned int knw  = kw;
            unsigned int kbw = neighborZ[kw];
            //unsigned int kte  = k;
            //unsigned int kbe  = kb;
            //unsigned int ktw  = kw;
            unsigned int kbs = neighborZ[ks];
            //unsigned int ktn  = k;
            //unsigned int kbn  = kb;
            //unsigned int kts  = ks;
            //unsigned int ktse = ks;
            //unsigned int kbnw = kbw;
            //unsigned int ktnw = kw;
            //unsigned int kbse = kbs;
            //unsigned int ktsw = ksw;
            //unsigned int kbne = kb;
            //unsigned int ktne = k;
            unsigned int kbsw = neighborZ[ksw];
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real fE = (D.f[dP00])[k];//ke
            real fW = (D.f[dM00])[kw];
            real fN = (D.f[d0P0])[k];//kn
            real fS = (D.f[d0M0])[ks];
            real fT = (D.f[d00P])[k];//kt
            real fB = (D.f[d00M])[kb];
            real fNE = (D.f[dPP0])[k];//kne
            real fSW = (D.f[dMM0])[ksw];
            real fSE = (D.f[dPM0])[ks];//kse
            real fNW = (D.f[dMP0])[kw];//knw
            real fTE = (D.f[dP0P])[k];//kte
            real fBW = (D.f[dM0M])[kbw];
            real fBE = (D.f[dP0M])[kb];//kbe
            real fTW = (D.f[dM0P])[kw];//ktw
            real fTN = (D.f[d0PP])[k];//ktn
            real fBS = (D.f[d0MM])[kbs];
            real fBN = (D.f[d0PM])[kb];//kbn
            real fTS = (D.f[d0MP])[ks];//kts
            real fZERO = (D.f[d000])[k];//kzero
            real fTNE = (D.f[dPPP])[k];//ktne
            real fTSW = (D.f[dMMP])[ksw];//ktsw
            real fTSE = (D.f[dPMP])[ks];//ktse
            real fTNW = (D.f[dMPP])[kw];//ktnw
            real fBNE = (D.f[dPPM])[kb];//kbne
            real fBSW = (D.f[dMMM])[kbsw];
            real fBSE = (D.f[dPMM])[kbs];//kbse
            real fBNW = (D.f[dMPM])[kbw];//kbnw
                                           ////////////////////////////////////////////////////////////////////////////////







                                           //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                           //BGK incomp
                                           //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real drho = (fTNE + fBSW) + (fTSW + fBNE) + (fTSE + fBNW) + (fTNW + fBSE) + (fNE + fSW) + (fNW + fSE) + (fTE + fBW) + (fBE + fTW) + (fTN + fBS) + (fBN + fTS) + (fE + fW) + (fN + fS) + (fT + fB) + fZERO;
            real vx1 = (fTNE - fBSW) + (fBNE - fTSW) + (fTSE - fBNW) + (fBSE - fTNW) + (fNE - fSW) + (fSE - fNW) + (fTE - fBW) + (fBE - fTW) + (fE - fW);
            real vx2 = (fTNE - fBSW) + (fBNE - fTSW) + (fBNW - fTSE) + (fTNW - fBSE) + (fNE - fSW) + (fNW - fSE) + (fTN - fBS) + (fBN - fTS) + (fN - fS);
            real vx3 = (fTNE - fBSW) + (fTSW - fBNE) + (fTSE - fBNW) + (fTNW - fBSE) + (fTE - fBW) + (fTW - fBE) + (fTN - fBS) + (fTS - fBN) + (fT - fB);
            real cusq = c3o2*(vx1*vx1 + vx2*vx2 + vx3*vx3);
            //////////////////////////////////////////////////////////////////////////                            
            fZERO = fZERO *(c1o1 + (-omega)) - (-omega)*   c8o27*  (drho - cusq);
            fE = fE    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(vx1)+c9o2*(vx1)*(vx1)-cusq);
            fW = fW    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(-vx1) + c9o2*(-vx1)*(-vx1) - cusq);
            fN = fN    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(vx2)+c9o2*(vx2)*(vx2)-cusq);
            fS = fS    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(-vx2) + c9o2*(-vx2)*(-vx2) - cusq);
            fT = fT    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(vx3)+c9o2*(vx3)*(vx3)-cusq);
            fB = fB    *(c1o1 + (-omega)) - (-omega)*   c2o27*  (drho + c3o1*(-vx3) + c9o2*(-vx3)*(-vx3) - cusq);
            fNE = fNE   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + c3o1*(vx1 + vx2) + c9o2*(vx1 + vx2)*(vx1 + vx2) - cusq);
            fSW = fSW   *(c1o1 + (-omega)) - (-omega)*   c1o54*  (drho + c3o1*(-vx1 - vx2) + c9o2*(-vx1 - vx2)*(-vx1 - vx2) - cusq);
            fSE = fSE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx1 - vx2) + c9o2*(vx1 - vx2)*(vx1 - vx2) - cusq);
            fNW = fNW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx1 + vx2) + c9o2*(-vx1 + vx2)*(-vx1 + vx2) - cusq);
            fTE = fTE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx1 + vx3) + c9o2*(vx1 + vx3)*(vx1 + vx3) - cusq);
            fBW = fBW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx1 - vx3) + c9o2*(-vx1 - vx3)*(-vx1 - vx3) - cusq);
            fBE = fBE   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx1 - vx3) + c9o2*(vx1 - vx3)*(vx1 - vx3) - cusq);
            fTW = fTW   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx1 + vx3) + c9o2*(-vx1 + vx3)*(-vx1 + vx3) - cusq);
            fTN = fTN   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx2 + vx3) + c9o2*(vx2 + vx3)*(vx2 + vx3) - cusq);
            fBS = fBS   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx2 - vx3) + c9o2*(-vx2 - vx3)*(-vx2 - vx3) - cusq);
            fBN = fBN   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(vx2 - vx3) + c9o2*(vx2 - vx3)*(vx2 - vx3) - cusq);
            fTS = fTS   *(c1o1 + (-omega)) - (-omega)*    c1o54* (drho + c3o1*(-vx2 + vx3) + c9o2*(-vx2 + vx3)*(-vx2 + vx3) - cusq);
            fTNE = fTNE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 + vx2 + vx3) + c9o2*(vx1 + vx2 + vx3)*(vx1 + vx2 + vx3) - cusq);
            fBSW = fBSW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 - vx2 - vx3) + c9o2*(-vx1 - vx2 - vx3)*(-vx1 - vx2 - vx3) - cusq);
            fBNE = fBNE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 + vx2 - vx3) + c9o2*(vx1 + vx2 - vx3)*(vx1 + vx2 - vx3) - cusq);
            fTSW = fTSW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 - vx2 + vx3) + c9o2*(-vx1 - vx2 + vx3)*(-vx1 - vx2 + vx3) - cusq);
            fTSE = fTSE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 - vx2 + vx3) + c9o2*(vx1 - vx2 + vx3)*(vx1 - vx2 + vx3) - cusq);
            fBNW = fBNW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 + vx2 - vx3) + c9o2*(-vx1 + vx2 - vx3)*(-vx1 + vx2 - vx3) - cusq);
            fBSE = fBSE  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(vx1 - vx2 - vx3) + c9o2*(vx1 - vx2 - vx3)*(vx1 - vx2 - vx3) - cusq);
            fTNW = fTNW  *(c1o1 + (-omega)) - (-omega)*    c1o216*(drho + c3o1*(-vx1 + vx2 + vx3) + c9o2*(-vx1 + vx2 + vx3)*(-vx1 + vx2 + vx3) - cusq);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







            //////////////////////////////////////////////////////////////////////////                            
            (D.f[dP00])[k] = fW;
            (D.f[dM00])[kw] = fE;
            (D.f[d0P0])[k] = fS;
            (D.f[d0M0])[ks] = fN;
            (D.f[d00P])[k] = fB;
            (D.f[d00M])[kb] = fT;
            (D.f[dPP0])[k] = fSW;
            (D.f[dMM0])[ksw] = fNE;
            (D.f[dPM0])[ks] = fNW;
            (D.f[dMP0])[kw] = fSE;
            (D.f[dP0P])[k] = fBW;
            (D.f[dM0M])[kbw] = fTE;
            (D.f[dP0M])[kb] = fTW;
            (D.f[dM0P])[kw] = fBE;
            (D.f[d0PP])[k] = fBS;
            (D.f[d0MM])[kbs] = fTN;
            (D.f[d0PM])[kb] = fTS;
            (D.f[d0MP])[ks] = fBN;
            (D.f[d000])[k] = fZERO;
            (D.f[dPPP])[k] = fBSW;
            (D.f[dPMP])[ks] = fBNW;
            (D.f[dPPM])[kb] = fTSW;
            (D.f[dPMM])[kbs] = fTNW;
            (D.f[dMPP])[kw] = fBSE;
            (D.f[dMMP])[ksw] = fBNE;
            (D.f[dMPM])[kbw] = fTSE;
            (D.f[dMMM])[kbsw] = fTNE;
            //////////////////////////////////////////////////////////////////////////                            
        }
    }
}