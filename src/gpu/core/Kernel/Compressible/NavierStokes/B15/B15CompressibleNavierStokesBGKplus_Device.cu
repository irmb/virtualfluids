#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void B15CompressibleNavierStokesBGKplus_Device(
    real omega,
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
            real mfcbb = (D.f[dP00])[k];//[ke   ];// +  c2over27 ;(D.f[dP00])[k  ];//ke
            real mfabb = (D.f[dM00])[kw];//[kw   ];// +  c2over27 ;(D.f[dM00])[kw ];
            real mfbcb = (D.f[d0P0])[k];//[kn   ];// +  c2over27 ;(D.f[d0P0])[k  ];//kn
            real mfbab = (D.f[d0M0])[ks];//[ks   ];// +  c2over27 ;(D.f[d0M0])[ks ];
            real mfbbc = (D.f[d00P])[k];//[kt   ];// +  c2over27 ;(D.f[d00P])[k  ];//kt
            real mfbba = (D.f[d00M])[kb];//[kb   ];// +  c2over27 ;(D.f[d00M])[kb ];
            real mfccb = (D.f[dPP0])[k];//[kne  ];// +  c1over54 ;(D.f[dPP0])[k  ];//kne
            real mfaab = (D.f[dMM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[dMM0])[ksw];
            real mfcab = (D.f[dPM0])[ks];//[kse  ];// +  c1over54 ;(D.f[dPM0])[ks ];//kse
            real mfacb = (D.f[dMP0])[kw];//[knw  ];// +  c1over54 ;(D.f[dMP0])[kw ];//knw
            real mfcbc = (D.f[dP0P])[k];//[kte  ];// +  c1over54 ;(D.f[dP0P])[k  ];//kte
            real mfaba = (D.f[dM0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[dM0M])[kbw];
            real mfcba = (D.f[dP0M])[kb];//[kbe  ];// +  c1over54 ;(D.f[dP0M])[kb ];//kbe
            real mfabc = (D.f[dM0P])[kw];//[ktw  ];// +  c1over54 ;(D.f[dM0P])[kw ];//ktw
            real mfbcc = (D.f[d0PP])[k];//[ktn  ];// +  c1over54 ;(D.f[d0PP])[k  ];//ktn
            real mfbaa = (D.f[d0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[d0MM])[kbs];
            real mfbca = (D.f[d0PM])[kb];//[kbn  ];// +  c1over54 ;(D.f[d0PM])[kb ];//kbn
            real mfbac = (D.f[d0MP])[ks];//[kts  ];// +  c1over54 ;(D.f[d0MP])[ks ];//kts
            real mfbbb = (D.f[d000])[k];//[kzero];// +  c8over27 ;(D.f[d000])[k  ];//kzero
            real mfccc = (D.f[dPPP])[k];//[ktne ];// +  c1over216;(D.f[dPPP])[k  ];//ktne
            real mfaac = (D.f[dMMP])[ksw];//[ktsw ];// +  c1over216;(D.f[dMMP])[ksw];//ktsw
            real mfcac = (D.f[dPMP])[ks];//[ktse ];// +  c1over216;(D.f[dPMP])[ks ];//ktse
            real mfacc = (D.f[dMPP])[kw];//[ktnw ];// +  c1over216;(D.f[dMPP])[kw ];//ktnw
            real mfcca = (D.f[dPPM])[kb];//[kbne ];// +  c1over216;(D.f[dPPM])[kb ];//kbne
            real mfaaa = (D.f[dMMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[dMMM])[kbsw];
            real mfcaa = (D.f[dPMM])[kbs];//[kbse ];// +  c1over216;(D.f[dPMM])[kbs];//kbse
            real mfaca = (D.f[dMPM])[kbw];//[kbnw ];// +  c1over216;(D.f[dMPM])[kbw];//kbnw
                                            ////////////////////////////////////////////////////////////////////////////////////
                                            //slow
                                            //real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
                                            //                       (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
                                            //                        ((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
                                            ////////////////////////////////////////////////////////////////////////////////////
            real rho = (mfccc + mfaaa + mfaca + mfcac + mfacc + mfcaa + mfaac + mfcca +
                mfbac + mfbca + mfbaa + mfbcc + mfabc + mfcba + mfaba + mfcbc + mfacb + mfcab + mfaab + mfccb +
                mfabb + mfcbb + mfbab + mfbcb + mfbba + mfbbc + mfbbb + c1o1);//!!!!Achtung + one
                                                                             ////////////////////////////////////////////////////////////////////////////////////
            real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                (mfcbb - mfabb)) / rho;
            real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                (mfbcb - mfbab)) / rho;
            real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                (mfbbc - mfbba)) / rho;
            ////////////////////////////////////////////////////////////////////////////////////
            real vx2 = vvx * vvx;
            real vy2 = vvy * vvy;
            real vz2 = vvz * vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            real m200 = (mfccc + mfaaa + mfaca + mfcac + mfacc + mfcaa + mfaac + mfcca +
                mfabc + mfcba + mfaba + mfcbc + mfacb + mfcab + mfaab + mfccb +
                mfabb + mfcbb);
            real m020 = (mfccc + mfaaa + mfaca + mfcac + mfacc + mfcaa + mfaac + mfcca +
                mfbac + mfbca + mfbaa + mfbcc + mfacb + mfcab + mfaab + mfccb +
                mfbab + mfbcb);
            real m002 = (mfccc + mfaaa + mfaca + mfcac + mfacc + mfcaa + mfaac + mfcca +
                mfbac + mfbca + mfbaa + mfbcc + mfabc + mfcba + mfaba + mfcbc +
                mfbba + mfbbc);
            ////////////////////////////////////////////////////////////////////////////////////
            //Galilei Korrektur
            real Gx = -c3o1 * vx2 * (-c1o2 * (c3o1 * m200 / rho + c1o1 / rho - c1o1 - c3o1 * vx2)) * (c1o1 - omega * c1o2);
            real Gy = -c3o1 * vy2 * (-c1o2 * (c3o1 * m020 / rho + c1o1 / rho - c1o1 - c3o1 * vy2)) * (c1o1 - omega * c1o2);
            real Gz = -c3o1 * vz2 * (-c1o2 * (c3o1 * m002 / rho + c1o1 / rho - c1o1 - c3o1 * vz2)) * (c1o1 - omega * c1o2);
            //real Gx     = zero;
            //real Gy     = zero;
            //real Gz     = zero;
            ////////////////////////////////////////////////////////////////////////////////////
            real XXb = -c2o3 + vx2 + Gx;
            real XXc = -c1o2 * (XXb + c1o1 + vvx);
            real XXa = XXc + vvx;
            real YYb = -c2o3 + vy2 + Gy;
            real YYc = -c1o2 * (YYb + c1o1 + vvy);
            real YYa = YYc + vvy;
            real ZZb = -c2o3 + vz2 + Gz;
            real ZZc = -c1o2 * (ZZb + c1o1 + vvz);
            real ZZa = ZZc + vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            mfcbb = mfcbb * (c1o1 - omega) + omega * (-rho * XXc * YYb * ZZb - c2o27);
            mfabb = mfabb * (c1o1 - omega) + omega * (-rho * XXa * YYb * ZZb - c2o27);
            mfbcb = mfbcb * (c1o1 - omega) + omega * (-rho * XXb * YYc * ZZb - c2o27);
            mfbab = mfbab * (c1o1 - omega) + omega * (-rho * XXb * YYa * ZZb - c2o27);
            mfbbc = mfbbc * (c1o1 - omega) + omega * (-rho * XXb * YYb * ZZc - c2o27);
            mfbba = mfbba * (c1o1 - omega) + omega * (-rho * XXb * YYb * ZZa - c2o27);
            mfccb = mfccb * (c1o1 - omega) + omega * (-rho * XXc * YYc * ZZb - c1o54);
            mfaab = mfaab * (c1o1 - omega) + omega * (-rho * XXa * YYa * ZZb - c1o54);
            mfcab = mfcab * (c1o1 - omega) + omega * (-rho * XXc * YYa * ZZb - c1o54);
            mfacb = mfacb * (c1o1 - omega) + omega * (-rho * XXa * YYc * ZZb - c1o54);
            mfcbc = mfcbc * (c1o1 - omega) + omega * (-rho * XXc * YYb * ZZc - c1o54);
            mfaba = mfaba * (c1o1 - omega) + omega * (-rho * XXa * YYb * ZZa - c1o54);
            mfcba = mfcba * (c1o1 - omega) + omega * (-rho * XXc * YYb * ZZa - c1o54);
            mfabc = mfabc * (c1o1 - omega) + omega * (-rho * XXa * YYb * ZZc - c1o54);
            mfbcc = mfbcc * (c1o1 - omega) + omega * (-rho * XXb * YYc * ZZc - c1o54);
            mfbaa = mfbaa * (c1o1 - omega) + omega * (-rho * XXb * YYa * ZZa - c1o54);
            mfbca = mfbca * (c1o1 - omega) + omega * (-rho * XXb * YYc * ZZa - c1o54);
            mfbac = mfbac * (c1o1 - omega) + omega * (-rho * XXb * YYa * ZZc - c1o54);
            mfbbb = mfbbb * (c1o1 - omega) + omega * (-rho * XXb * YYb * ZZb - c8o27);
            mfccc = mfccc * (c1o1 - omega) + omega * (-rho * XXc * YYc * ZZc - c1o216);
            mfaac = mfaac * (c1o1 - omega) + omega * (-rho * XXa * YYa * ZZc - c1o216);
            mfcac = mfcac * (c1o1 - omega) + omega * (-rho * XXc * YYa * ZZc - c1o216);
            mfacc = mfacc * (c1o1 - omega) + omega * (-rho * XXa * YYc * ZZc - c1o216);
            mfcca = mfcca * (c1o1 - omega) + omega * (-rho * XXc * YYc * ZZa - c1o216);
            mfaaa = mfaaa * (c1o1 - omega) + omega * (-rho * XXa * YYa * ZZa - c1o216);
            mfcaa = mfcaa * (c1o1 - omega) + omega * (-rho * XXc * YYa * ZZa - c1o216);
            mfaca = mfaca * (c1o1 - omega) + omega * (-rho * XXa * YYc * ZZa - c1o216);
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            //fast
            //            real oMdrho = one; //comp special
            //            //real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
            //            //                       mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
            //            //                       mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb + one);//fehlt mfbbb nicht mehr !!!!Achtung + one
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            real m0, m1, m2;    
            //            real vx2;
            //            real vy2;
            //            real vz2;
            //            vx2=vvx*vvx;
            //            vy2=vvy*vvy;
            //            vz2=vvz*vvz;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            real wadjust;
            //            real qudricLimit = 0.01f;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            //Hin
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // Z - Dir
            //            m2    = mfaaa    + mfaac;
            //            m1    = mfaac    - mfaaa;
            //            m0    = m2        + mfaab;
            //            mfaaa = m0;
            //            m0   += c1o36 * oMdrho;    
            //            mfaab = m1 ;
            //            mfaac = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaba  + mfabc;
            //            m1    = mfabc  - mfaba;
            //            m0    = m2        + mfabb;
            //            mfaba = m0;
            //            m0   += c1o9 * oMdrho;
            //            mfabb = m1 ;
            //            mfabc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaca  + mfacc;
            //            m1    = mfacc  - mfaca;
            //            m0    = m2        + mfacb;
            //            mfaca = m0;
            //            m0   += c1o36 * oMdrho;
            //            mfacb = m1 ;
            //            mfacc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfbaa    + mfbac;
            //            m1    = mfbac    - mfbaa;
            //            m0    = m2        + mfbab;
            //            mfbaa = m0;
            //            m0   += c1o9 * oMdrho;
            //            mfbab = m1 ;
            //            mfbac = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfbba  + mfbbc;
            //            m1    = mfbbc  - mfbba;
            //            m0    = m2        + mfbbb;
            //            mfbba = m0;
            //            m0   += c4o9 * oMdrho;
            //            mfbbb = m1 ;
            //            mfbbc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfbca  + mfbcc;
            //            m1    = mfbcc  - mfbca;
            //            m0    = m2        + mfbcb;
            //            mfbca = m0;
            //            m0   += c1o9 * oMdrho;
            //            mfbcb = m1 ;
            //            mfbcc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfcaa    + mfcac;
            //            m1    = mfcac    - mfcaa;
            //            m0    = m2        + mfcab;
            //            mfcaa = m0;
            //            m0   += c1o36 * oMdrho;
            //            mfcab = m1 ;
            //            mfcac = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfcba  + mfcbc;
            //            m1    = mfcbc  - mfcba;
            //            m0    = m2        + mfcbb;
            //            mfcba = m0;
            //            m0   += c1o9 * oMdrho;
            //            mfcbb = m1 ;
            //            mfcbc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfcca  + mfccc;
            //            m1    = mfccc  - mfcca;
            //            m0    = m2        + mfccb;
            //            mfcca = m0;
            //            m0   += c1o36 * oMdrho;
            //            mfccb = m1 ;
            //            mfccc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // Y - Dir
            //            m2    = mfaaa    + mfaca;
            //            m1    = mfaca    - mfaaa;
            //            m0    = m2        + mfaba;
            //            mfaaa = m0;
            //            m0   += c1o6 * oMdrho;
            //            mfaba = m1 ;
            //            mfaca = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaab  + mfacb;
            //            m1    = mfacb  - mfaab;
            //            m0    = m2        + mfabb;
            //            mfaab = m0;
            //            mfabb = m1 ;
            //            mfacb = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaac  + mfacc;
            //            m1    = mfacc  - mfaac;
            //            m0    = m2        + mfabc;
            //            mfaac = m0;
            //            m0   += c1o18 * oMdrho;
            //            mfabc = m1 ;
            //            mfacc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfbaa    + mfbca;
            //            m1    = mfbca    - mfbaa;
            //            m0    = m2        + mfbba;
            //            mfbaa = m0;
            //            m0   += c2o3 * oMdrho;
            //            mfbba = m1 ;
            //            mfbca = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfbab  + mfbcb;
            //            m1    = mfbcb  - mfbab;
            //            m0    = m2        + mfbbb;
            //            mfbab = m0;
            //            mfbbb = m1 ;
            //            mfbcb = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfbac  + mfbcc;
            //            m1    = mfbcc  - mfbac;
            //            m0    = m2        + mfbbc;
            //            mfbac = m0;
            //            m0   += c2o9 * oMdrho;
            //            mfbbc = m1 ;
            //            mfbcc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfcaa    + mfcca;
            //            m1    = mfcca    - mfcaa;
            //            m0    = m2        + mfcba;
            //            mfcaa = m0;
            //            m0   += c1o6 * oMdrho;
            //            mfcba = m1 ;
            //            mfcca = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfcab  + mfccb;
            //            m1    = mfccb  - mfcab;
            //            m0    = m2        + mfcbb;
            //            mfcab = m0;
            //            mfcbb = m1 ;
            //            mfccb = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfcac  + mfccc;
            //            m1    = mfccc  - mfcac;
            //            m0    = m2        + mfcbc;
            //            mfcac = m0;
            //            m0   += c1o18 * oMdrho;
            //            mfcbc = m1 ;
            //            mfccc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9        Konditionieren
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // X - Dir
            //            m2    = mfaaa    + mfcaa;
            //            m1    = mfcaa    - mfaaa;
            //            m0    = m2        + mfbaa;
            //            mfaaa = m0;
            //            m0   += one* oMdrho;
            //            mfbaa = m1 ;
            //            mfcaa = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaba  + mfcba;
            //            m1    = mfcba  - mfaba;
            //            m0    = m2        + mfbba;
            //            mfaba = m0;
            //            mfbba = m1 ;
            //            mfcba = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaca  + mfcca;
            //            m1    = mfcca  - mfaca;
            //            m0    = m2        + mfbca;
            //            mfaca = m0;
            //            m0   += c1o3 * oMdrho;
            //            mfbca = m1 ;
            //            mfcca = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaab    + mfcab;
            //            m1    = mfcab    - mfaab;
            //            m0    = m2        + mfbab;
            //            mfaab = m0;
            //            mfbab = m1 ;
            //            mfcab = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfabb  + mfcbb;
            //            m1    = mfcbb  - mfabb;
            //            m0    = m2        + mfbbb;
            //            mfabb = m0;
            //            mfbbb = m1 ;
            //            mfcbb = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfacb  + mfccb;
            //            m1    = mfccb  - mfacb;
            //            m0    = m2        + mfbcb;
            //            mfacb = m0;
            //            mfbcb = m1 ;
            //            mfccb = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfaac    + mfcac;
            //            m1    = mfcac    - mfaac;
            //            m0    = m2        + mfbac;
            //            mfaac = m0;
            //            m0   += c1o3 * oMdrho;
            //            mfbac = m1 ;
            //            mfcac = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfabc  + mfcbc;
            //            m1    = mfcbc  - mfabc;
            //            m0    = m2        + mfbbc;
            //            mfabc = m0;
            //            mfbbc = m1 ;
            //            mfcbc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m2    = mfacc  + mfccc;
            //            m1    = mfccc  - mfacc;
            //            m0    = m2        + mfbcc;
            //            mfacc = m0;
            //            m0   += c1o9 * oMdrho;
            //            mfbcc = m1 ;
            //            mfccc = m2 ;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //
            //
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // BGK
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            real OxxPyyPzz = omega;
            //            real OxyyPxzz  = omega;//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
            //            real OxyyMxzz  = omega;//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
            //            real O4        = omega;
            //            real O5        = omega;
            //            real O6        = omega;
            //
            //            real mxxPyyPzz = mfcaa + mfaca + mfaac;
            //            real mxxMyy    = mfcaa - mfaca;
            //            real mxxMzz      = mfcaa - mfaac;
            //
            //            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //            //incl. correction
            //            {
            //                real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz+(-two*vx2+vy2+vz2)*rho) + c1o2 * OxxPyyPzz * (mfaaa+(vx2+vy2+vz2)*rho - mxxPyyPzz);
            //                real dyuy = dxux + omega * c3o2 * (mxxMyy+(-vx2+vy2)*rho);
            //                real dzuz = dxux + omega * c3o2 * (mxxMzz+(-vx2+vz2)*rho);
            //
            //                //relax
            //                mxxPyyPzz += OxxPyyPzz*(mfaaa +(vx2+vy2+vz2)*rho - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
            //                mxxMyy    += omega * ((vx2-vy2)*rho-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
            //                mxxMzz    += omega * ((vx2-vz2)*rho-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);
            //            }
            //            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //
            ////             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////             //no correction
            ////             mxxPyyPzz += OxxPyyPzz*(mfaaa+(vx2+vy2+vz2)*rho-mxxPyyPzz);
            ////             mxxMyy    += -(-omega) * ((vx2-vy2)*rho-mxxMyy);
            ////             mxxMzz    += -(-omega) * ((vx2-vz2)*rho-mxxMzz);
            ////             ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //            mfabb     += omega * ((vvy*vvz)*rho-mfabb);
            //            mfbab     += omega * ((vvx*vvz)*rho-mfbab);
            //            mfbba     += omega * ((vvx*vvy)*rho-mfbba);
            //
            //            // linear combinations back
            //            mfcaa = c1o3 * (       mxxMyy +      mxxMzz + mxxPyyPzz);
            //            mfaca = c1o3 * (-two*  mxxMyy +      mxxMzz + mxxPyyPzz);
            //            mfaac = c1o3 * (       mxxMyy - two* mxxMzz + mxxPyyPzz);
            //
            //            //3.
            //            // linear combinations
            //
            //            real mxxyPyzz = mfcba + mfabc;
            //            real mxxyMyzz = mfcba - mfabc;
            //
            //            real mxxzPyyz = mfcab + mfacb;
            //            real mxxzMyyz = mfcab - mfacb;
            //
            //            real mxyyPxzz = mfbca + mfbac;
            //            real mxyyMxzz = mfbca - mfbac;
            //
            //            mxxyMyzz += OxyyMxzz*((vx2-vz2)*vvy*rho-mxxyMyzz);
            //            mxxzMyyz += OxyyMxzz*((vx2-vy2)*vvz*rho-mxxzMyyz);
            //            mxyyMxzz += OxyyMxzz*((vy2-vz2)*vvx*rho-mxyyMxzz);
            //
            //            mxxyPyzz += OxyyPxzz*((c2o3+vx2+vz2)*vvy*rho-mxxyPyzz);
            //            mxxzPyyz += OxyyPxzz*((c2o3+vx2+vy2)*vvz*rho-mxxzPyyz);
            //            mxyyPxzz += OxyyPxzz*((c2o3+vy2+vz2)*vvx*rho-mxyyPxzz);
            //
            //            mfbbb += OxyyMxzz * (vvx*vvy*vvz*rho - mfbbb);
            //            
            //            mfcba = ( mxxyMyzz + mxxyPyzz) * c1o2;
            //            mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
            //            mfcab = ( mxxzMyyz + mxxzPyyz) * c1o2;
            //            mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
            //            mfbca = ( mxyyMxzz + mxyyPxzz) * c1o2;
            //            mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;
            //
            //            //4.
            //            //mfacc += O4*((c1o3+vy2)*(c1o3+vz2)*rho+c1o9*(mfaaa-one)-mfacc);
            //            //mfcac += O4*((c1o3+vx2)*(c1o3+vz2)*rho+c1o9*(mfaaa-one)-mfcac);
            //            //mfcca += O4*((c1o3+vx2)*(c1o3+vy2)*rho+c1o9*(mfaaa-one)-mfcca);
            //            mfacc += O4*((c1o3+vy2)*(c1o3+vz2)*rho-c1o9-mfacc);
            //            mfcac += O4*((c1o3+vx2)*(c1o3+vz2)*rho-c1o9-mfcac);
            //            mfcca += O4*((c1o3+vx2)*(c1o3+vy2)*rho-c1o9-mfcca);
            //            
            //            mfcbb += O4*((c1o3+vx2)*vvy*vvz*rho-mfcbb);
            //            mfbcb += O4*((c1o3+vy2)*vvx*vvz*rho-mfbcb);
            //            mfbbc += O4*((c1o3+vz2)*vvx*vvy*rho-mfbbc);
            //
            //            //5.
            //            mfbcc += O5*((c1o3+vy2)*(c1o3+vz2)*vvx*rho-mfbcc);
            //            mfcbc += O5*((c1o3+vx2)*(c1o3+vz2)*vvy*rho-mfcbc);
            //            mfccb += O5*((c1o3+vx2)*(c1o3+vy2)*vvz*rho-mfccb);
            //
            //            //6.
            //            mfccc += O6*((c1o3+vx2)*(c1o3+vy2)*(c1o3+vz2)*rho-c1o27-mfccc);
            //
            //
            //            //bad fix
            //            vvx = zero;
            //            vvy = zero;
            //            vvz = zero;
            //            vx2 = zero;
            //            vy2 = zero;
            //            vz2 = zero;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            //back
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // Z - Dir
            //            m0 =  mfaac * c1o2 +      mfaab * (vvz - c1o2) + (mfaaa + one* oMdrho) * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfaac        - two* mfaab *  vvz         +  mfaaa                * (one- vz2)              - one* oMdrho * vz2; 
            //            m2 =  mfaac * c1o2 +      mfaab * (vvz + c1o2) + (mfaaa + one* oMdrho) * (     vz2 + vvz) * c1o2;
            //            mfaaa = m0;
            //            mfaab = m1;
            //            mfaac = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfabc * c1o2 +      mfabb * (vvz - c1o2) + mfaba * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfabc        - two* mfabb *  vvz         + mfaba * (one- vz2); 
            //            m2 =  mfabc * c1o2 +      mfabb * (vvz + c1o2) + mfaba * (     vz2 + vvz) * c1o2;
            //            mfaba = m0;
            //            mfabb = m1;
            //            mfabc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfacc * c1o2 +      mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfacc        - two* mfacb *  vvz         +  mfaca                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
            //            m2 =  mfacc * c1o2 +      mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
            //            mfaca = m0;
            //            mfacb = m1;
            //            mfacc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfbac * c1o2 +      mfbab * (vvz - c1o2) + mfbaa * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfbac        - two* mfbab *  vvz         + mfbaa * (one- vz2); 
            //            m2 =  mfbac * c1o2 +      mfbab * (vvz + c1o2) + mfbaa * (     vz2 + vvz) * c1o2;
            //            mfbaa = m0;
            //            mfbab = m1;
            //            mfbac = m2;
            //            /////////b//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfbbc * c1o2 +      mfbbb * (vvz - c1o2) + mfbba * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfbbc        - two* mfbbb *  vvz         + mfbba * (one- vz2); 
            //            m2 =  mfbbc * c1o2 +      mfbbb * (vvz + c1o2) + mfbba * (     vz2 + vvz) * c1o2;
            //            mfbba = m0;
            //            mfbbb = m1;
            //            mfbbc = m2;
            //            /////////b//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfbcc * c1o2 +      mfbcb * (vvz - c1o2) + mfbca * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfbcc        - two* mfbcb *  vvz         + mfbca * (one- vz2); 
            //            m2 =  mfbcc * c1o2 +      mfbcb * (vvz + c1o2) + mfbca * (     vz2 + vvz) * c1o2;
            //            mfbca = m0;
            //            mfbcb = m1;
            //            mfbcc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcac * c1o2 +      mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfcac        - two* mfcab *  vvz         +  mfcaa                  * (one- vz2)              - c1o3 * oMdrho * vz2; 
            //            m2 =  mfcac * c1o2 +      mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (     vz2 + vvz) * c1o2;
            //            mfcaa = m0;
            //            mfcab = m1;
            //            mfcac = m2;
            //            /////////c//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcbc * c1o2 +      mfcbb * (vvz - c1o2) + mfcba * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfcbc        - two* mfcbb *  vvz         + mfcba * (one- vz2); 
            //            m2 =  mfcbc * c1o2 +      mfcbb * (vvz + c1o2) + mfcba * (     vz2 + vvz) * c1o2;
            //            mfcba = m0;
            //            mfcbb = m1;
            //            mfcbc = m2;
            //            /////////c//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfccc * c1o2 +      mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 - vvz) * c1o2; 
            //            m1 = -mfccc        - two* mfccb *  vvz         +  mfcca                  * (one- vz2)              - c1o9 * oMdrho * vz2; 
            //            m2 =  mfccc * c1o2 +      mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (     vz2 + vvz) * c1o2;
            //            mfcca = m0;
            //            mfccb = m1;
            //            mfccc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // Y - Dir
            //            m0 =  mfaca * c1o2 +      mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfaca        - two* mfaba *  vvy         +  mfaaa                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
            //            m2 =  mfaca * c1o2 +      mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
            //            mfaaa = m0;
            //            mfaba = m1;
            //            mfaca = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfacb * c1o2 +      mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfacb        - two* mfabb *  vvy         +  mfaab                  * (one- vy2)              - c2o3 * oMdrho * vy2; 
            //            m2 =  mfacb * c1o2 +      mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (     vy2 + vvy) * c1o2;
            //            mfaab = m0;
            //            mfabb = m1;
            //            mfacb = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfacc * c1o2 +      mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfacc        - two* mfabc *  vvy         +  mfaac                  * (one- vy2)              - c1o6 * oMdrho * vy2; 
            //            m2 =  mfacc * c1o2 +      mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (     vy2 + vvy) * c1o2;
            //            mfaac = m0;
            //            mfabc = m1;
            //            mfacc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfbca * c1o2 +      mfbba * (vvy - c1o2) + mfbaa * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfbca        - two* mfbba *  vvy         + mfbaa * (one- vy2); 
            //            m2 =  mfbca * c1o2 +      mfbba * (vvy + c1o2) + mfbaa * (     vy2 + vvy) * c1o2;
            //            mfbaa = m0;
            //            mfbba = m1;
            //            mfbca = m2;
            //            /////////b//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfbcb * c1o2 +      mfbbb * (vvy - c1o2) + mfbab * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfbcb        - two* mfbbb *  vvy         + mfbab * (one- vy2); 
            //            m2 =  mfbcb * c1o2 +      mfbbb * (vvy + c1o2) + mfbab * (     vy2 + vvy) * c1o2;
            //            mfbab = m0;
            //            mfbbb = m1;
            //            mfbcb = m2;
            //            /////////b//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfbcc * c1o2 +      mfbbc * (vvy - c1o2) + mfbac * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfbcc        - two* mfbbc *  vvy         + mfbac * (one- vy2); 
            //            m2 =  mfbcc * c1o2 +      mfbbc * (vvy + c1o2) + mfbac * (     vy2 + vvy) * c1o2;
            //            mfbac = m0;
            //            mfbbc = m1;
            //            mfbcc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcca * c1o2 +      mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfcca        - two* mfcba *  vvy         +  mfcaa                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
            //            m2 =  mfcca * c1o2 +      mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
            //            mfcaa = m0;
            //            mfcba = m1;
            //            mfcca = m2;
            //            /////////c//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfccb * c1o2 +      mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfccb        - two* mfcbb *  vvy         +  mfcab                  * (one- vy2)              - c2o9 * oMdrho * vy2; 
            //            m2 =  mfccb * c1o2 +      mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (     vy2 + vvy) * c1o2;
            //            mfcab = m0;
            //            mfcbb = m1;
            //            mfccb = m2;
            //            /////////c//////////////////////////////////////////////////////////////////////////
            //            m0 =  mfccc * c1o2 +      mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 - vvy) * c1o2; 
            //            m1 = -mfccc        - two* mfcbc *  vvy         +  mfcac                   * (one- vy2)              - c1o18 * oMdrho * vy2; 
            //            m2 =  mfccc * c1o2 +      mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (     vy2 + vvy) * c1o2;
            //            mfcac = m0;
            //            mfcbc = m1;
            //            mfccc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            // X - Dir
            //            m0 =  mfcaa * c1o2 +      mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfcaa        - two* mfbaa *  vvx         +  mfaaa                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
            //            m2 =  mfcaa * c1o2 +      mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfaaa = m0;
            //            mfbaa = m1;
            //            mfcaa = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcba * c1o2 +      mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfcba        - two* mfbba *  vvx         +  mfaba                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
            //            m2 =  mfcba * c1o2 +      mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfaba = m0;
            //            mfbba = m1;
            //            mfcba = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcca * c1o2 +      mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfcca        - two* mfbca *  vvx         +  mfaca                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
            //            m2 =  mfcca * c1o2 +      mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfaca = m0;
            //            mfbca = m1;
            //            mfcca = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcab * c1o2 +      mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfcab        - two* mfbab *  vvx         +  mfaab                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
            //            m2 =  mfcab * c1o2 +      mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfaab = m0;
            //            mfbab = m1;
            //            mfcab = m2;
            //            ///////////b////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcbb * c1o2 +      mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfcbb        - two* mfbbb *  vvx         +  mfabb                  * (one- vx2)              - c4o9 * oMdrho * vx2; 
            //            m2 =  mfcbb * c1o2 +      mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfabb = m0;
            //            mfbbb = m1;
            //            mfcbb = m2;
            //            ///////////b////////////////////////////////////////////////////////////////////////
            //            m0 =  mfccb * c1o2 +      mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfccb        - two* mfbcb *  vvx         +  mfacb                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
            //            m2 =  mfccb * c1o2 +      mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfacb = m0;
            //            mfbcb = m1;
            //            mfccb = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            ////////////////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcac * c1o2 +      mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfcac        - two* mfbac *  vvx         +  mfaac                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
            //            m2 =  mfcac * c1o2 +      mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfaac = m0;
            //            mfbac = m1;
            //            mfcac = m2;
            //            ///////////c////////////////////////////////////////////////////////////////////////
            //            m0 =  mfcbc * c1o2 +      mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfcbc        - two* mfbbc *  vvx         +  mfabc                  * (one- vx2)              - c1o9 * oMdrho * vx2; 
            //            m2 =  mfcbc * c1o2 +      mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfabc = m0;
            //            mfbbc = m1;
            //            mfcbc = m2;
            //            ///////////c////////////////////////////////////////////////////////////////////////
            //            m0 =  mfccc * c1o2 +      mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 - vvx) * c1o2; 
            //            m1 = -mfccc        - two* mfbcc *  vvx         +  mfacc                   * (one- vx2)              - c1o36 * oMdrho * vx2; 
            //            m2 =  mfccc * c1o2 +      mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (     vx2 + vvx) * c1o2;
            //            mfacc = m0;
            //            mfbcc = m1;
            //            mfccc = m2;
            //            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            (D.f[dP00])[k] = mfabb;//(D.f[ dP00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ dP00   ])[k   ]                                                                     
            (D.f[dM00])[kw] = mfcbb;//(D.f[ dM00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ dM00   ])[kw  ]                                                                   
            (D.f[d0P0])[k] = mfbab;//(D.f[ d0P0   ])[kn   ] = mfbab;// -  c2over27 ;     (D.f[ d0P0   ])[k   ]
            (D.f[d0M0])[ks] = mfbcb;//(D.f[ d0M0   ])[ks   ] = mfbcb;// -  c2over27 ;     (D.f[ d0M0   ])[ks  ]
            (D.f[d00P])[k] = mfbba;//(D.f[ d00P   ])[kt   ] = mfbba;// -  c2over27 ;     (D.f[ d00P   ])[k   ]
            (D.f[d00M])[kb] = mfbbc;//(D.f[ d00M   ])[kb   ] = mfbbc;// -  c2over27 ;     (D.f[ d00M   ])[kb  ]
            (D.f[dPP0])[k] = mfaab;//(D.f[ dPP0  ])[kne  ] = mfaab;// -  c1over54 ;     (D.f[ dPP0  ])[k   ]
            (D.f[dMM0])[ksw] = mfccb;//(D.f[ dMM0  ])[ksw  ] = mfccb;// -  c1over54 ;     (D.f[ dMM0  ])[ksw ]
            (D.f[dPM0])[ks] = mfacb;//(D.f[ dPM0  ])[kse  ] = mfacb;// -  c1over54 ;     (D.f[ dPM0  ])[ks  ]
            (D.f[dMP0])[kw] = mfcab;//(D.f[ dMP0  ])[knw  ] = mfcab;// -  c1over54 ;     (D.f[ dMP0  ])[kw  ]
            (D.f[dP0P])[k] = mfaba;//(D.f[ dP0P  ])[kte  ] = mfaba;// -  c1over54 ;     (D.f[ dP0P  ])[k   ]
            (D.f[dM0M])[kbw] = mfcbc;//(D.f[ dM0M  ])[kbw  ] = mfcbc;// -  c1over54 ;     (D.f[ dM0M  ])[kbw ]
            (D.f[dP0M])[kb] = mfabc;//(D.f[ dP0M  ])[kbe  ] = mfabc;// -  c1over54 ;     (D.f[ dP0M  ])[kb  ]
            (D.f[dM0P])[kw] = mfcba;//(D.f[ dM0P  ])[ktw  ] = mfcba;// -  c1over54 ;     (D.f[ dM0P  ])[kw  ]
            (D.f[d0PP])[k] = mfbaa;//(D.f[ d0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;     (D.f[ d0PP  ])[k   ]
            (D.f[d0MM])[kbs] = mfbcc;//(D.f[ d0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;     (D.f[ d0MM  ])[kbs ]
            (D.f[d0PM])[kb] = mfbac;//(D.f[ d0PM  ])[kbn  ] = mfbac;// -  c1over54 ;     (D.f[ d0PM  ])[kb  ]
            (D.f[d0MP])[ks] = mfbca;//(D.f[ d0MP  ])[kts  ] = mfbca;// -  c1over54 ;     (D.f[ d0MP  ])[ks  ]
            (D.f[d000])[k] = mfbbb;//(D.f[ d000])[kzero] = mfbbb;// -  c8over27 ;     (D.f[ d000])[k   ]
            (D.f[dPPP])[k] = mfaaa;//(D.f[ dPPP ])[ktne ] = mfaaa;// -  c1over216;     (D.f[ dPPP ])[k   ]
            (D.f[dPMP])[ks] = mfaca;//(D.f[ dPMP ])[ktse ] = mfaca;// -  c1over216;     (D.f[ dPMP ])[ks  ]
            (D.f[dPPM])[kb] = mfaac;//(D.f[ dPPM ])[kbne ] = mfaac;// -  c1over216;     (D.f[ dPPM ])[kb  ]
            (D.f[dPMM])[kbs] = mfacc;//(D.f[ dPMM ])[kbse ] = mfacc;// -  c1over216;     (D.f[ dPMM ])[kbs ]
            (D.f[dMPP])[kw] = mfcaa;//(D.f[ dMPP ])[ktnw ] = mfcaa;// -  c1over216;     (D.f[ dMPP ])[kw  ]
            (D.f[dMMP])[ksw] = mfcca;//(D.f[ dMMP ])[ktsw ] = mfcca;// -  c1over216;     (D.f[ dMMP ])[ksw ]
            (D.f[dMPM])[kbw] = mfcac;//(D.f[ dMPM ])[kbnw ] = mfcac;// -  c1over216;     (D.f[ dMPM ])[kbw ]
            (D.f[dMMM])[kbsw] = mfccc;//(D.f[ dMMM ])[kbsw ] = mfccc;// -  c1over216;     (D.f[ dMMM ])[kbsw]
                                        ////////////////////////////////////////////////////////////////////////////////////
        }
    }
}