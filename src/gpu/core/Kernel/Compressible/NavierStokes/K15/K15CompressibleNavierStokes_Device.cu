//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_Kernel Kernel
//! \ingroup gpu_core core
//! \{
#include "Calculation/Calculation.h" 
#include "lbm/constants/D3Q27.h"
#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
#include "math.h"

__global__ void K15CompressibleNavierStokes_Device(real omega,
    unsigned int* bcMatD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* DDStart,
    int size_Mat,
    int level,
    real* forces,
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

    if (k<size_Mat)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if (BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/)
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

            //unsigned int kzero= k;
            //unsigned int ke   = k;
            //unsigned int kw   = neighborX[k];
            //unsigned int kn   = k;
            //unsigned int ks   = neighborY[k];
            //unsigned int kt   = k;
            //unsigned int kb   = neighborZ[k];
            //unsigned int ksw  = neighborY[kw];
            //unsigned int kne  = k;
            //unsigned int kse  = ks;
            //unsigned int knw  = kw;
            //unsigned int kbw  = neighborZ[kw];
            //unsigned int kte  = k;
            //unsigned int kbe  = kb;
            //unsigned int ktw  = kw;
            //unsigned int kbs  = neighborZ[ks];
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
            //unsigned int kbsw = neighborZ[ksw];
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
            real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

            real rho = c1o1 + drho;
            ////////////////////////////////////////////////////////////////////////////////////
            //slow
            //real oMdrho = one - ((((mfccc+mfaaa) + (mfaca+mfcac)) + ((mfacc+mfcaa) + (mfaac+mfcca))) + 
            //                       (((mfbac+mfbca) + (mfbaa+mfbcc)) + ((mfabc+mfcba) + (mfaba+mfcbc)) + ((mfacb+mfcab) + (mfaab+mfccb))) +
            //                        ((mfabb+mfcbb) + (mfbab+mfbcb)  +  (mfbba+mfbbc)));//fehlt mfbbb
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
            //the force be with you
            real fx = forces[0] / (pow((double)c2o1, (double)level)); //zero;//0.0032653/(pow(two,level)); //0.000000005;//(two/1600000.0) / 120.0; //
            real fy = forces[1] / (pow((double)c2o1, (double)level)); //zero;
            real fz = forces[2] / (pow((double)c2o1, (double)level)); //zero;
            vvx += fx*c1o2;
            vvy += fy*c1o2;
            vvz += fz*c1o2;
            ////////////////////////////////////////////////////////////////////////////////////
            //real omega = omega_in;
            ////////////////////////////////////////////////////////////////////////////////////
            //fast
            real oMdrho = c1o1; // comp special
                               //real oMdrho = one - (mfccc+mfaaa + mfaca+mfcac + mfacc+mfcaa + mfaac+mfcca + 
                               //                       mfbac+mfbca + mfbaa+mfbcc + mfabc+mfcba + mfaba+mfcbc + mfacb+mfcab + mfaab+mfccb +
                               //                       mfabb+mfcbb + mfbab+mfbcb + mfbba+mfbbc + mfbbb);//fehlt mfbbb nicht mehr
                               //real vvx    =mfccc-mfaaa + mfcac-mfaca + mfcaa-mfacc + mfcca-mfaac + 
                               //                mfcba-mfabc + mfcbc-mfaba + mfcab-mfacb + mfccb-mfaab +
                               //                mfcbb-mfabb;
                               //real vvy    =mfccc-mfaaa + mfaca-mfcac + mfacc-mfcaa + mfcca-mfaac + 
                               //                mfbca-mfbac + mfbcc-mfbaa + mfacb-mfcab + mfccb-mfaab +
                               //                mfbcb-mfbab;
                               //real vvz    =mfccc-mfaaa + mfcac-mfaca + mfacc-mfcaa + mfaac-mfcca + 
                               //                mfbac-mfbca + mfbcc-mfbaa + mfabc-mfcba + mfcbc-mfaba +
                               //                mfbbc-mfbba;
                               ////////////////////////////////////////////////////////////////////////////////////
                               // oMdrho assembler style -------> faaaaaastaaaa
                               // or much sloooowaaaa ... it dep�ndssssss on sadaku
            real m0, m1, m2;
            //real oMdrho;
            //{
            //    oMdrho=mfccc+mfaaa;
            //    m0=mfaca+mfcac;
            //    m1=mfacc+mfcaa;
            //    m2=mfaac+mfcca;
            //    oMdrho+=m0;
            //    m1+=m2;
            //    oMdrho+=m1;
            //    m0=mfbac+mfbca;
            //    m1=mfbaa+mfbcc;
            //    m0+=m1;
            //    m1=mfabc+mfcba;
            //    m2=mfaba+mfcbc;
            //    m1+=m2;
            //    m0+=m1;
            //    m1=mfacb+mfcab;
            //    m2=mfaab+mfccb;
            //    m1+=m2;
            //    m0+=m1;
            //    oMdrho+=m0;
            //    m0=mfabb+mfcbb;
            //    m1=mfbab+mfbcb;
            //    m2=mfbba+mfbbc;
            //    m0+=m1+m2;
            //    m0+=mfbbb; //hat gefehlt
            //    oMdrho = one - (oMdrho + m0);
            //}
            //real vvx;
            real vx2;
            //{
            //    vvx = mfccc-mfaaa;
            //    m0  = mfcac-mfaca;
            //    m1  = mfcaa-mfacc;
            //    m2  = mfcca-mfaac;
            //    vvx+= m0;
            //    m1 += m2;
            //    vvx+= m1;
            //    vx2 = mfcba-mfabc;
            //    m0  = mfcbc-mfaba;
            //    m1  = mfcab-mfacb;
            //    m2  = mfccb-mfaab;
            //    vx2+= m0;
            //    m1 += m2;
            //    vx2+= m1;
            //    vvx+= vx2;
            //    vx2 = mfcbb-mfabb;
            //    vvx+= vx2;
            //}
            //real vvy;
            real vy2;
            //{
            //    vvy = mfccc-mfaaa;
            //    m0  = mfaca-mfcac;
            //    m1  = mfacc-mfcaa;
            //    m2  = mfcca-mfaac;
            //    vvy+= m0;
            //    m1 += m2;
            //    vvy+= m1;
            //    vy2 = mfbca-mfbac;
            //    m0  = mfbcc-mfbaa;
            //    m1  = mfacb-mfcab;
            //    m2  = mfccb-mfaab;
            //    vy2+= m0;
            //    m1 += m2;
            //    vy2+= m1;
            //    vvy+= vy2;
            //    vy2 = mfbcb-mfbab;
            //    vvy+= vy2;
            //}
            //real vvz;
            real vz2;
            //{
            //    vvz = mfccc-mfaaa;
            //    m0  = mfcac-mfaca;
            //    m1  = mfacc-mfcaa;
            //    m2  = mfaac-mfcca;
            //    vvz+= m0;
            //    m1 += m2;
            //    vvz+= m1;
            //    vz2 = mfbac-mfbca;
            //    m0  = mfbcc-mfbaa;
            //    m1  = mfabc-mfcba;
            //    m2  = mfcbc-mfaba;
            //    vz2+= m0;
            //    m1 += m2;
            //    vz2+= m1;
            //    vvz+= vz2;
            //    vz2 = mfbbc-mfbba;
            //    vvz+= vz2;
            //}
            vx2 = vvx*vvx;
            vy2 = vvy*vvy;
            vz2 = vvz*vvz;
            //////////////////////////////////////////////////////////////////////////////////////
            //// test rundungsfehler....
            //mfcbb = (mfcbb - c2over27* (rho-one))/(rho);
            //mfabb = (mfabb - c2over27* (rho-one))/(rho);
            //mfbcb = (mfbcb - c2over27* (rho-one))/(rho);
            //mfbab = (mfbab - c2over27* (rho-one))/(rho);
            //mfbbc = (mfbbc - c2over27* (rho-one))/(rho);
            //mfbba = (mfbba - c2over27* (rho-one))/(rho);
            //mfccb = (mfccb - c1over54* (rho-one))/(rho);
            //mfaab = (mfaab - c1over54* (rho-one))/(rho);
            //mfcab = (mfcab - c1over54* (rho-one))/(rho);
            //mfacb = (mfacb - c1over54* (rho-one))/(rho);
            //mfcbc = (mfcbc - c1over54* (rho-one))/(rho);
            //mfaba = (mfaba - c1over54* (rho-one))/(rho);
            //mfcba = (mfcba - c1over54* (rho-one))/(rho);
            //mfabc = (mfabc - c1over54* (rho-one))/(rho);
            //mfbcc = (mfbcc - c1over54* (rho-one))/(rho);
            //mfbaa = (mfbaa - c1over54* (rho-one))/(rho);
            //mfbca = (mfbca - c1over54* (rho-one))/(rho);
            //mfbac = (mfbac - c1over54* (rho-one))/(rho);
            //mfbbb = (mfbbb - c8over27* (rho-one))/(rho);
            //mfccc = (mfccc - c1over216*(rho-one))/(rho);
            //mfaac = (mfaac - c1over216*(rho-one))/(rho);
            //mfcac = (mfcac - c1over216*(rho-one))/(rho);
            //mfacc = (mfacc - c1over216*(rho-one))/(rho);
            //mfcca = (mfcca - c1over216*(rho-one))/(rho);
            //mfaaa = (mfaaa - c1over216*(rho-one))/(rho);
            //mfcaa = (mfcaa - c1over216*(rho-one))/(rho);
            //mfaca = (mfaca - c1over216*(rho-one))/(rho);
            ////////////////////////////////////////////////////////////////////////////////////
            real wadjust;
            real qudricLimitP = c1o100;// * 0.0001f;
            real qudricLimitM = c1o100;// * 0.0001f;
            real qudricLimitD = c1o100;// * 0.001f;
                                       //real s9 = minusomega;
                                       //test
                                       //s9 = 0.;
                                       ////////////////////////////////////////////////////////////////////////////////////
                                       //Hin
                                       ////////////////////////////////////////////////////////////////////////////////////
                                       // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
                                       ////////////////////////////////////////////////////////////////////////////////////
                                       // Z - Dir
            m2 = mfaaa + mfaac;
            m1 = mfaac - mfaaa;
            m0 = m2 + mfaab;
            mfaaa = m0;
            m0 += c1o36 * oMdrho;
            mfaab = m1 - m0 * vvz;
            mfaac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaba + mfabc;
            m1 = mfabc - mfaba;
            m0 = m2 + mfabb;
            mfaba = m0;
            m0 += c1o9 * oMdrho;
            mfabb = m1 - m0 * vvz;
            mfabc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaca + mfacc;
            m1 = mfacc - mfaca;
            m0 = m2 + mfacb;
            mfaca = m0;
            m0 += c1o36 * oMdrho;
            mfacb = m1 - m0 * vvz;
            mfacc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfbaa + mfbac;
            m1 = mfbac - mfbaa;
            m0 = m2 + mfbab;
            mfbaa = m0;
            m0 += c1o9 * oMdrho;
            mfbab = m1 - m0 * vvz;
            mfbac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfbba + mfbbc;
            m1 = mfbbc - mfbba;
            m0 = m2 + mfbbb;
            mfbba = m0;
            m0 += c4o9 * oMdrho;
            mfbbb = m1 - m0 * vvz;
            mfbbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfbca + mfbcc;
            m1 = mfbcc - mfbca;
            m0 = m2 + mfbcb;
            mfbca = m0;
            m0 += c1o9 * oMdrho;
            mfbcb = m1 - m0 * vvz;
            mfbcc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfcaa + mfcac;
            m1 = mfcac - mfcaa;
            m0 = m2 + mfcab;
            mfcaa = m0;
            m0 += c1o36 * oMdrho;
            mfcab = m1 - m0 * vvz;
            mfcac = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfcba + mfcbc;
            m1 = mfcbc - mfcba;
            m0 = m2 + mfcbb;
            mfcba = m0;
            m0 += c1o9 * oMdrho;
            mfcbb = m1 - m0 * vvz;
            mfcbc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfcca + mfccc;
            m1 = mfccc - mfcca;
            m0 = m2 + mfccb;
            mfcca = m0;
            m0 += c1o36 * oMdrho;
            mfccb = m1 - m0 * vvz;
            mfccc = m2 - c2o1*    m1 * vvz + vz2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            m2 = mfaaa + mfaca;
            m1 = mfaca - mfaaa;
            m0 = m2 + mfaba;
            mfaaa = m0;
            m0 += c1o6 * oMdrho;
            mfaba = m1 - m0 * vvy;
            mfaca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaab + mfacb;
            m1 = mfacb - mfaab;
            m0 = m2 + mfabb;
            mfaab = m0;
            mfabb = m1 - m0 * vvy;
            mfacb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaac + mfacc;
            m1 = mfacc - mfaac;
            m0 = m2 + mfabc;
            mfaac = m0;
            m0 += c1o18 * oMdrho;
            mfabc = m1 - m0 * vvy;
            mfacc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfbaa + mfbca;
            m1 = mfbca - mfbaa;
            m0 = m2 + mfbba;
            mfbaa = m0;
            m0 += c2o3 * oMdrho;
            mfbba = m1 - m0 * vvy;
            mfbca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfbab + mfbcb;
            m1 = mfbcb - mfbab;
            m0 = m2 + mfbbb;
            mfbab = m0;
            mfbbb = m1 - m0 * vvy;
            mfbcb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfbac + mfbcc;
            m1 = mfbcc - mfbac;
            m0 = m2 + mfbbc;
            mfbac = m0;
            m0 += c2o9 * oMdrho;
            mfbbc = m1 - m0 * vvy;
            mfbcc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfcaa + mfcca;
            m1 = mfcca - mfcaa;
            m0 = m2 + mfcba;
            mfcaa = m0;
            m0 += c1o6 * oMdrho;
            mfcba = m1 - m0 * vvy;
            mfcca = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfcab + mfccb;
            m1 = mfccb - mfcab;
            m0 = m2 + mfcbb;
            mfcab = m0;
            mfcbb = m1 - m0 * vvy;
            mfccb = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfcac + mfccc;
            m1 = mfccc - mfcac;
            m0 = m2 + mfcbc;
            mfcac = m0;
            m0 += c1o18 * oMdrho;
            mfcbc = m1 - m0 * vvy;
            mfccc = m2 - c2o1*    m1 * vvy + vy2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9        Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // X - Dir
            m2 = mfaaa + mfcaa;
            m1 = mfcaa - mfaaa;
            m0 = m2 + mfbaa;
            mfaaa = m0;
            m0 += c1o1* oMdrho;
            mfbaa = m1 - m0 * vvx;
            mfcaa = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaba + mfcba;
            m1 = mfcba - mfaba;
            m0 = m2 + mfbba;
            mfaba = m0;
            mfbba = m1 - m0 * vvx;
            mfcba = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaca + mfcca;
            m1 = mfcca - mfaca;
            m0 = m2 + mfbca;
            mfaca = m0;
            m0 += c1o3 * oMdrho;
            mfbca = m1 - m0 * vvx;
            mfcca = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaab + mfcab;
            m1 = mfcab - mfaab;
            m0 = m2 + mfbab;
            mfaab = m0;
            mfbab = m1 - m0 * vvx;
            mfcab = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfabb + mfcbb;
            m1 = mfcbb - mfabb;
            m0 = m2 + mfbbb;
            mfabb = m0;
            mfbbb = m1 - m0 * vvx;
            mfcbb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfacb + mfccb;
            m1 = mfccb - mfacb;
            m0 = m2 + mfbcb;
            mfacb = m0;
            mfbcb = m1 - m0 * vvx;
            mfccb = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfaac + mfcac;
            m1 = mfcac - mfaac;
            m0 = m2 + mfbac;
            mfaac = m0;
            m0 += c1o3 * oMdrho;
            mfbac = m1 - m0 * vvx;
            mfcac = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfabc + mfcbc;
            m1 = mfcbc - mfabc;
            m0 = m2 + mfbbc;
            mfabc = m0;
            mfbbc = m1 - m0 * vvx;
            mfcbc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            m2 = mfacc + mfccc;
            m1 = mfccc - mfacc;
            m0 = m2 + mfbcc;
            mfacc = m0;
            m0 += c1o9 * oMdrho;
            mfbcc = m1 - m0 * vvx;
            mfccc = m2 - c2o1*    m1 * vvx + vx2 * m0;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            real OxxPyyPzz = c1o1; //omega; // one;    //set the bulk viscosity one is high / two is very low and zero is (too) high

                                  ////////////////////////////////////////////////////////////
                                  //3.
                                  //////////////////////////////
            real OxyyPxzz = c1o1;//three  * (two - omega) / (three  - omega);//one;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//two-omega;//eight*(two-omega)/(eight -omega);//one;//omega;//two-omega;//
            real OxyyMxzz = c1o1;//six    * (two - omega) / (six    - omega);//one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
            real Oxyz = c1o1;//twelve * (two - omega) / (twelve + omega);//one;//two-omega;//(1000.*(-2. + omega))/(-1000. + 439.*omega);//(eight * (omega - two)) / (omega - eight);//omega;//one;//eight*(two-omega)/(eight -omega);//one;//two-omega;//one;// 
                            //////////////////////////////
                            //real OxyyPxzz  = two-omega;//
                            //real OxyyMxzz  = two-omega;//
                            //////////////////////////////
                            //real OxyyPxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
                            //real OxyyMxzz  = (eight * (omega - two)) / (omega - eight);//Ginzburg
                            //////////////////////////////
                            //real OxyyPxzz  = omega;//BGK
                            //real OxyyMxzz  = omega;//BGK
                            //////////////////////////////
                            //real OxyyPxzz  = (one + omega) / two;//1P5
                            //real OxyyMxzz  = (one + omega) / two;//1P5
                            //////////////////////////////
                            //real OxyyPxzz  = (three - omega) / two;//0P5
                            //real OxyyMxzz  = (three - omega) / two;//0P5
                            //////////////////////////////
                            //real OxyyPxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
                            //real OxyyMxzz  = (one + (eight * (omega - two)) / (omega - eight)) / two;//one + Ginzburg / two ... Car
                            ////////////////////////////////////////////////////////////
                            //4.
                            //////////////////////////////
            real O4 = c1o1;
            //////////////////////////////
            //real O4        = omega;//TRT
            ////////////////////////////////////////////////////////////
            //5.
            //////////////////////////////
            real O5 = c1o1;
            ////////////////////////////////////////////////////////////
            //6.
            //////////////////////////////
            real O6 = c1o1;
            ////////////////////////////////////////////////////////////


            //central moments to cumulants
            //4.
            //real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab) / rho;  //bis 15.05.2015 verwendet
            //real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb) / rho;  //bis 15.05.2015 verwendet
            //real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb) / rho;  //bis 15.05.2015 verwendet
            real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;    //ab 15.05.2015 verwendet
            real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho; //ab 15.05.2015 verwendet
            real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho; //ab 15.05.2015 verwendet

            real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
            real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
            real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

            //5.
            real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
            real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
            real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

            //6.

            real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb
                - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                - c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                - c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
                + (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                    + c2o1 * (mfcaa * mfaca * mfaac)
                    + c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
                - c1o3 * (mfacc + mfcac + mfcca) / rho
                - c1o9 * (mfcaa + mfaca + mfaac) / rho
                + (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                    + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
                + c1o27*((drho * drho - drho) / (rho*rho)));
            //+ c1o27*(one -three/rho +two/(rho*rho)));




            ////Cum 4.
            //real CUMcbb = mfcbb - ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab)/rho; // 
            //real CUMbcb = mfbcb - ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb)/rho; // 
            //real CUMbbc = mfbbc - ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb)/rho; // 

            //real CUMcca = mfcca - ((mfcaa * mfaca + two * mfbba * mfbba) /rho + c1o3 * (mfcaa + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho);
            //real CUMcac = mfcac - ((mfcaa * mfaac + two * mfbab * mfbab) /rho + c1o3 * (mfcaa + mfaac) * oMdrho + c1o9*(oMdrho-one)*oMdrho);
            //real CUMacc = mfacc - ((mfaac * mfaca + two * mfabb * mfabb) /rho + c1o3 * (mfaac + mfaca) * oMdrho + c1o9*(oMdrho-one)*oMdrho);

            ////Cum 5.
            //real CUMbcc = mfbcc - (mfaac * mfbca + mfaca * mfbac + four * mfabb * mfbbb + two * (mfbab * mfacb + mfbba * mfabc)) / rho - c1o3 * (mfbca + mfbac) * oMdrho;
            //real CUMcbc = mfcbc - (mfaac * mfcba + mfcaa * mfabc + four * mfbab * mfbbb + two * (mfabb * mfcab + mfbba * mfbac)) / rho - c1o3 * (mfcba + mfabc) * oMdrho;
            //real CUMccb = mfccb - (mfcaa * mfacb + mfaca * mfcab + four * mfbba * mfbbb + two * (mfbab * mfbca + mfabb * mfcba)) / rho - c1o3 * (mfacb + mfcab) * oMdrho;

            ////Cum 6.
            //real CUMccc = mfccc  +((-four *  mfbbb * mfbbb  
            //                -           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
            //                -    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
            //                -     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
            //                +(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
            //                +     two * (mfcaa * mfaca * mfaac)
            //                + sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
            //                -    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
            //                -    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
            //                +(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
            //                +           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) / rho * c2o3*oMdrho) +c1o27*oMdrho;





            //2.
            // linear combinations
            real mxxPyyPzz = mfcaa + mfaca + mfaac;
            real mxxMyy = mfcaa - mfaca;
            real mxxMzz = mfcaa - mfaac;

            //////////////////////////////////////////////////////////////////////////
            //             real magicBulk=(CUMacc+CUMcac+CUMcca)*(one/OxxPyyPzz-c1o2)*c3o2*8.;

            //////////////////////////////////////////////////////////////////////////
            //limiter-Scheise Teil 1
            //real oxxyy,oxxzz,oxy,oxz,oyz;
            //real smag=0.001;
            //oxxyy    = omega+(one-omega)*abs(mxxMyy)/(abs(mxxMyy)+smag);
            //oxxzz    = omega+(one-omega)*abs(mxxMzz)/(abs(mxxMzz)+smag);
            //oxy      = omega+(one-omega)*abs(mfbba)/(abs(mfbba)+smag);
            //oxz      = omega+(one-omega)*abs(mfbab)/(abs(mfbab)+smag);
            //oyz      = omega+(one-omega)*abs(mfabb)/(abs(mfabb)+smag);

            ////////////////////////////////////////////////////////////////////////////
            ////Teil 1b
            //real constante = 1000.0;
            //real nuEddi = constante * abs(mxxPyyPzz);
            //real omegaLimit = one / (one / omega + three * nuEddi);

            //{
            //    real dxux = c1o2 * (-omegaLimit) *(mxxMyy + mxxMzz) +  OxxPyyPzz * (mfaaa - mxxPyyPzz);
            //    real dyuy = dxux + omegaLimit * c3o2 * mxxMyy;
            //    real dzuz = dxux + omegaLimit * c3o2 * mxxMzz;

            ////relax
            //mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- three * (one - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);
            //mxxMyy    += omegaLimit * (-mxxMyy) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vy2 * dyuy);
            //mxxMzz    += omegaLimit * (-mxxMzz) - three * (one + c1o2 * (-omegaLimit)) * (vx2 * dxux + vz2 * dzuz);

            //}
            //mfabb     += omegaLimit * (-mfabb);
            //mfbab     += omegaLimit * (-mfbab);
            //mfbba     += omegaLimit * (-mfbba);
            ////////////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //incl. correction        (hat noch nicht so gut funktioniert...Optimierungsbedarf??)
            {
                real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
                real dyuy = dxux + omega * c3o2 * mxxMyy;
                real dzuz = dxux + omega * c3o2 * mxxMzz;

                //relax
                mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
                mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
                mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

                //////////////////////////////////////////////////////////////////////////
                //limiter-Scheise Teil 2
                //mxxMyy    += oxxyy * (-mxxMyy) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vy2 * dyuy);
                //mxxMzz    += oxxzz * (-mxxMzz) - three * (one + c1o2 * (-omega)) * (vx2 * dxux + vz2 * dzuz);
                //////////////////////////////////////////////////////////////////////////

            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////no correction
            //mxxPyyPzz += OxxPyyPzz*(mfaaa-mxxPyyPzz);//-magicBulk*OxxPyyPzz;
            //mxxMyy    += -(-omega) * (-mxxMyy);
            //mxxMzz    += -(-omega) * (-mxxMzz);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            mfabb += omega * (-mfabb);
            mfbab += omega * (-mfbab);
            mfbba += omega * (-mfbba);

            //////////////////////////////////////////////////////////////////////////
            //limiter-Scheise Teil 3
            //mfabb     += oyz * (-mfabb);
            //mfbab     += oxz * (-mfbab);
            //mfbba     += oxy * (-mfbba);
            //////////////////////////////////////////////////////////////////////////

            // linear combinations back
            mfcaa = c1o3 * (mxxMyy + mxxMzz + mxxPyyPzz);
            mfaca = c1o3 * (-c2o1*  mxxMyy + mxxMzz + mxxPyyPzz);
            mfaac = c1o3 * (mxxMyy - c2o1* mxxMzz + mxxPyyPzz);

            //3.
            // linear combinations

            real mxxyPyzz = mfcba + mfabc;
            real mxxyMyzz = mfcba - mfabc;

            real mxxzPyyz = mfcab + mfacb;
            real mxxzMyyz = mfcab - mfacb;

            real mxyyPxzz = mfbca + mfbac;
            real mxyyMxzz = mfbca - mfbac;

            //relax
            //////////////////////////////////////////////////////////////////////////
            //das ist der limiter
            wadjust = Oxyz + (c1o1 - Oxyz)*abs(mfbbb) / (abs(mfbbb) + qudricLimitD);
            mfbbb += wadjust * (-mfbbb);
            wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxxyPyzz) / (abs(mxxyPyzz) + qudricLimitP);
            mxxyPyzz += wadjust * (-mxxyPyzz);
            wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxxyMyzz) / (abs(mxxyMyzz) + qudricLimitM);
            mxxyMyzz += wadjust * (-mxxyMyzz);
            wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxxzPyyz) / (abs(mxxzPyyz) + qudricLimitP);
            mxxzPyyz += wadjust * (-mxxzPyyz);
            wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxxzMyyz) / (abs(mxxzMyyz) + qudricLimitM);
            mxxzMyyz += wadjust * (-mxxzMyyz);
            wadjust = OxyyPxzz + (c1o1 - OxyyPxzz)*abs(mxyyPxzz) / (abs(mxyyPxzz) + qudricLimitP);
            mxyyPxzz += wadjust * (-mxyyPxzz);
            wadjust = OxyyMxzz + (c1o1 - OxyyMxzz)*abs(mxyyMxzz) / (abs(mxyyMxzz) + qudricLimitM);
            mxyyMxzz += wadjust * (-mxyyMxzz);
            //////////////////////////////////////////////////////////////////////////
            //ohne limiter
            //mfbbb     += OxyyMxzz * (-mfbbb);
            //mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
            //mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
            //mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
            //mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
            //mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
            //mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
            //////////////////////////////////////////////////////////////////////////

            //// linear combinations back
            ////generic
            //mfcba =  zero;
            //mfabc =  zero;
            //mfcab =  zero;
            //mfacb =  zero;
            //mfbca =  zero;
            //mfbac =  zero;

            mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
            mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
            mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
            mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
            mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
            mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

            //4.
            //CUMacc =  zero; 
            //CUMcac =  zero; 
            //CUMcca =  zero; 
            //       
            //CUMbbc =  zero; 
            //CUMbcb =  zero; 
            //CUMcbb =  zero; 
            //      
            ////5.   
            //CUMbcc =  zero;
            //CUMcbc =  zero;
            //CUMccb =  zero;
            //       
            ////6.   
            //CUMccc =  zero;

            //4.
            //////////////////////////////////////////////////////////////////////////
            //mit limiter
            //    wadjust    = O4+(one-O4)*abs(CUMacc)/(abs(CUMacc)+qudricLimit);
            //CUMacc    += wadjust * (-CUMacc);
            //    wadjust    = O4+(one-O4)*abs(CUMcac)/(abs(CUMcac)+qudricLimit);
            //CUMcac    += wadjust * (-CUMcac); 
            //    wadjust    = O4+(one-O4)*abs(CUMcca)/(abs(CUMcca)+qudricLimit);
            //CUMcca    += wadjust * (-CUMcca); 

            //    wadjust    = O4+(one-O4)*abs(CUMbbc)/(abs(CUMbbc)+qudricLimit);
            //CUMbbc    += wadjust * (-CUMbbc); 
            //    wadjust    = O4+(one-O4)*abs(CUMbcb)/(abs(CUMbcb)+qudricLimit);
            //CUMbcb    += wadjust * (-CUMbcb); 
            //    wadjust    = O4+(one-O4)*abs(CUMcbb)/(abs(CUMcbb)+qudricLimit);
            //CUMcbb    += wadjust * (-CUMcbb); 
            //////////////////////////////////////////////////////////////////////////
            //ohne limiter
            CUMacc += O4 * (-CUMacc);
            CUMcac += O4 * (-CUMcac);
            CUMcca += O4 * (-CUMcca);

            CUMbbc += O4 * (-CUMbbc);
            CUMbcb += O4 * (-CUMbcb);
            CUMcbb += O4 * (-CUMcbb);
            //////////////////////////////////////////////////////////////////////////


            //5.
            CUMbcc += O5 * (-CUMbcc);
            CUMcbc += O5 * (-CUMcbc);
            CUMccb += O5 * (-CUMccb);

            //6.
            CUMccc += O6 * (-CUMccc);



            //back cumulants to central moments
            //4.
            //mfcbb = CUMcbb + ((mfcaa + c1o3 * oMdrho) * mfabb + two * mfbba * mfbab) / rho; //bis 15.05.2015 verwendet
            //mfbcb = CUMbcb + ((mfaca + c1o3 * oMdrho) * mfbab + two * mfbba * mfabb) / rho; //bis 15.05.2015 verwendet
            //mfbbc = CUMbbc + ((mfaac + c1o3 * oMdrho) * mfbba + two * mfbab * mfabb) / rho; //bis 15.05.2015 verwendet
            mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho; //ab 15.05.2015 verwendet
            mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho; //ab 15.05.2015 verwendet
            mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho; //ab 15.05.2015 verwendet

            mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));//(one/rho-one));
            mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));//(one/rho-one));
            mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));//(one/rho-one));

                                                                                                                          //5.
            mfbcc = CUMbcc + ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) / rho;
            mfcbc = CUMcbc + ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) / rho;
            mfccb = CUMccb + ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) / rho;

            //6.

            mfccc = CUMccc - ((-c4o1 *  mfbbb * mfbbb
                - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                - c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                - c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
                + (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                    + c2o1 * (mfcaa * mfaca * mfaac)
                    + c16o1 *  mfbba * mfbab * mfabb) / (rho * rho)
                - c1o3 * (mfacc + mfcac + mfcca) / rho
                - c1o9 * (mfcaa + mfaca + mfaac) / rho
                + (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                    + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) / (rho * rho) * c2o3
                + c1o27*((drho * drho - drho) / (rho*rho)));
            //+ c1o27*(one -three/rho +two/(rho*rho)));


            //mfccc = CUMccc - ((-four *  mfbbb * mfbbb  
            //                -           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
            //                -    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
            //                -     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
            //                +(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
            //                +     two * (mfcaa * mfaca * mfaac)
            //                + sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
            //                -    c1o3 * (mfacc + mfcac + mfcca) * oMdrho  -c1o9*oMdrho*oMdrho
            //                -    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho)
            //                +(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
            //                +           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) / rho * c2o3*oMdrho) -c1o27*oMdrho;





            //mfccc = CUMccc  -(( -four *  mfbbb * mfbbb  
            //                -           (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
            //                -    four * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
            //                -     two * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) / rho
            //                +(   four * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
            //                +     two * (mfcaa * mfaca * mfaac)
            //                + sixteen *  mfbba * mfbab * mfabb) / (rho * rho)
            //                +    (c1o3 * (mfacc + mfcac + mfcca) )/rho
            //                +     (c1o9 * (mfcaa + mfaca + mfaac) )/rho
            //                - (four*c1o3*(mfabb * mfabb + mfbab * mfbab + mfbba * mfbba) + two * ( c1o3 *(mfcaa * mfaca + mfcaa * mfaac + mfaca * mfaac) + c1o9* (mfcaa + mfaca + mfaac )))/(rho*rho) 
            //                + c1o27*(three/rho - two/(rho * rho) - one);
            //                //-    c1o9 * (mfcaa + mfaca + mfaac) * oMdrho*(one-two* oMdrho)- c1o27* oMdrho * oMdrho*(-two* oMdrho) //????
            //                //+(    two * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
            //                //+           (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa)) / rho * c2o3*oMdrho) -c1o27*oMdrho;

            //Test Cascade
            //4.
            //mfacc = mfaaa * c1o9 * O4 + (one - O4) * mfacc; 
            //mfcac = mfaaa * c1o9 * O4 + (one - O4) * mfcac; 
            //mfcca = mfaaa * c1o9 * O4 + (one - O4) * mfcca; 
            //
            //mfbbc += O4 * (-mfbbc); 
            //mfbcb += O4 * (-mfbcb); 
            //mfcbb += O4 * (-mfcbb); 

            ////5.
            //mfbcc += O5 * (-mfbcc);
            //mfcbc += O5 * (-mfcbc);
            //mfccb += O5 * (-mfccb);

            //6.
            //mfccc = mfaaa * c1o27 * O6 + (one - O6) * mfccc;

            ////////////////////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////////////////////
            //the force be with you
            mfbaa = -mfbaa;
            mfaba = -mfaba;
            mfaab = -mfaab;
            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            //back
            ////////////////////////////////////////////////////////////////////////////////////
            //mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Z - Dir
            m0 = mfaac * c1o2 + mfaab * (vvz - c1o2) + (mfaaa + c1o1* oMdrho) * (vz2 - vvz) * c1o2;
            m1 = -mfaac - c2o1* mfaab *  vvz + mfaaa                * (c1o1 - vz2) - c1o1* oMdrho * vz2;
            m2 = mfaac * c1o2 + mfaab * (vvz + c1o2) + (mfaaa + c1o1* oMdrho) * (vz2 + vvz) * c1o2;
            mfaaa = m0;
            mfaab = m1;
            mfaac = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfabc * c1o2 + mfabb * (vvz - c1o2) + mfaba * (vz2 - vvz) * c1o2;
            m1 = -mfabc - c2o1* mfabb *  vvz + mfaba * (c1o1 - vz2);
            m2 = mfabc * c1o2 + mfabb * (vvz + c1o2) + mfaba * (vz2 + vvz) * c1o2;
            mfaba = m0;
            mfabb = m1;
            mfabc = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfacc * c1o2 + mfacb * (vvz - c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
            m1 = -mfacc - c2o1* mfacb *  vvz + mfaca                  * (c1o1 - vz2) - c1o3 * oMdrho * vz2;
            m2 = mfacc * c1o2 + mfacb * (vvz + c1o2) + (mfaca + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
            mfaca = m0;
            mfacb = m1;
            mfacc = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfbac * c1o2 + mfbab * (vvz - c1o2) + mfbaa * (vz2 - vvz) * c1o2;
            m1 = -mfbac - c2o1* mfbab *  vvz + mfbaa * (c1o1 - vz2);
            m2 = mfbac * c1o2 + mfbab * (vvz + c1o2) + mfbaa * (vz2 + vvz) * c1o2;
            mfbaa = m0;
            mfbab = m1;
            mfbac = m2;
            /////////b//////////////////////////////////////////////////////////////////////////
            m0 = mfbbc * c1o2 + mfbbb * (vvz - c1o2) + mfbba * (vz2 - vvz) * c1o2;
            m1 = -mfbbc - c2o1* mfbbb *  vvz + mfbba * (c1o1 - vz2);
            m2 = mfbbc * c1o2 + mfbbb * (vvz + c1o2) + mfbba * (vz2 + vvz) * c1o2;
            mfbba = m0;
            mfbbb = m1;
            mfbbc = m2;
            /////////b//////////////////////////////////////////////////////////////////////////
            m0 = mfbcc * c1o2 + mfbcb * (vvz - c1o2) + mfbca * (vz2 - vvz) * c1o2;
            m1 = -mfbcc - c2o1* mfbcb *  vvz + mfbca * (c1o1 - vz2);
            m2 = mfbcc * c1o2 + mfbcb * (vvz + c1o2) + mfbca * (vz2 + vvz) * c1o2;
            mfbca = m0;
            mfbcb = m1;
            mfbcc = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfcac * c1o2 + mfcab * (vvz - c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 - vvz) * c1o2;
            m1 = -mfcac - c2o1* mfcab *  vvz + mfcaa                  * (c1o1 - vz2) - c1o3 * oMdrho * vz2;
            m2 = mfcac * c1o2 + mfcab * (vvz + c1o2) + (mfcaa + c1o3 * oMdrho) * (vz2 + vvz) * c1o2;
            mfcaa = m0;
            mfcab = m1;
            mfcac = m2;
            /////////c//////////////////////////////////////////////////////////////////////////
            m0 = mfcbc * c1o2 + mfcbb * (vvz - c1o2) + mfcba * (vz2 - vvz) * c1o2;
            m1 = -mfcbc - c2o1* mfcbb *  vvz + mfcba * (c1o1 - vz2);
            m2 = mfcbc * c1o2 + mfcbb * (vvz + c1o2) + mfcba * (vz2 + vvz) * c1o2;
            mfcba = m0;
            mfcbb = m1;
            mfcbc = m2;
            /////////c//////////////////////////////////////////////////////////////////////////
            m0 = mfccc * c1o2 + mfccb * (vvz - c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 - vvz) * c1o2;
            m1 = -mfccc - c2o1* mfccb *  vvz + mfcca                  * (c1o1 - vz2) - c1o9 * oMdrho * vz2;
            m2 = mfccc * c1o2 + mfccb * (vvz + c1o2) + (mfcca + c1o9 * oMdrho) * (vz2 + vvz) * c1o2;
            mfcca = m0;
            mfccb = m1;
            mfccc = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            m0 = mfaca * c1o2 + mfaba * (vvy - c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
            m1 = -mfaca - c2o1* mfaba *  vvy + mfaaa                  * (c1o1 - vy2) - c1o6 * oMdrho * vy2;
            m2 = mfaca * c1o2 + mfaba * (vvy + c1o2) + (mfaaa + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
            mfaaa = m0;
            mfaba = m1;
            mfaca = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfacb * c1o2 + mfabb * (vvy - c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 - vvy) * c1o2;
            m1 = -mfacb - c2o1* mfabb *  vvy + mfaab                  * (c1o1 - vy2) - c2o3 * oMdrho * vy2;
            m2 = mfacb * c1o2 + mfabb * (vvy + c1o2) + (mfaab + c2o3 * oMdrho) * (vy2 + vvy) * c1o2;
            mfaab = m0;
            mfabb = m1;
            mfacb = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfacc * c1o2 + mfabc * (vvy - c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 - vvy) * c1o2;
            m1 = -mfacc - c2o1* mfabc *  vvy + mfaac                  * (c1o1 - vy2) - c1o6 * oMdrho * vy2;
            m2 = mfacc * c1o2 + mfabc * (vvy + c1o2) + (mfaac + c1o6 * oMdrho) * (vy2 + vvy) * c1o2;
            mfaac = m0;
            mfabc = m1;
            mfacc = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfbca * c1o2 + mfbba * (vvy - c1o2) + mfbaa * (vy2 - vvy) * c1o2;
            m1 = -mfbca - c2o1* mfbba *  vvy + mfbaa * (c1o1 - vy2);
            m2 = mfbca * c1o2 + mfbba * (vvy + c1o2) + mfbaa * (vy2 + vvy) * c1o2;
            mfbaa = m0;
            mfbba = m1;
            mfbca = m2;
            /////////b//////////////////////////////////////////////////////////////////////////
            m0 = mfbcb * c1o2 + mfbbb * (vvy - c1o2) + mfbab * (vy2 - vvy) * c1o2;
            m1 = -mfbcb - c2o1* mfbbb *  vvy + mfbab * (c1o1 - vy2);
            m2 = mfbcb * c1o2 + mfbbb * (vvy + c1o2) + mfbab * (vy2 + vvy) * c1o2;
            mfbab = m0;
            mfbbb = m1;
            mfbcb = m2;
            /////////b//////////////////////////////////////////////////////////////////////////
            m0 = mfbcc * c1o2 + mfbbc * (vvy - c1o2) + mfbac * (vy2 - vvy) * c1o2;
            m1 = -mfbcc - c2o1* mfbbc *  vvy + mfbac * (c1o1 - vy2);
            m2 = mfbcc * c1o2 + mfbbc * (vvy + c1o2) + mfbac * (vy2 + vvy) * c1o2;
            mfbac = m0;
            mfbbc = m1;
            mfbcc = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfcca * c1o2 + mfcba * (vvy - c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
            m1 = -mfcca - c2o1* mfcba *  vvy + mfcaa                   * (c1o1 - vy2) - c1o18 * oMdrho * vy2;
            m2 = mfcca * c1o2 + mfcba * (vvy + c1o2) + (mfcaa + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
            mfcaa = m0;
            mfcba = m1;
            mfcca = m2;
            /////////c//////////////////////////////////////////////////////////////////////////
            m0 = mfccb * c1o2 + mfcbb * (vvy - c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 - vvy) * c1o2;
            m1 = -mfccb - c2o1* mfcbb *  vvy + mfcab                  * (c1o1 - vy2) - c2o9 * oMdrho * vy2;
            m2 = mfccb * c1o2 + mfcbb * (vvy + c1o2) + (mfcab + c2o9 * oMdrho) * (vy2 + vvy) * c1o2;
            mfcab = m0;
            mfcbb = m1;
            mfccb = m2;
            /////////c//////////////////////////////////////////////////////////////////////////
            m0 = mfccc * c1o2 + mfcbc * (vvy - c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 - vvy) * c1o2;
            m1 = -mfccc - c2o1* mfcbc *  vvy + mfcac                   * (c1o1 - vy2) - c1o18 * oMdrho * vy2;
            m2 = mfccc * c1o2 + mfcbc * (vvy + c1o2) + (mfcac + c1o18 * oMdrho) * (vy2 + vvy) * c1o2;
            mfcac = m0;
            mfcbc = m1;
            mfccc = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // X - Dir
            m0 = mfcaa * c1o2 + mfbaa * (vvx - c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfcaa - c2o1* mfbaa *  vvx + mfaaa                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
            m2 = mfcaa * c1o2 + mfbaa * (vvx + c1o2) + (mfaaa + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
            mfaaa = m0;
            mfbaa = m1;
            mfcaa = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfcba * c1o2 + mfbba * (vvx - c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfcba - c2o1* mfbba *  vvx + mfaba                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
            m2 = mfcba * c1o2 + mfbba * (vvx + c1o2) + (mfaba + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
            mfaba = m0;
            mfbba = m1;
            mfcba = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfcca * c1o2 + mfbca * (vvx - c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfcca - c2o1* mfbca *  vvx + mfaca                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
            m2 = mfcca * c1o2 + mfbca * (vvx + c1o2) + (mfaca + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
            mfaca = m0;
            mfbca = m1;
            mfcca = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfcab * c1o2 + mfbab * (vvx - c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfcab - c2o1* mfbab *  vvx + mfaab                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
            m2 = mfcab * c1o2 + mfbab * (vvx + c1o2) + (mfaab + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
            mfaab = m0;
            mfbab = m1;
            mfcab = m2;
            ///////////b////////////////////////////////////////////////////////////////////////
            m0 = mfcbb * c1o2 + mfbbb * (vvx - c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfcbb - c2o1* mfbbb *  vvx + mfabb                  * (c1o1 - vx2) - c4o9 * oMdrho * vx2;
            m2 = mfcbb * c1o2 + mfbbb * (vvx + c1o2) + (mfabb + c4o9 * oMdrho) * (vx2 + vvx) * c1o2;
            mfabb = m0;
            mfbbb = m1;
            mfcbb = m2;
            ///////////b////////////////////////////////////////////////////////////////////////
            m0 = mfccb * c1o2 + mfbcb * (vvx - c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfccb - c2o1* mfbcb *  vvx + mfacb                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
            m2 = mfccb * c1o2 + mfbcb * (vvx + c1o2) + (mfacb + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
            mfacb = m0;
            mfbcb = m1;
            mfccb = m2;
            ////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////
            m0 = mfcac * c1o2 + mfbac * (vvx - c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfcac - c2o1* mfbac *  vvx + mfaac                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
            m2 = mfcac * c1o2 + mfbac * (vvx + c1o2) + (mfaac + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
            mfaac = m0;
            mfbac = m1;
            mfcac = m2;
            ///////////c////////////////////////////////////////////////////////////////////////
            m0 = mfcbc * c1o2 + mfbbc * (vvx - c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfcbc - c2o1* mfbbc *  vvx + mfabc                  * (c1o1 - vx2) - c1o9 * oMdrho * vx2;
            m2 = mfcbc * c1o2 + mfbbc * (vvx + c1o2) + (mfabc + c1o9 * oMdrho) * (vx2 + vvx) * c1o2;
            mfabc = m0;
            mfbbc = m1;
            mfcbc = m2;
            ///////////c////////////////////////////////////////////////////////////////////////
            m0 = mfccc * c1o2 + mfbcc * (vvx - c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 - vvx) * c1o2;
            m1 = -mfccc - c2o1* mfbcc *  vvx + mfacc                   * (c1o1 - vx2) - c1o36 * oMdrho * vx2;
            m2 = mfccc * c1o2 + mfbcc * (vvx + c1o2) + (mfacc + c1o36 * oMdrho) * (vx2 + vvx) * c1o2;
            mfacc = m0;
            mfbcc = m1;
            mfccc = m2;
            ////////////////////////////////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////////////////////////////////
            //// test rundungsfehler....
            //mfcbb = (mfcbb * rho + c2over27* (rho-one));
            //mfabb = (mfabb * rho + c2over27* (rho-one));
            //mfbcb = (mfbcb * rho + c2over27* (rho-one));
            //mfbab = (mfbab * rho + c2over27* (rho-one));
            //mfbbc = (mfbbc * rho + c2over27* (rho-one));
            //mfbba = (mfbba * rho + c2over27* (rho-one));
            //mfccb = (mfccb * rho + c1over54* (rho-one));
            //mfaab = (mfaab * rho + c1over54* (rho-one));
            //mfcab = (mfcab * rho + c1over54* (rho-one));
            //mfacb = (mfacb * rho + c1over54* (rho-one));
            //mfcbc = (mfcbc * rho + c1over54* (rho-one));
            //mfaba = (mfaba * rho + c1over54* (rho-one));
            //mfcba = (mfcba * rho + c1over54* (rho-one));
            //mfabc = (mfabc * rho + c1over54* (rho-one));
            //mfbcc = (mfbcc * rho + c1over54* (rho-one));
            //mfbaa = (mfbaa * rho + c1over54* (rho-one));
            //mfbca = (mfbca * rho + c1over54* (rho-one));
            //mfbac = (mfbac * rho + c1over54* (rho-one));
            //mfbbb = (mfbbb * rho + c8over27* (rho-one));
            //mfccc = (mfccc * rho + c1over216*(rho-one));
            //mfaac = (mfaac * rho + c1over216*(rho-one));
            //mfcac = (mfcac * rho + c1over216*(rho-one));
            //mfacc = (mfacc * rho + c1over216*(rho-one));
            //mfcca = (mfcca * rho + c1over216*(rho-one));
            //mfaaa = (mfaaa * rho + c1over216*(rho-one));
            //mfcaa = (mfcaa * rho + c1over216*(rho-one));
            //mfaca = (mfaca * rho + c1over216*(rho-one));
            real drhoPost =
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                    ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
            mfbbb += drho - drhoPost;
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
//! \}
