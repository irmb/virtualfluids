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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Cumulant27chim.cu
//! \ingroup GPU
//! \author Martin Schoenherr
//=======================================================================================
/* Device code */
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

#include "math.h"

#include "lbm/Chimera.h"


////////////////////////////////////////////////////////////////////////////////
__global__ void Cumulant_One_preconditioned_errorDiffusion_chim_Comp_SP_27(
    real omega,
    unsigned int* bcMatD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* DDStart,
    unsigned long long numberOfLBnodes,
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

    if (k<numberOfLBnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if (BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/)
        {
            Distributions27 D;
            if (EvenOrOdd == true)
            {
                D.f[DIR_P00] = &DDStart[DIR_P00 * numberOfLBnodes];
                D.f[DIR_M00] = &DDStart[DIR_M00 * numberOfLBnodes];
                D.f[DIR_0P0] = &DDStart[DIR_0P0 * numberOfLBnodes];
                D.f[DIR_0M0] = &DDStart[DIR_0M0 * numberOfLBnodes];
                D.f[DIR_00P] = &DDStart[DIR_00P * numberOfLBnodes];
                D.f[DIR_00M] = &DDStart[DIR_00M * numberOfLBnodes];
                D.f[DIR_PP0] = &DDStart[DIR_PP0 * numberOfLBnodes];
                D.f[DIR_MM0] = &DDStart[DIR_MM0 * numberOfLBnodes];
                D.f[DIR_PM0] = &DDStart[DIR_PM0 * numberOfLBnodes];
                D.f[DIR_MP0] = &DDStart[DIR_MP0 * numberOfLBnodes];
                D.f[DIR_P0P] = &DDStart[DIR_P0P * numberOfLBnodes];
                D.f[DIR_M0M] = &DDStart[DIR_M0M * numberOfLBnodes];
                D.f[DIR_P0M] = &DDStart[DIR_P0M * numberOfLBnodes];
                D.f[DIR_M0P] = &DDStart[DIR_M0P * numberOfLBnodes];
                D.f[DIR_0PP] = &DDStart[DIR_0PP * numberOfLBnodes];
                D.f[DIR_0MM] = &DDStart[DIR_0MM * numberOfLBnodes];
                D.f[DIR_0PM] = &DDStart[DIR_0PM * numberOfLBnodes];
                D.f[DIR_0MP] = &DDStart[DIR_0MP * numberOfLBnodes];
                D.f[DIR_000] = &DDStart[DIR_000 * numberOfLBnodes];
                D.f[DIR_PPP] = &DDStart[DIR_PPP * numberOfLBnodes];
                D.f[DIR_MMP] = &DDStart[DIR_MMP * numberOfLBnodes];
                D.f[DIR_PMP] = &DDStart[DIR_PMP * numberOfLBnodes];
                D.f[DIR_MPP] = &DDStart[DIR_MPP * numberOfLBnodes];
                D.f[DIR_PPM] = &DDStart[DIR_PPM * numberOfLBnodes];
                D.f[DIR_MMM] = &DDStart[DIR_MMM * numberOfLBnodes];
                D.f[DIR_PMM] = &DDStart[DIR_PMM * numberOfLBnodes];
                D.f[DIR_MPM] = &DDStart[DIR_MPM * numberOfLBnodes];
            }
            else
            {
                D.f[DIR_M00] = &DDStart[DIR_P00 * numberOfLBnodes];
                D.f[DIR_P00] = &DDStart[DIR_M00 * numberOfLBnodes];
                D.f[DIR_0M0] = &DDStart[DIR_0P0 * numberOfLBnodes];
                D.f[DIR_0P0] = &DDStart[DIR_0M0 * numberOfLBnodes];
                D.f[DIR_00M] = &DDStart[DIR_00P * numberOfLBnodes];
                D.f[DIR_00P] = &DDStart[DIR_00M * numberOfLBnodes];
                D.f[DIR_MM0] = &DDStart[DIR_PP0 * numberOfLBnodes];
                D.f[DIR_PP0] = &DDStart[DIR_MM0 * numberOfLBnodes];
                D.f[DIR_MP0] = &DDStart[DIR_PM0 * numberOfLBnodes];
                D.f[DIR_PM0] = &DDStart[DIR_MP0 * numberOfLBnodes];
                D.f[DIR_M0M] = &DDStart[DIR_P0P * numberOfLBnodes];
                D.f[DIR_P0P] = &DDStart[DIR_M0M * numberOfLBnodes];
                D.f[DIR_M0P] = &DDStart[DIR_P0M * numberOfLBnodes];
                D.f[DIR_P0M] = &DDStart[DIR_M0P * numberOfLBnodes];
                D.f[DIR_0MM] = &DDStart[DIR_0PP * numberOfLBnodes];
                D.f[DIR_0PP] = &DDStart[DIR_0MM * numberOfLBnodes];
                D.f[DIR_0MP] = &DDStart[DIR_0PM * numberOfLBnodes];
                D.f[DIR_0PM] = &DDStart[DIR_0MP * numberOfLBnodes];
                D.f[DIR_000] = &DDStart[DIR_000 * numberOfLBnodes];
                D.f[DIR_MMM] = &DDStart[DIR_PPP * numberOfLBnodes];
                D.f[DIR_PPM] = &DDStart[DIR_MMP * numberOfLBnodes];
                D.f[DIR_MPM] = &DDStart[DIR_PMP * numberOfLBnodes];
                D.f[DIR_PMM] = &DDStart[DIR_MPP * numberOfLBnodes];
                D.f[DIR_MMP] = &DDStart[DIR_PPM * numberOfLBnodes];
                D.f[DIR_PPP] = &DDStart[DIR_MMM * numberOfLBnodes];
                D.f[DIR_MPP] = &DDStart[DIR_PMM * numberOfLBnodes];
                D.f[DIR_PMP] = &DDStart[DIR_MPM * numberOfLBnodes];
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
            real mfcbb = (D.f[DIR_P00])[k];//[ke   ];// +  c2over27 ;(D.f[DIR_P00])[k  ];//ke
            real mfabb = (D.f[DIR_M00])[kw];//[kw   ];// +  c2over27 ;(D.f[DIR_M00])[kw ];
            real mfbcb = (D.f[DIR_0P0])[k];//[kn   ];// +  c2over27 ;(D.f[DIR_0P0])[k  ];//kn
            real mfbab = (D.f[DIR_0M0])[ks];//[ks   ];// +  c2over27 ;(D.f[DIR_0M0])[ks ];
            real mfbbc = (D.f[DIR_00P])[k];//[kt   ];// +  c2over27 ;(D.f[DIR_00P])[k  ];//kt
            real mfbba = (D.f[DIR_00M])[kb];//[kb   ];// +  c2over27 ;(D.f[DIR_00M])[kb ];
            real mfccb = (D.f[DIR_PP0])[k];//[kne  ];// +  c1over54 ;(D.f[DIR_PP0])[k  ];//kne
            real mfaab = (D.f[DIR_MM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[DIR_MM0])[ksw];
            real mfcab = (D.f[DIR_PM0])[ks];//[kse  ];// +  c1over54 ;(D.f[DIR_PM0])[ks ];//kse
            real mfacb = (D.f[DIR_MP0])[kw];//[knw  ];// +  c1over54 ;(D.f[DIR_MP0])[kw ];//knw
            real mfcbc = (D.f[DIR_P0P])[k];//[kte  ];// +  c1over54 ;(D.f[DIR_P0P])[k  ];//kte
            real mfaba = (D.f[DIR_M0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[DIR_M0M])[kbw];
            real mfcba = (D.f[DIR_P0M])[kb];//[kbe  ];// +  c1over54 ;(D.f[DIR_P0M])[kb ];//kbe
            real mfabc = (D.f[DIR_M0P])[kw];//[ktw  ];// +  c1over54 ;(D.f[DIR_M0P])[kw ];//ktw
            real mfbcc = (D.f[DIR_0PP])[k];//[ktn  ];// +  c1over54 ;(D.f[DIR_0PP])[k  ];//ktn
            real mfbaa = (D.f[DIR_0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[DIR_0MM])[kbs];
            real mfbca = (D.f[DIR_0PM])[kb];//[kbn  ];// +  c1over54 ;(D.f[DIR_0PM])[kb ];//kbn
            real mfbac = (D.f[DIR_0MP])[ks];//[kts  ];// +  c1over54 ;(D.f[DIR_0MP])[ks ];//kts
            real mfbbb = (D.f[DIR_000])[k];//[kzero];// +  c8over27 ;(D.f[DIR_000])[k  ];//kzero
            real mfccc = (D.f[DIR_PPP])[k];//[ktne ];// +  c1over216;(D.f[DIR_PPP])[k  ];//ktne
            real mfaac = (D.f[DIR_MMP])[ksw];//[ktsw ];// +  c1over216;(D.f[DIR_MMP])[ksw];//ktsw
            real mfcac = (D.f[DIR_PMP])[ks];//[ktse ];// +  c1over216;(D.f[DIR_PMP])[ks ];//ktse
            real mfacc = (D.f[DIR_MPP])[kw];//[ktnw ];// +  c1over216;(D.f[DIR_MPP])[kw ];//ktnw
            real mfcca = (D.f[DIR_PPM])[kb];//[kbne ];// +  c1over216;(D.f[DIR_PPM])[kb ];//kbne
            real mfaaa = (D.f[DIR_MMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[DIR_MMM])[kbsw];
            real mfcaa = (D.f[DIR_PMM])[kbs];//[kbse ];// +  c1over216;(D.f[DIR_PMM])[kbs];//kbse
            real mfaca = (D.f[DIR_MPM])[kbw];//[kbnw ];// +  c1over216;(D.f[DIR_MPM])[kbw];//kbnw
                                               ////////////////////////////////////////////////////////////////////////////////////
            real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

            real rho = c1o1 + drho;
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
            //real oMdrho = c1o1; // comp special
            //real m0, m1, m2;
            real vx2;
            real vy2;
            real vz2;
            vx2 = vvx*vvx;
            vy2 = vvy*vvy;
            vz2 = vvz*vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            //real wadjust;
            //real qudricLimitP = c1o100;// * 0.0001f;
            //real qudricLimitM = c1o100;// * 0.0001f;
            //real qudricLimitD = c1o100;// * 0.001f;
            //real s9 = minusomega;
            //test
            //s9 = 0.;


            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real EQcbb = c0o1;
            real EQabb = c0o1;
            real EQbcb = c0o1;
            real EQbab = c0o1;
            real EQbbc = c0o1;
            real EQbba = c0o1;
            real EQccb = c0o1;
            real EQaab = c0o1;
            real EQcab = c0o1;
            real EQacb = c0o1;
            real EQcbc = c0o1;
            real EQaba = c0o1;
            real EQcba = c0o1;
            real EQabc = c0o1;
            real EQbcc = c0o1;
            real EQbaa = c0o1;
            real EQbca = c0o1;
            real EQbac = c0o1;
            real EQbbb = c0o1;
            real EQccc = drho * c1o27;
            real EQaac = drho * c1o3;
            real EQcac = drho * c1o9;
            real EQacc = drho * c1o9;
            real EQcca = drho * c1o9;
            real EQaaa = drho;
            real EQcaa = drho * c1o3;
            real EQaca = drho * c1o3;
            ////////////////////////////////////////////////////////////////////////////////////
            vf::lbm::backwardChimeraWithK(EQaaa, EQaab, EQaac, vvz, vz2, c1o1);
            vf::lbm::backwardChimeraWithK(EQaca, EQacb, EQacc, vvz, vz2, c1o3);
            ///////////////////////////////////////////////////////////
            EQcaa = EQaca; EQcab = EQacb; EQcac = EQacc;
            ///////////////////////////////////////////////////////////
            vf::lbm::backwardChimeraWithK(EQcca, EQccb, EQccc, vvz, vz2, c1o9);

            vf::lbm::backwardChimeraWithK(EQaaa, EQaba, EQaca, vvy, vy2, c1o6);
            vf::lbm::backwardChimeraWithK(EQaab, EQabb, EQacb, vvy, vy2, c2o3);
            vf::lbm::backwardChimeraWithK(EQaac, EQabc, EQacc, vvy, vy2, c1o6);
            vf::lbm::backwardChimeraWithK(EQcaa, EQcba, EQcca, vvy, vy2, c1o18);
            vf::lbm::backwardChimeraWithK(EQcab, EQcbb, EQccb, vvy, vy2, c2o9);
            vf::lbm::backwardChimeraWithK(EQcac, EQcbc, EQccc, vvy, vy2, c1o18);

            vf::lbm::backwardChimeraWithK(EQaaa, EQbaa, EQcaa, vvx, vx2, c1o36);
            vf::lbm::backwardChimeraWithK(EQaab, EQbab, EQcab, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQaac, EQbac, EQcac, vvx, vx2, c1o36);
            vf::lbm::backwardChimeraWithK(EQaba, EQbba, EQcba, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQabb, EQbbb, EQcbb, vvx, vx2, c4o9);
            vf::lbm::backwardChimeraWithK(EQabc, EQbbc, EQcbc, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQaca, EQbca, EQcca, vvx, vx2, c1o36);
            vf::lbm::backwardChimeraWithK(EQacb, EQbcb, EQccb, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQacc, EQbcc, EQccc, vvx, vx2, c1o36);

            ////////////////////////////////////////////////////////////////////////////////////
            //Pre-condition
            mfcbb -= EQcbb;
            mfabb -= EQabb;
            mfbcb -= EQbcb;
            mfbab -= EQbab;
            mfbbc -= EQbbc;
            mfbba -= EQbba;
            mfccb -= EQccb;
            mfaab -= EQaab;
            mfcab -= EQcab;
            mfacb -= EQacb;
            mfcbc -= EQcbc;
            mfaba -= EQaba;
            mfcba -= EQcba;
            mfabc -= EQabc;
            mfbcc -= EQbcc;
            mfbaa -= EQbaa;
            mfbca -= EQbca;
            mfbac -= EQbac;
            mfbbb -= EQbbb;
            mfccc -= EQccc;
            mfaac -= EQaac;
            mfcac -= EQcac;
            mfacc -= EQacc;
            mfcca -= EQcca;
            mfaaa -= EQaaa;
            mfcaa -= EQcaa;
            mfaca -= EQaca;

            ////////////////////////////////////////////////////////////////////////////////////
            //Hin
            ////////////////////////////////////////////////////////////////////////////////////
            vf::lbm::forwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
            vf::lbm::forwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
            vf::lbm::forwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
            vf::lbm::forwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
            vf::lbm::forwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
            vf::lbm::forwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
            vf::lbm::forwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
            vf::lbm::forwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
            vf::lbm::forwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

            vf::lbm::forwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
            vf::lbm::forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
            vf::lbm::forwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
            vf::lbm::forwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
            vf::lbm::forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
            vf::lbm::forwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
            vf::lbm::forwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
            vf::lbm::forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
            vf::lbm::forwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

            vf::lbm::forwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
            vf::lbm::forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
            vf::lbm::forwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
            vf::lbm::forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
            vf::lbm::forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
            vf::lbm::forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
            vf::lbm::forwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
            vf::lbm::forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
            vf::lbm::forwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

            //////////////////////////////////////////////////////////////////////////////////////
            ////Hin
            //////////////////////////////////////////////////////////////////////////////////////
            //// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Z - Dir
            //forwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c4o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Y - Dir
            //forwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o18);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c2o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c2o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// X - Dir
            //forwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, one);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            real OxxPyyPzz = c1o1; //omega; // one;	//set the bulk viscosity one is high / two is very low and zero is (too) high

            ////////////////////////////////////////////////////////////
            //3.
            //////////////////////////////
            real OxyyPxzz = c1o1;
            real OxyyMxzz = c1o1;
            //real Oxyz = c1o1;
            ////////////////////////////////////////////////////////////
            //4.
            //////////////////////////////
            real O4 = c1o1;
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
            real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
            real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
            real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

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




            //2.
            // linear combinations
            real mxxPyyPzz = mfcaa + mfaca + mfaac;
            real mxxMyy = mfcaa - mfaca;
            real mxxMzz = mfcaa - mfaac;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
            {
                real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
                real dyuy = dxux + omega * c3o2 * mxxMyy;
                real dzuz = dxux + omega * c3o2 * mxxMzz;

                //relax
                mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
                mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
                mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////no correction
            //mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);//-magicBulk*OxxPyyPzz;
            //mxxMyy += -(-omega) * (-mxxMyy);
            //mxxMzz += -(-omega) * (-mxxMzz);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            mfabb += omega * (-mfabb);
            mfbab += omega * (-mfbab);
            mfbba += omega * (-mfbba);

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
            mfbbb += OxyyMxzz * (-mfbbb);
            mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
            mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
            mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
            mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
            mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
            mxyyMxzz += OxyyMxzz * (-mxyyMxzz);
            //////////////////////////////////////////////////////////////////////////

            mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
            mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
            mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
            mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
            mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
            mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

            //4.
            //////////////////////////////////////////////////////////////////////////
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
            mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
            mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
            mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

            mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
            mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
            mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

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

            ////////////////////////////////////////////////////////////////////////////////////
            //the force be with you
            mfbaa = -mfbaa;
            mfaba = -mfaba;
            mfaab = -mfaab;
            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            //back
            ////////////////////////////////////////////////////////////////////////////////////
            vf::lbm::backwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
            vf::lbm::backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
            vf::lbm::backwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
            vf::lbm::backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
            vf::lbm::backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
            vf::lbm::backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
            vf::lbm::backwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
            vf::lbm::backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
            vf::lbm::backwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

            vf::lbm::backwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
            vf::lbm::backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
            vf::lbm::backwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
            vf::lbm::backwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
            vf::lbm::backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
            vf::lbm::backwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
            vf::lbm::backwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
            vf::lbm::backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
            vf::lbm::backwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

            vf::lbm::backwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
            vf::lbm::backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
            vf::lbm::backwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
            vf::lbm::backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
            vf::lbm::backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
            vf::lbm::backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
            vf::lbm::backwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
            vf::lbm::backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
            vf::lbm::backwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

            ////////////////////////////////////////////////////////////////////////////////////
            //mfcbb += EQcbb;
            //mfabb += EQabb;
            //mfbcb += EQbcb;
            //mfbab += EQbab;
            //mfbbc += EQbbc;
            //mfbba += EQbba;
            //mfccb += EQccb;
            //mfaab += EQaab;
            //mfcab += EQcab;
            //mfacb += EQacb;
            //mfcbc += EQcbc;
            //mfaba += EQaba;
            //mfcba += EQcba;
            //mfabc += EQabc;
            //mfbcc += EQbcc;
            //mfbaa += EQbaa;
            //mfbca += EQbca;
            //mfbac += EQbac;
            //mfbbb += EQbbb;
            //mfccc += EQccc;
            //mfaac += EQaac;
            //mfcac += EQcac;
            //mfacc += EQacc;
            //mfcca += EQcca;
            //mfaaa += EQaaa;
            //mfcaa += EQcaa;
            //mfaca += EQaca;
            ////////////////////////////////////////////////////////////////////////////////////
            ////Error diffusion
            real fTEMP = mfbbb + EQbbb;
            real delta0 = mfbbb - (fTEMP - EQbbb);
            delta0 *= c1o4;
            mfbbb = fTEMP;


            fTEMP = mfcbb + EQcbb;
            real deltacbb = mfcbb - (fTEMP - EQcbb);
            mfcbb = fTEMP;
            //mfcbb+=EQcbb;

            fTEMP = mfabb + EQabb;
            real deltaabb = mfabb - (fTEMP - EQabb);
            mfabb = fTEMP;
            //mfabb+=EQabb;

            fTEMP = mfbcb + EQbcb;
            real deltabcb = mfbcb - (fTEMP - EQbcb);
            mfbcb = fTEMP;
            //mfbcb+=EQbcb;

            fTEMP = mfbab + EQbab;
            real deltabab = mfbab - (fTEMP - EQbab);
            mfbab = fTEMP;
            //mfbab+=EQbab;

            fTEMP = mfbbc + EQbbc;
            real deltabbc = mfbbc - (fTEMP - EQbbc);
            mfbbc = fTEMP;
            //mfbbc+=EQbbc;

            fTEMP = mfbba + EQbba;
            real deltabba = mfbba - (fTEMP - EQbba);
            mfbba = fTEMP;
            //mfbba+=EQbba;

            EQccb += (delta0 + c1o2*(deltacbb + deltabcb));
            fTEMP = mfccb + EQccb;
            real deltaccb = mfccb - (fTEMP - EQccb);
            mfccb = fTEMP;
            //mfccb+=EQccb+(delta0+c1o2*(deltacbb+deltabcb));

            EQaab += (delta0 + c1o2*(deltaabb + deltabab));
            fTEMP = mfaab + EQaab;
            real deltaaab = mfaab - (fTEMP - EQaab);
            mfaab = fTEMP;
            //mfaab+=EQaab+(delta0+c1o2*(deltaabb+deltabab));

            EQcab += (delta0 + c1o2*(deltacbb + deltabab));
            fTEMP = mfcab + EQcab;
            real deltacab = mfcab - (fTEMP - EQcab);
            mfcab = fTEMP;
            //mfcab+=EQcab+(delta0+c1o2*(deltacbb+deltabab));

            EQacb += (delta0 + c1o2*(deltaabb + deltabcb));
            fTEMP = mfacb + EQacb;
            real deltaacb = mfacb - (fTEMP - EQacb);
            mfacb = fTEMP;
            //mfacb+=EQacb+(delta0+c1o2*(deltaabb+deltabcb));

            EQcbc += (delta0 + c1o2*(deltacbb + deltabbc));
            fTEMP = mfcbc + EQcbc;
            real deltacbc = mfcbc - (fTEMP - EQcbc);
            mfcbc = fTEMP;
            //mfcbc+=EQcbc+(delta0+c1o2*(deltacbb+deltabbc));

            EQaba += (delta0 + c1o2*(deltaabb + deltabba));
            fTEMP = mfaba + EQaba;
            real deltaaba = mfaba - (fTEMP - EQaba);
            mfaba = fTEMP;
            //mfaba+=EQaba+(delta0+c1o2*(deltaabb+deltabba));

            EQcba += (delta0 + c1o2*(deltacbb + deltabba));
            fTEMP = mfcba + EQcba;
            real deltacba = mfcba - (fTEMP - EQcba);
            mfcba = fTEMP;
            //mfcba+=EQcba+(delta0+c1o2*(deltacbb+deltabba));

            EQabc += (delta0 + c1o2*(deltaabb + deltabbc));
            fTEMP = mfabc + EQabc;
            real deltaabc = mfabc - (fTEMP - EQabc);
            mfabc = fTEMP;
            //mfabc+=EQabc+(delta0+c1o2*(deltaabb+deltabbc));

            EQbcc += (delta0 + c1o2*(deltabcb + deltabbc));
            fTEMP = mfbcc + EQbcc;
            real deltabcc = mfbcc - (fTEMP - EQbcc);
            mfbcc = fTEMP;
            //mfbcc+=EQbcc+(delta0+c1o2*(deltabcb+deltabbc));

            EQbaa += (delta0 + c1o2*(deltabab + deltabba));
            fTEMP = mfbaa + EQbaa;
            real deltabaa = mfbaa - (fTEMP - EQbaa);
            mfbaa = fTEMP;
            //mfbaa+=EQbaa+(delta0+c1o2*(deltabab+deltabba));

            EQbca += (delta0 + c1o2*(deltabcb + deltabba));
            fTEMP = mfbca + EQbca;
            real deltabca = mfbca - (fTEMP - EQbca);
            mfbca = fTEMP;
            //mfbca+=EQbca+(delta0+c1o2*(deltabcb+deltabba));

            EQbac += (delta0 + c1o2*(deltabab + deltabbc));
            fTEMP = mfbac + EQbac;
            real deltabac = mfbac - (fTEMP - EQbac);
            mfbac = fTEMP;
            //mfbac+=EQbac+(delta0+c1o2*(deltabab+deltabbc));

            mfccc += EQccc - (delta0 + c1o4*(deltacbb + deltabcb + deltabbc) - c1o2*(deltabcc + deltacbc + deltaccb));
            mfaac += EQaac - (delta0 + c1o4*(deltaabb + deltabab + deltabbc) - c1o2*(deltabac + deltaabc + deltaaab));
            mfcac += EQcac - (delta0 + c1o4*(deltacbb + deltabab + deltabbc) - c1o2*(deltabac + deltacbc + deltacab));
            mfacc += EQacc - (delta0 + c1o4*(deltaabb + deltabcb + deltabbc) - c1o2*(deltabcc + deltaabc + deltaacb));
            mfcca += EQcca - (delta0 + c1o4*(deltacbb + deltabcb + deltabba) - c1o2*(deltabca + deltacba + deltaccb));
            mfaaa += EQaaa - (delta0 + c1o4*(deltaabb + deltabab + deltabba) - c1o2*(deltabaa + deltaaba + deltaaab));
            mfcaa += EQcaa - (delta0 + c1o4*(deltacbb + deltabab + deltabba) - c1o2*(deltabaa + deltacba + deltacab));
            mfaca += EQaca - (delta0 + c1o4*(deltaabb + deltabcb + deltabba) - c1o2*(deltabca + deltaaba + deltaacb));



            //////////////////////////////////////////////////////////////////////////////////////
            ////back
            //////////////////////////////////////////////////////////////////////////////////////
            ////mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Z - Dir
            //backwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, one);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o3);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            ////mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Y - Dir
            //backwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaab, mfabb, mfacb, vvy, vy2, c2o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbaa, mfbba, mfbca, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbab, mfbbb, mfbcb, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbac, mfbbc, mfbcc, vvz, vz2);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o18);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcab, mfcbb, mfccb, vvy, vy2, c2o9);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            ////mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// X - Dir
            //backwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaba, mfbba, mfcba, vvx, vx2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaab, mfbab, mfcab, vvx, vx2, c1o9);
            /////////////b////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfabb, mfbbb, mfcbb, vvx, vx2, c4o9);
            /////////////b////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfacb, mfbcb, mfccb, vvx, vx2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o36);
            /////////////c////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfabc, mfbbc, mfcbc, vvx, vx2, c1o9);
            /////////////c////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////////////////////////
            //real drhoPost =
            //	((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
            //	(((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
            //		((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
            //mfbbb += drho - drhoPost;
            ////////////////////////////////////////////////////////////////////////////////////
            (D.f[DIR_P00])[k] = mfabb;//(D.f[ DIR_P00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ DIR_P00   ])[k   ]                                                                     
            (D.f[DIR_M00])[kw] = mfcbb;//(D.f[ DIR_M00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ DIR_M00   ])[kw  ]                                                                   
            (D.f[DIR_0P0])[k] = mfbab;//(D.f[ DIR_0P0   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ DIR_0P0   ])[k   ]
            (D.f[DIR_0M0])[ks] = mfbcb;//(D.f[ DIR_0M0   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ DIR_0M0   ])[ks  ]
            (D.f[DIR_00P])[k] = mfbba;//(D.f[ DIR_00P   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ DIR_00P   ])[k   ]
            (D.f[DIR_00M])[kb] = mfbbc;//(D.f[ DIR_00M   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ DIR_00M   ])[kb  ]
            (D.f[DIR_PP0])[k] = mfaab;//(D.f[ DIR_PP0  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ DIR_PP0  ])[k   ]
            (D.f[DIR_MM0])[ksw] = mfccb;//(D.f[ DIR_MM0  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ DIR_MM0  ])[ksw ]
            (D.f[DIR_PM0])[ks] = mfacb;//(D.f[ DIR_PM0  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ DIR_PM0  ])[ks  ]
            (D.f[DIR_MP0])[kw] = mfcab;//(D.f[ DIR_MP0  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ DIR_MP0  ])[kw  ]
            (D.f[DIR_P0P])[k] = mfaba;//(D.f[ DIR_P0P  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ DIR_P0P  ])[k   ]
            (D.f[DIR_M0M])[kbw] = mfcbc;//(D.f[ DIR_M0M  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ DIR_M0M  ])[kbw ]
            (D.f[DIR_P0M])[kb] = mfabc;//(D.f[ DIR_P0M  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ DIR_P0M  ])[kb  ]
            (D.f[DIR_M0P])[kw] = mfcba;//(D.f[ DIR_M0P  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ DIR_M0P  ])[kw  ]
            (D.f[DIR_0PP])[k] = mfbaa;//(D.f[ DIR_0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ DIR_0PP  ])[k   ]
            (D.f[DIR_0MM])[kbs] = mfbcc;//(D.f[ DIR_0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ DIR_0MM  ])[kbs ]
            (D.f[DIR_0PM])[kb] = mfbac;//(D.f[ DIR_0PM  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ DIR_0PM  ])[kb  ]
            (D.f[DIR_0MP])[ks] = mfbca;//(D.f[ DIR_0MP  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ DIR_0MP  ])[ks  ]
            (D.f[DIR_000])[k] = mfbbb;//(D.f[ DIR_000])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ DIR_000])[k   ]
            (D.f[DIR_PPP])[k] = mfaaa;//(D.f[ DIR_PPP ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ DIR_PPP ])[k   ]
            (D.f[DIR_PMP])[ks] = mfaca;//(D.f[ DIR_PMP ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ DIR_PMP ])[ks  ]
            (D.f[DIR_PPM])[kb] = mfaac;//(D.f[ DIR_PPM ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ DIR_PPM ])[kb  ]
            (D.f[DIR_PMM])[kbs] = mfacc;//(D.f[ DIR_PMM ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ DIR_PMM ])[kbs ]
            (D.f[DIR_MPP])[kw] = mfcaa;//(D.f[ DIR_MPP ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ DIR_MPP ])[kw  ]
            (D.f[DIR_MMP])[ksw] = mfcca;//(D.f[ DIR_MMP ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ DIR_MMP ])[ksw ]
            (D.f[DIR_MPM])[kbw] = mfcac;//(D.f[ DIR_MPM ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ DIR_MPM ])[kbw ]
            (D.f[DIR_MMM])[kbsw] = mfccc;//(D.f[ DIR_MMM ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ DIR_MMM ])[kbsw]
                                        ////////////////////////////////////////////////////////////////////////////////////
        }
    }
}
////////////////////////////////////////////////////////////////////////////////








































////////////////////////////////////////////////////////////////////////////////
__global__ void Cumulant_One_preconditioned_chim_Comp_SP_27(
    real omega,
    unsigned int* bcMatD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* DDStart,
    unsigned long long numberOfLBnodes,
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

    if (k<numberOfLBnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if (BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/)
        {
            Distributions27 D;
            if (EvenOrOdd == true)
            {
                D.f[DIR_P00] = &DDStart[DIR_P00 * numberOfLBnodes];
                D.f[DIR_M00] = &DDStart[DIR_M00 * numberOfLBnodes];
                D.f[DIR_0P0] = &DDStart[DIR_0P0 * numberOfLBnodes];
                D.f[DIR_0M0] = &DDStart[DIR_0M0 * numberOfLBnodes];
                D.f[DIR_00P] = &DDStart[DIR_00P * numberOfLBnodes];
                D.f[DIR_00M] = &DDStart[DIR_00M * numberOfLBnodes];
                D.f[DIR_PP0] = &DDStart[DIR_PP0 * numberOfLBnodes];
                D.f[DIR_MM0] = &DDStart[DIR_MM0 * numberOfLBnodes];
                D.f[DIR_PM0] = &DDStart[DIR_PM0 * numberOfLBnodes];
                D.f[DIR_MP0] = &DDStart[DIR_MP0 * numberOfLBnodes];
                D.f[DIR_P0P] = &DDStart[DIR_P0P * numberOfLBnodes];
                D.f[DIR_M0M] = &DDStart[DIR_M0M * numberOfLBnodes];
                D.f[DIR_P0M] = &DDStart[DIR_P0M * numberOfLBnodes];
                D.f[DIR_M0P] = &DDStart[DIR_M0P * numberOfLBnodes];
                D.f[DIR_0PP] = &DDStart[DIR_0PP * numberOfLBnodes];
                D.f[DIR_0MM] = &DDStart[DIR_0MM * numberOfLBnodes];
                D.f[DIR_0PM] = &DDStart[DIR_0PM * numberOfLBnodes];
                D.f[DIR_0MP] = &DDStart[DIR_0MP * numberOfLBnodes];
                D.f[DIR_000] = &DDStart[DIR_000 * numberOfLBnodes];
                D.f[DIR_PPP] = &DDStart[DIR_PPP * numberOfLBnodes];
                D.f[DIR_MMP] = &DDStart[DIR_MMP * numberOfLBnodes];
                D.f[DIR_PMP] = &DDStart[DIR_PMP * numberOfLBnodes];
                D.f[DIR_MPP] = &DDStart[DIR_MPP * numberOfLBnodes];
                D.f[DIR_PPM] = &DDStart[DIR_PPM * numberOfLBnodes];
                D.f[DIR_MMM] = &DDStart[DIR_MMM * numberOfLBnodes];
                D.f[DIR_PMM] = &DDStart[DIR_PMM * numberOfLBnodes];
                D.f[DIR_MPM] = &DDStart[DIR_MPM * numberOfLBnodes];
            }
            else
            {
                D.f[DIR_M00] = &DDStart[DIR_P00 * numberOfLBnodes];
                D.f[DIR_P00] = &DDStart[DIR_M00 * numberOfLBnodes];
                D.f[DIR_0M0] = &DDStart[DIR_0P0 * numberOfLBnodes];
                D.f[DIR_0P0] = &DDStart[DIR_0M0 * numberOfLBnodes];
                D.f[DIR_00M] = &DDStart[DIR_00P * numberOfLBnodes];
                D.f[DIR_00P] = &DDStart[DIR_00M * numberOfLBnodes];
                D.f[DIR_MM0] = &DDStart[DIR_PP0 * numberOfLBnodes];
                D.f[DIR_PP0] = &DDStart[DIR_MM0 * numberOfLBnodes];
                D.f[DIR_MP0] = &DDStart[DIR_PM0 * numberOfLBnodes];
                D.f[DIR_PM0] = &DDStart[DIR_MP0 * numberOfLBnodes];
                D.f[DIR_M0M] = &DDStart[DIR_P0P * numberOfLBnodes];
                D.f[DIR_P0P] = &DDStart[DIR_M0M * numberOfLBnodes];
                D.f[DIR_M0P] = &DDStart[DIR_P0M * numberOfLBnodes];
                D.f[DIR_P0M] = &DDStart[DIR_M0P * numberOfLBnodes];
                D.f[DIR_0MM] = &DDStart[DIR_0PP * numberOfLBnodes];
                D.f[DIR_0PP] = &DDStart[DIR_0MM * numberOfLBnodes];
                D.f[DIR_0MP] = &DDStart[DIR_0PM * numberOfLBnodes];
                D.f[DIR_0PM] = &DDStart[DIR_0MP * numberOfLBnodes];
                D.f[DIR_000] = &DDStart[DIR_000 * numberOfLBnodes];
                D.f[DIR_MMM] = &DDStart[DIR_PPP * numberOfLBnodes];
                D.f[DIR_PPM] = &DDStart[DIR_MMP * numberOfLBnodes];
                D.f[DIR_MPM] = &DDStart[DIR_PMP * numberOfLBnodes];
                D.f[DIR_PMM] = &DDStart[DIR_MPP * numberOfLBnodes];
                D.f[DIR_MMP] = &DDStart[DIR_PPM * numberOfLBnodes];
                D.f[DIR_PPP] = &DDStart[DIR_MMM * numberOfLBnodes];
                D.f[DIR_MPP] = &DDStart[DIR_PMM * numberOfLBnodes];
                D.f[DIR_PMP] = &DDStart[DIR_MPM * numberOfLBnodes];
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
            real mfcbb = (D.f[DIR_P00])[k];//[ke   ];// +  c2over27 ;(D.f[DIR_P00])[k  ];//ke
            real mfabb = (D.f[DIR_M00])[kw];//[kw   ];// +  c2over27 ;(D.f[DIR_M00])[kw ];
            real mfbcb = (D.f[DIR_0P0])[k];//[kn   ];// +  c2over27 ;(D.f[DIR_0P0])[k  ];//kn
            real mfbab = (D.f[DIR_0M0])[ks];//[ks   ];// +  c2over27 ;(D.f[DIR_0M0])[ks ];
            real mfbbc = (D.f[DIR_00P])[k];//[kt   ];// +  c2over27 ;(D.f[DIR_00P])[k  ];//kt
            real mfbba = (D.f[DIR_00M])[kb];//[kb   ];// +  c2over27 ;(D.f[DIR_00M])[kb ];
            real mfccb = (D.f[DIR_PP0])[k];//[kne  ];// +  c1over54 ;(D.f[DIR_PP0])[k  ];//kne
            real mfaab = (D.f[DIR_MM0])[ksw];//[ksw  ];// +  c1over54 ;(D.f[DIR_MM0])[ksw];
            real mfcab = (D.f[DIR_PM0])[ks];//[kse  ];// +  c1over54 ;(D.f[DIR_PM0])[ks ];//kse
            real mfacb = (D.f[DIR_MP0])[kw];//[knw  ];// +  c1over54 ;(D.f[DIR_MP0])[kw ];//knw
            real mfcbc = (D.f[DIR_P0P])[k];//[kte  ];// +  c1over54 ;(D.f[DIR_P0P])[k  ];//kte
            real mfaba = (D.f[DIR_M0M])[kbw];//[kbw  ];// +  c1over54 ;(D.f[DIR_M0M])[kbw];
            real mfcba = (D.f[DIR_P0M])[kb];//[kbe  ];// +  c1over54 ;(D.f[DIR_P0M])[kb ];//kbe
            real mfabc = (D.f[DIR_M0P])[kw];//[ktw  ];// +  c1over54 ;(D.f[DIR_M0P])[kw ];//ktw
            real mfbcc = (D.f[DIR_0PP])[k];//[ktn  ];// +  c1over54 ;(D.f[DIR_0PP])[k  ];//ktn
            real mfbaa = (D.f[DIR_0MM])[kbs];//[kbs  ];// +  c1over54 ;(D.f[DIR_0MM])[kbs];
            real mfbca = (D.f[DIR_0PM])[kb];//[kbn  ];// +  c1over54 ;(D.f[DIR_0PM])[kb ];//kbn
            real mfbac = (D.f[DIR_0MP])[ks];//[kts  ];// +  c1over54 ;(D.f[DIR_0MP])[ks ];//kts
            real mfbbb = (D.f[DIR_000])[k];//[kzero];// +  c8over27 ;(D.f[DIR_000])[k  ];//kzero
            real mfccc = (D.f[DIR_PPP])[k];//[ktne ];// +  c1over216;(D.f[DIR_PPP])[k  ];//ktne
            real mfaac = (D.f[DIR_MMP])[ksw];//[ktsw ];// +  c1over216;(D.f[DIR_MMP])[ksw];//ktsw
            real mfcac = (D.f[DIR_PMP])[ks];//[ktse ];// +  c1over216;(D.f[DIR_PMP])[ks ];//ktse
            real mfacc = (D.f[DIR_MPP])[kw];//[ktnw ];// +  c1over216;(D.f[DIR_MPP])[kw ];//ktnw
            real mfcca = (D.f[DIR_PPM])[kb];//[kbne ];// +  c1over216;(D.f[DIR_PPM])[kb ];//kbne
            real mfaaa = (D.f[DIR_MMM])[kbsw];//[kbsw ];// +  c1over216;(D.f[DIR_MMM])[kbsw];
            real mfcaa = (D.f[DIR_PMM])[kbs];//[kbse ];// +  c1over216;(D.f[DIR_PMM])[kbs];//kbse
            real mfaca = (D.f[DIR_MPM])[kbw];//[kbnw ];// +  c1over216;(D.f[DIR_MPM])[kbw];//kbnw
                                               ////////////////////////////////////////////////////////////////////////////////////
            real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

            real rho = c1o1 + drho;
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
            //real oMdrho = c1o1; // comp special
            //real m0, m1, m2;
            real vx2;
            real vy2;
            real vz2;
            vx2 = vvx*vvx;
            vy2 = vvy*vvy;
            vz2 = vvz*vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            //real wadjust;
            //real qudricLimitP = c1o100;// * 0.0001f;
            //real qudricLimitM = c1o100;// * 0.0001f;
            //real qudricLimitD = c1o100;// * 0.001f;
            //real s9 = minusomega;
            //test
            //s9 = 0.;


            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            real EQcbb = c0o1;
            real EQabb = c0o1;
            real EQbcb = c0o1;
            real EQbab = c0o1;
            real EQbbc = c0o1;
            real EQbba = c0o1;
            real EQccb = c0o1;
            real EQaab = c0o1;
            real EQcab = c0o1;
            real EQacb = c0o1;
            real EQcbc = c0o1;
            real EQaba = c0o1;
            real EQcba = c0o1;
            real EQabc = c0o1;
            real EQbcc = c0o1;
            real EQbaa = c0o1;
            real EQbca = c0o1;
            real EQbac = c0o1;
            real EQbbb = c0o1;
            real EQccc = drho * c1o27;
            real EQaac = drho * c1o3;
            real EQcac = drho * c1o9;
            real EQacc = drho * c1o9;
            real EQcca = drho * c1o9;
            real EQaaa = drho;
            real EQcaa = drho * c1o3;
            real EQaca = drho * c1o3;
            ////////////////////////////////////////////////////////////////////////////////////
            vf::lbm::backwardChimeraWithK(EQaaa, EQaab, EQaac, vvz, vz2, c1o1);
            vf::lbm::backwardChimeraWithK(EQaca, EQacb, EQacc, vvz, vz2, c1o3);
            ///////////////////////////////////////////////////////////
            EQcaa = EQaca; EQcab = EQacb; EQcac = EQacc;
            ///////////////////////////////////////////////////////////
            vf::lbm::backwardChimeraWithK(EQcca, EQccb, EQccc, vvz, vz2, c1o9);

            vf::lbm::backwardChimeraWithK(EQaaa, EQaba, EQaca, vvy, vy2, c1o6);
            vf::lbm::backwardChimeraWithK(EQaab, EQabb, EQacb, vvy, vy2, c2o3);
            vf::lbm::backwardChimeraWithK(EQaac, EQabc, EQacc, vvy, vy2, c1o6);
            vf::lbm::backwardChimeraWithK(EQcaa, EQcba, EQcca, vvy, vy2, c1o18);
            vf::lbm::backwardChimeraWithK(EQcab, EQcbb, EQccb, vvy, vy2, c2o9);
            vf::lbm::backwardChimeraWithK(EQcac, EQcbc, EQccc, vvy, vy2, c1o18);

            vf::lbm::backwardChimeraWithK(EQaaa, EQbaa, EQcaa, vvx, vx2, c1o36);
            vf::lbm::backwardChimeraWithK(EQaab, EQbab, EQcab, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQaac, EQbac, EQcac, vvx, vx2, c1o36);
            vf::lbm::backwardChimeraWithK(EQaba, EQbba, EQcba, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQabb, EQbbb, EQcbb, vvx, vx2, c4o9);
            vf::lbm::backwardChimeraWithK(EQabc, EQbbc, EQcbc, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQaca, EQbca, EQcca, vvx, vx2, c1o36);
            vf::lbm::backwardChimeraWithK(EQacb, EQbcb, EQccb, vvx, vx2, c1o9);
            vf::lbm::backwardChimeraWithK(EQacc, EQbcc, EQccc, vvx, vx2, c1o36);

            ////////////////////////////////////////////////////////////////////////////////////
            //Pre-condition
            mfcbb -= EQcbb;
            mfabb -= EQabb;
            mfbcb -= EQbcb;
            mfbab -= EQbab;
            mfbbc -= EQbbc;
            mfbba -= EQbba;
            mfccb -= EQccb;
            mfaab -= EQaab;
            mfcab -= EQcab;
            mfacb -= EQacb;
            mfcbc -= EQcbc;
            mfaba -= EQaba;
            mfcba -= EQcba;
            mfabc -= EQabc;
            mfbcc -= EQbcc;
            mfbaa -= EQbaa;
            mfbca -= EQbca;
            mfbac -= EQbac;
            mfbbb -= EQbbb;
            mfccc -= EQccc;
            mfaac -= EQaac;
            mfcac -= EQcac;
            mfacc -= EQacc;
            mfcca -= EQcca;
            mfaaa -= EQaaa;
            mfcaa -= EQcaa;
            mfaca -= EQaca;

            ////////////////////////////////////////////////////////////////////////////////////
            //Hin
            ////////////////////////////////////////////////////////////////////////////////////
            vf::lbm::forwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
            vf::lbm::forwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
            vf::lbm::forwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
            vf::lbm::forwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
            vf::lbm::forwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
            vf::lbm::forwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
            vf::lbm::forwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
            vf::lbm::forwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
            vf::lbm::forwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

            vf::lbm::forwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
            vf::lbm::forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
            vf::lbm::forwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
            vf::lbm::forwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
            vf::lbm::forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
            vf::lbm::forwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
            vf::lbm::forwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
            vf::lbm::forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
            vf::lbm::forwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

            vf::lbm::forwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
            vf::lbm::forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
            vf::lbm::forwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
            vf::lbm::forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
            vf::lbm::forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
            vf::lbm::forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
            vf::lbm::forwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
            vf::lbm::forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
            vf::lbm::forwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

            //////////////////////////////////////////////////////////////////////////////////////
            ////Hin
            //////////////////////////////////////////////////////////////////////////////////////
            //// mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Z - Dir
            //forwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, c4o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //// mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Y - Dir
            //forwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o18);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, c2o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, c2o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //// mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// X - Dir
            //forwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, one);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
            //////////////////////////////////////////////////////////////////////////////////////
            //forwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            real OxxPyyPzz = c1o1; //omega; // one;	//set the bulk viscosity one is high / two is very low and zero is (too) high

            ////////////////////////////////////////////////////////////
            //3.
            //////////////////////////////
            real OxyyPxzz = c1o1;
            real OxyyMxzz = c1o1;
            //real Oxyz = c1o1;
            ////////////////////////////////////////////////////////////
            //4.
            //////////////////////////////
            real O4 = c1o1;
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
            real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
            real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
            real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

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




            //2.
            // linear combinations
            real mxxPyyPzz = mfcaa + mfaca + mfaac;
            real mxxMyy = mfcaa - mfaca;
            real mxxMzz = mfcaa - mfaac;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
            {
                real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
                real dyuy = dxux + omega * c3o2 * mxxMyy;
                real dzuz = dxux + omega * c3o2 * mxxMzz;

                //relax
                mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz) - c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
                mxxMyy += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
                mxxMzz += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////no correction
            //mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);//-magicBulk*OxxPyyPzz;
            //mxxMyy += -(-omega) * (-mxxMyy);
            //mxxMzz += -(-omega) * (-mxxMzz);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            mfabb += omega * (-mfabb);
            mfbab += omega * (-mfbab);
            mfbba += omega * (-mfbba);

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
            mfbbb += OxyyMxzz * (-mfbbb);
            mxxyPyzz += OxyyPxzz * (-mxxyPyzz);
            mxxyMyzz += OxyyMxzz * (-mxxyMyzz);
            mxxzPyyz += OxyyPxzz * (-mxxzPyyz);
            mxxzMyyz += OxyyMxzz * (-mxxzMyyz);
            mxyyPxzz += OxyyPxzz * (-mxyyPxzz);
            mxyyMxzz += OxyyMxzz * (-mxyyMxzz);
            //////////////////////////////////////////////////////////////////////////

            mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
            mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
            mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
            mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
            mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
            mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

            //4.
            //////////////////////////////////////////////////////////////////////////
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
            mfcbb = CUMcbb + ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) / rho;
            mfbcb = CUMbcb + ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) / rho;
            mfbbc = CUMbbc + ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) / rho;

            mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) / rho - c1o9*(drho / rho));
            mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) / rho - c1o9*(drho / rho));
            mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) / rho - c1o9*(drho / rho));

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

            ////////////////////////////////////////////////////////////////////////////////////
            //the force be with you
            mfbaa = -mfbaa;
            mfaba = -mfaba;
            mfaab = -mfaab;
            ////////////////////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////////////////////
            //back
            ////////////////////////////////////////////////////////////////////////////////////
            vf::lbm::backwardChimera(mfaaa, mfaab, mfaac, vvz, vz2);
            vf::lbm::backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
            vf::lbm::backwardChimera(mfaca, mfacb, mfacc, vvz, vz2);
            vf::lbm::backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
            vf::lbm::backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
            vf::lbm::backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
            vf::lbm::backwardChimera(mfcaa, mfcab, mfcac, vvz, vz2);
            vf::lbm::backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
            vf::lbm::backwardChimera(mfcca, mfccb, mfccc, vvz, vz2);

            vf::lbm::backwardChimera(mfaaa, mfaba, mfaca, vvy, vy2);
            vf::lbm::backwardChimera(mfaab, mfabb, mfacb, vvy, vy2);
            vf::lbm::backwardChimera(mfaac, mfabc, mfacc, vvy, vy2);
            vf::lbm::backwardChimera(mfbaa, mfbba, mfbca, vvy, vy2);
            vf::lbm::backwardChimera(mfbab, mfbbb, mfbcb, vvy, vy2);
            vf::lbm::backwardChimera(mfbac, mfbbc, mfbcc, vvy, vy2);
            vf::lbm::backwardChimera(mfcaa, mfcba, mfcca, vvy, vy2);
            vf::lbm::backwardChimera(mfcab, mfcbb, mfccb, vvy, vy2);
            vf::lbm::backwardChimera(mfcac, mfcbc, mfccc, vvy, vy2);

            vf::lbm::backwardChimera(mfaaa, mfbaa, mfcaa, vvx, vx2);
            vf::lbm::backwardChimera(mfaab, mfbab, mfcab, vvx, vx2);
            vf::lbm::backwardChimera(mfaac, mfbac, mfcac, vvx, vx2);
            vf::lbm::backwardChimera(mfaba, mfbba, mfcba, vvx, vx2);
            vf::lbm::backwardChimera(mfabb, mfbbb, mfcbb, vvx, vx2);
            vf::lbm::backwardChimera(mfabc, mfbbc, mfcbc, vvx, vx2);
            vf::lbm::backwardChimera(mfaca, mfbca, mfcca, vvx, vx2);
            vf::lbm::backwardChimera(mfacb, mfbcb, mfccb, vvx, vx2);
            vf::lbm::backwardChimera(mfacc, mfbcc, mfccc, vvx, vx2);

            ////////////////////////////////////////////////////////////////////////////////////
            mfcbb+=EQcbb;
            mfabb+=EQabb;
            mfbcb+=EQbcb;
            mfbab+=EQbab;
            mfbbc+=EQbbc;
            mfbba+=EQbba;
            mfccb+=EQccb;
            mfaab+=EQaab;
            mfcab+=EQcab;
            mfacb+=EQacb;
            mfcbc+=EQcbc;
            mfaba+=EQaba;
            mfcba+=EQcba;
            mfabc+=EQabc;
            mfbcc+=EQbcc;
            mfbaa+=EQbaa;
            mfbca+=EQbca;
            mfbac+=EQbac;
            mfbbb+=EQbbb;
            mfccc+=EQccc;
            mfaac+=EQaac;
            mfcac+=EQcac;
            mfacc+=EQacc;
            mfcca+=EQcca;
            mfaaa+=EQaaa;
            mfcaa+=EQcaa;
            mfaca+=EQaca;


            //////////////////////////////////////////////////////////////////////////////////////
            ////back
            //////////////////////////////////////////////////////////////////////////////////////
            ////mit 1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9   Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Z - Dir
            //backwardChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, one);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfaba, mfabb, mfabc, vvz, vz2);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, c1o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbaa, mfbab, mfbac, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbba, mfbbb, mfbbc, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbca, mfbcb, mfbcc, vvz, vz2);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, c1o3);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfcba, mfcbb, mfcbc, vvz, vz2);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            ////mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// Y - Dir
            //backwardChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaab, mfabb, mfacb, vvy, vy2, c2o3);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, c1o6);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbaa, mfbba, mfbca, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbab, mfbbb, mfbcb, vvz, vz2);
            ///////////b//////////////////////////////////////////////////////////////////////////
            //backwardChimera(mfbac, mfbbc, mfbcc, vvz, vz2);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, c1o18);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcab, mfcbb, mfccb, vvy, vy2, c2o9);
            ///////////c//////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, c1o18);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            ////mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
            //////////////////////////////////////////////////////////////////////////////////////
            //// X - Dir
            //backwardChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaba, mfbba, mfcba, vvx, vx2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaab, mfbab, mfcab, vvx, vx2, c1o9);
            /////////////b////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfabb, mfbbb, mfcbb, vvx, vx2, c4o9);
            /////////////b////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfacb, mfbcb, mfccb, vvx, vx2, c1o9);
            //////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, c1o36);
            /////////////c////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfabc, mfbbc, mfcbc, vvx, vx2, c1o9);
            /////////////c////////////////////////////////////////////////////////////////////////
            //backwardChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, c1o36);
            //////////////////////////////////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////////////////////////////////
            real drhoPost =
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                    ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
            mfbbb += drho - drhoPost;
            ////////////////////////////////////////////////////////////////////////////////////
            (D.f[DIR_P00])[k] = mfabb;//(D.f[ DIR_P00   ])[ke   ] = mfabb;// -  c2over27 ;  (D.f[ DIR_P00   ])[k   ]                                                                     
            (D.f[DIR_M00])[kw] = mfcbb;//(D.f[ DIR_M00   ])[kw   ] = mfcbb;// -  c2over27 ;  (D.f[ DIR_M00   ])[kw  ]                                                                   
            (D.f[DIR_0P0])[k] = mfbab;//(D.f[ DIR_0P0   ])[kn   ] = mfbab;// -  c2over27 ;	 (D.f[ DIR_0P0   ])[k   ]
            (D.f[DIR_0M0])[ks] = mfbcb;//(D.f[ DIR_0M0   ])[ks   ] = mfbcb;// -  c2over27 ;	 (D.f[ DIR_0M0   ])[ks  ]
            (D.f[DIR_00P])[k] = mfbba;//(D.f[ DIR_00P   ])[kt   ] = mfbba;// -  c2over27 ;	 (D.f[ DIR_00P   ])[k   ]
            (D.f[DIR_00M])[kb] = mfbbc;//(D.f[ DIR_00M   ])[kb   ] = mfbbc;// -  c2over27 ;	 (D.f[ DIR_00M   ])[kb  ]
            (D.f[DIR_PP0])[k] = mfaab;//(D.f[ DIR_PP0  ])[kne  ] = mfaab;// -  c1over54 ;	 (D.f[ DIR_PP0  ])[k   ]
            (D.f[DIR_MM0])[ksw] = mfccb;//(D.f[ DIR_MM0  ])[ksw  ] = mfccb;// -  c1over54 ;	 (D.f[ DIR_MM0  ])[ksw ]
            (D.f[DIR_PM0])[ks] = mfacb;//(D.f[ DIR_PM0  ])[kse  ] = mfacb;// -  c1over54 ;	 (D.f[ DIR_PM0  ])[ks  ]
            (D.f[DIR_MP0])[kw] = mfcab;//(D.f[ DIR_MP0  ])[knw  ] = mfcab;// -  c1over54 ;	 (D.f[ DIR_MP0  ])[kw  ]
            (D.f[DIR_P0P])[k] = mfaba;//(D.f[ DIR_P0P  ])[kte  ] = mfaba;// -  c1over54 ;	 (D.f[ DIR_P0P  ])[k   ]
            (D.f[DIR_M0M])[kbw] = mfcbc;//(D.f[ DIR_M0M  ])[kbw  ] = mfcbc;// -  c1over54 ;	 (D.f[ DIR_M0M  ])[kbw ]
            (D.f[DIR_P0M])[kb] = mfabc;//(D.f[ DIR_P0M  ])[kbe  ] = mfabc;// -  c1over54 ;	 (D.f[ DIR_P0M  ])[kb  ]
            (D.f[DIR_M0P])[kw] = mfcba;//(D.f[ DIR_M0P  ])[ktw  ] = mfcba;// -  c1over54 ;	 (D.f[ DIR_M0P  ])[kw  ]
            (D.f[DIR_0PP])[k] = mfbaa;//(D.f[ DIR_0PP  ])[ktn  ] = mfbaa;// -  c1over54 ;	 (D.f[ DIR_0PP  ])[k   ]
            (D.f[DIR_0MM])[kbs] = mfbcc;//(D.f[ DIR_0MM  ])[kbs  ] = mfbcc;// -  c1over54 ;	 (D.f[ DIR_0MM  ])[kbs ]
            (D.f[DIR_0PM])[kb] = mfbac;//(D.f[ DIR_0PM  ])[kbn  ] = mfbac;// -  c1over54 ;	 (D.f[ DIR_0PM  ])[kb  ]
            (D.f[DIR_0MP])[ks] = mfbca;//(D.f[ DIR_0MP  ])[kts  ] = mfbca;// -  c1over54 ;	 (D.f[ DIR_0MP  ])[ks  ]
            (D.f[DIR_000])[k] = mfbbb;//(D.f[ DIR_000])[kzero] = mfbbb;// -  c8over27 ;	 (D.f[ DIR_000])[k   ]
            (D.f[DIR_PPP])[k] = mfaaa;//(D.f[ DIR_PPP ])[ktne ] = mfaaa;// -  c1over216;	 (D.f[ DIR_PPP ])[k   ]
            (D.f[DIR_PMP])[ks] = mfaca;//(D.f[ DIR_PMP ])[ktse ] = mfaca;// -  c1over216;	 (D.f[ DIR_PMP ])[ks  ]
            (D.f[DIR_PPM])[kb] = mfaac;//(D.f[ DIR_PPM ])[kbne ] = mfaac;// -  c1over216;	 (D.f[ DIR_PPM ])[kb  ]
            (D.f[DIR_PMM])[kbs] = mfacc;//(D.f[ DIR_PMM ])[kbse ] = mfacc;// -  c1over216;	 (D.f[ DIR_PMM ])[kbs ]
            (D.f[DIR_MPP])[kw] = mfcaa;//(D.f[ DIR_MPP ])[ktnw ] = mfcaa;// -  c1over216;	 (D.f[ DIR_MPP ])[kw  ]
            (D.f[DIR_MMP])[ksw] = mfcca;//(D.f[ DIR_MMP ])[ktsw ] = mfcca;// -  c1over216;	 (D.f[ DIR_MMP ])[ksw ]
            (D.f[DIR_MPM])[kbw] = mfcac;//(D.f[ DIR_MPM ])[kbnw ] = mfcac;// -  c1over216;	 (D.f[ DIR_MPM ])[kbw ]
            (D.f[DIR_MMM])[kbsw] = mfccc;//(D.f[ DIR_MMM ])[kbsw ] = mfccc;// -  c1over216;	 (D.f[ DIR_MMM ])[kbsw]
            ////////////////////////////////////////////////////////////////////////////////////
        }
    }
}
////////////////////////////////////////////////////////////////////////////////








































////////////////////////////////////////////////////////////////////////////////
__global__ void Cumulant_One_chim_Comp_SP_27(
    real omega,
    unsigned int* bcMatD,
    unsigned int* neighborX,
    unsigned int* neighborY,
    unsigned int* neighborZ,
    real* DDStart,
    unsigned long long numberOfLBnodes,
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

    if (k<numberOfLBnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        unsigned int BC;
        BC = bcMatD[k];

        if (BC >= GEO_FLUID/*(BC != GEO_SOLID) && (BC != GEO_VOID)*/)
        {
            Distributions27 D;
            if (EvenOrOdd == true)
            {
                D.f[DIR_P00] = &DDStart[DIR_P00 * numberOfLBnodes];
                D.f[DIR_M00] = &DDStart[DIR_M00 * numberOfLBnodes];
                D.f[DIR_0P0] = &DDStart[DIR_0P0 * numberOfLBnodes];
                D.f[DIR_0M0] = &DDStart[DIR_0M0 * numberOfLBnodes];
                D.f[DIR_00P] = &DDStart[DIR_00P * numberOfLBnodes];
                D.f[DIR_00M] = &DDStart[DIR_00M * numberOfLBnodes];
                D.f[DIR_PP0] = &DDStart[DIR_PP0 * numberOfLBnodes];
                D.f[DIR_MM0] = &DDStart[DIR_MM0 * numberOfLBnodes];
                D.f[DIR_PM0] = &DDStart[DIR_PM0 * numberOfLBnodes];
                D.f[DIR_MP0] = &DDStart[DIR_MP0 * numberOfLBnodes];
                D.f[DIR_P0P] = &DDStart[DIR_P0P * numberOfLBnodes];
                D.f[DIR_M0M] = &DDStart[DIR_M0M * numberOfLBnodes];
                D.f[DIR_P0M] = &DDStart[DIR_P0M * numberOfLBnodes];
                D.f[DIR_M0P] = &DDStart[DIR_M0P * numberOfLBnodes];
                D.f[DIR_0PP] = &DDStart[DIR_0PP * numberOfLBnodes];
                D.f[DIR_0MM] = &DDStart[DIR_0MM * numberOfLBnodes];
                D.f[DIR_0PM] = &DDStart[DIR_0PM * numberOfLBnodes];
                D.f[DIR_0MP] = &DDStart[DIR_0MP * numberOfLBnodes];
                D.f[DIR_000] = &DDStart[DIR_000 * numberOfLBnodes];
                D.f[DIR_PPP] = &DDStart[DIR_PPP * numberOfLBnodes];
                D.f[DIR_MMP] = &DDStart[DIR_MMP * numberOfLBnodes];
                D.f[DIR_PMP] = &DDStart[DIR_PMP * numberOfLBnodes];
                D.f[DIR_MPP] = &DDStart[DIR_MPP * numberOfLBnodes];
                D.f[DIR_PPM] = &DDStart[DIR_PPM * numberOfLBnodes];
                D.f[DIR_MMM] = &DDStart[DIR_MMM * numberOfLBnodes];
                D.f[DIR_PMM] = &DDStart[DIR_PMM * numberOfLBnodes];
                D.f[DIR_MPM] = &DDStart[DIR_MPM * numberOfLBnodes];
            }
            else
            {
                D.f[DIR_M00] = &DDStart[DIR_P00 * numberOfLBnodes];
                D.f[DIR_P00] = &DDStart[DIR_M00 * numberOfLBnodes];
                D.f[DIR_0M0] = &DDStart[DIR_0P0 * numberOfLBnodes];
                D.f[DIR_0P0] = &DDStart[DIR_0M0 * numberOfLBnodes];
                D.f[DIR_00M] = &DDStart[DIR_00P * numberOfLBnodes];
                D.f[DIR_00P] = &DDStart[DIR_00M * numberOfLBnodes];
                D.f[DIR_MM0] = &DDStart[DIR_PP0 * numberOfLBnodes];
                D.f[DIR_PP0] = &DDStart[DIR_MM0 * numberOfLBnodes];
                D.f[DIR_MP0] = &DDStart[DIR_PM0 * numberOfLBnodes];
                D.f[DIR_PM0] = &DDStart[DIR_MP0 * numberOfLBnodes];
                D.f[DIR_M0M] = &DDStart[DIR_P0P * numberOfLBnodes];
                D.f[DIR_P0P] = &DDStart[DIR_M0M * numberOfLBnodes];
                D.f[DIR_M0P] = &DDStart[DIR_P0M * numberOfLBnodes];
                D.f[DIR_P0M] = &DDStart[DIR_M0P * numberOfLBnodes];
                D.f[DIR_0MM] = &DDStart[DIR_0PP * numberOfLBnodes];
                D.f[DIR_0PP] = &DDStart[DIR_0MM * numberOfLBnodes];
                D.f[DIR_0MP] = &DDStart[DIR_0PM * numberOfLBnodes];
                D.f[DIR_0PM] = &DDStart[DIR_0MP * numberOfLBnodes];
                D.f[DIR_000] = &DDStart[DIR_000 * numberOfLBnodes];
                D.f[DIR_MMM] = &DDStart[DIR_PPP * numberOfLBnodes];
                D.f[DIR_PPM] = &DDStart[DIR_MMP * numberOfLBnodes];
                D.f[DIR_MPM] = &DDStart[DIR_PMP * numberOfLBnodes];
                D.f[DIR_PMM] = &DDStart[DIR_MPP * numberOfLBnodes];
                D.f[DIR_MMP] = &DDStart[DIR_PPM * numberOfLBnodes];
                D.f[DIR_PPP] = &DDStart[DIR_MMM * numberOfLBnodes];
                D.f[DIR_MPP] = &DDStart[DIR_PMM * numberOfLBnodes];
                D.f[DIR_PMP] = &DDStart[DIR_MPM * numberOfLBnodes];
            }
            ////////////////////////////////////////////////////////////////////////////////
            //index
            unsigned int kw = neighborX[k];
            unsigned int ks = neighborY[k];
            unsigned int kb = neighborZ[k];
            unsigned int ksw = neighborY[kw];
            unsigned int kbw = neighborZ[kw];
            unsigned int kbs = neighborZ[ks];
            unsigned int kbsw = neighborZ[ksw];
            ////////////////////////////////////////////////////////////////////////////////////
            real mfcbb = (D.f[DIR_P00])[k   ];
            real mfabb = (D.f[DIR_M00])[kw  ];
            real mfbcb = (D.f[DIR_0P0])[k   ];
            real mfbab = (D.f[DIR_0M0])[ks  ];
            real mfbbc = (D.f[DIR_00P])[k   ];
            real mfbba = (D.f[DIR_00M])[kb  ];
            real mfccb = (D.f[DIR_PP0])[k   ];
            real mfaab = (D.f[DIR_MM0])[ksw ];
            real mfcab = (D.f[DIR_PM0])[ks  ];
            real mfacb = (D.f[DIR_MP0])[kw  ];
            real mfcbc = (D.f[DIR_P0P])[k   ];
            real mfaba = (D.f[DIR_M0M])[kbw ];
            real mfcba = (D.f[DIR_P0M])[kb  ];
            real mfabc = (D.f[DIR_M0P])[kw  ];
            real mfbcc = (D.f[DIR_0PP])[k   ];
            real mfbaa = (D.f[DIR_0MM])[kbs ];
            real mfbca = (D.f[DIR_0PM])[kb  ];
            real mfbac = (D.f[DIR_0MP])[ks  ];
            real mfbbb = (D.f[DIR_000])[k   ];
            real mfccc = (D.f[DIR_PPP])[k   ];
            real mfaac = (D.f[DIR_MMP])[ksw ];
            real mfcac = (D.f[DIR_PMP])[ks  ];
            real mfacc = (D.f[DIR_MPP])[kw  ];
            real mfcca = (D.f[DIR_PPM])[kb  ];
            real mfaaa = (D.f[DIR_MMM])[kbsw];
            real mfcaa = (D.f[DIR_PMM])[kbs ];
            real mfaca = (D.f[DIR_MPM])[kbw ];
            ////////////////////////////////////////////////////////////////////////////////////
            real drho = ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;

            real rho = c1o1 + drho;
            real OOrho = c1o1 / rho;
            ////////////////////////////////////////////////////////////////////////////////////
            real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) +
                (mfcbb - mfabb)) * OOrho;
            real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) +
                (mfbcb - mfbab)) * OOrho;
            real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) +
                (mfbbc - mfbba)) * OOrho;
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
            //real oMdrho = c1o1; // comp special
            //real m0, m1, m2;
            real vx2;
            real vy2;
            real vz2;
            vx2 = vvx*vvx;
            vy2 = vvy*vvy;
            vz2 = vvz*vvz;
            ////////////////////////////////////////////////////////////////////////////////////
            //real wadjust;
            //real qudricLimitP = c1o100;// * 0.0001f;
            //real qudricLimitM = c1o100;// * 0.0001f;
            //real qudricLimitD = c1o100;// * 0.001f;
            ////////////////////////////////////////////////////////////////////////////////////
            //Hin
            ////////////////////////////////////////////////////////////////////////////////////
            // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Z - Dir
            vf::lbm::forwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, 36.0f, c1o36);
            vf::lbm::forwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::forwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, 36.0f, c1o36);
            vf::lbm::forwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::forwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, 2.25f, c4o9 );
            vf::lbm::forwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::forwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, 36.0f, c1o36);
            vf::lbm::forwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::forwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, 36.0f, c1o36);

            ////////////////////////////////////////////////////////////////////////////////////
            // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            vf::lbm::forwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, 6.0f , c1o6 );
            vf::lbm::forwardChimera(     mfaab, mfabb, mfacb, vvy, vy2);
            vf::lbm::forwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, 18.0f, c1o18);
            vf::lbm::forwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, 1.5f , c2o3 );
            vf::lbm::forwardChimera(     mfbab, mfbbb, mfbcb, vvy, vy2);
            vf::lbm::forwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, 4.5f , c2o9 );
            vf::lbm::forwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, 6.0f , c1o6 );
            vf::lbm::forwardChimera(     mfcab, mfcbb, mfccb, vvy, vy2);
            vf::lbm::forwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, 18.0f, c1o18);

            ////////////////////////////////////////////////////////////////////////////////////
            // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9		Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // X - Dir
            vf::lbm::forwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
            vf::lbm::forwardChimera(     mfaba, mfbba, mfcba, vvx, vx2);
            vf::lbm::forwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, 3.0f, c1o3);
            vf::lbm::forwardChimera(     mfaab, mfbab, mfcab, vvx, vx2);
            vf::lbm::forwardChimera(     mfabb, mfbbb, mfcbb, vvx, vx2);
            vf::lbm::forwardChimera(     mfacb, mfbcb, mfccb, vvx, vx2);
            vf::lbm::forwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, 3.0f, c1o3);
            vf::lbm::forwardChimera(     mfabc, mfbbc, mfcbc, vvx, vx2);
            vf::lbm::forwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, 9.0f, c1o9);

            ////////////////////////////////////////////////////////////////////////////////////
            // Cumulants
            ////////////////////////////////////////////////////////////////////////////////////
            real OxxPyyPzz = c1o1;
            ////////////////////////////////////////////////////////////
            //3.
            //////////////////////////////
            real OxyyPxzz = c1o1;
            real OxyyMxzz = c1o1;
            //real Oxyz = c1o1;
            ////////////////////////////////////////////////////////////
            //4.
            //////////////////////////////
            real O4 = c1o1;
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
            real CUMcbb = mfcbb - ((mfcaa + c1o3) * mfabb + c2o1 * mfbba * mfbab) * OOrho;
            real CUMbcb = mfbcb - ((mfaca + c1o3) * mfbab + c2o1 * mfbba * mfabb) * OOrho;
            real CUMbbc = mfbbc - ((mfaac + c1o3) * mfbba + c2o1 * mfbab * mfabb) * OOrho;

            real CUMcca = mfcca - (((mfcaa * mfaca + c2o1 * mfbba * mfbba) + c1o3 * (mfcaa + mfaca)) * OOrho - c1o9*(drho * OOrho));
            real CUMcac = mfcac - (((mfcaa * mfaac + c2o1 * mfbab * mfbab) + c1o3 * (mfcaa + mfaac)) * OOrho - c1o9*(drho * OOrho));
            real CUMacc = mfacc - (((mfaac * mfaca + c2o1 * mfabb * mfabb) + c1o3 * (mfaac + mfaca)) * OOrho - c1o9*(drho * OOrho));

            //5.
            real CUMbcc = mfbcc - ((mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + c1o3 * (mfbca + mfbac)) * OOrho;
            real CUMcbc = mfcbc - ((mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + c1o3 * (mfcba + mfabc)) * OOrho;
            real CUMccb = mfccb - ((mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) + c1o3 * (mfacb + mfcab)) * OOrho;

            //6.
            real CUMccc = mfccc + ((-c4o1 *  mfbbb * mfbbb
                - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                - c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                - c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) * OOrho
                + (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                + c2o1 * (mfcaa * mfaca * mfaac)
                + c16o1 *  mfbba * mfbab * mfabb) * OOrho * OOrho
                - c1o3 * (mfacc + mfcac + mfcca) * OOrho
                - c1o9 * (mfcaa + mfaca + mfaac) * OOrho
                + (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) * OOrho * OOrho  * c2o3
                + c1o27*((drho * drho - drho) * OOrho * OOrho ));


            //2.
            // linear combinations
            real mxxPyyPzz = mfcaa + mfaca + mfaac;
            real mxxMyy = mfcaa - mfaca;
            real mxxMzz = mfcaa - mfaac;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //incl. correction		(hat noch nicht so gut funktioniert...Optimierungsbedarf??)
            {
                real dxux = c1o2 * (-omega) *(mxxMyy + mxxMzz) + c1o2 *  OxxPyyPzz * (mfaaa - mxxPyyPzz);
                real dyuy = dxux + omega * c3o2 * mxxMyy;
                real dzuz = dxux + omega * c3o2 * mxxMzz;

                //relax
                mxxPyyPzz += OxxPyyPzz*(mfaaa  - mxxPyyPzz)- c3o1 * (c1o1 - c1o2 * OxxPyyPzz) * (vx2 * dxux + vy2 * dyuy + vz2 * dzuz);//-magicBulk*OxxPyyPzz;
                mxxMyy    += omega * (-mxxMyy) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vy2 * dyuy);
                mxxMzz    += omega * (-mxxMzz) - c3o1 * (c1o1 + c1o2 * (-omega)) * (vx2 * dxux - vz2 * dzuz);

            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////no correction
            //mxxPyyPzz += OxxPyyPzz*(mfaaa - mxxPyyPzz);//-magicBulk*OxxPyyPzz;
            //mxxMyy += -(-omega) * (-mxxMyy);
            //mxxMzz += -(-omega) * (-mxxMzz);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            mfabb += omega * (-mfabb);
            mfbab += omega * (-mfbab);
            mfbba += omega * (-mfbba);

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
            mfbbb     += OxyyMxzz * (-mfbbb);
            mxxyPyzz  += OxyyPxzz * (-mxxyPyzz);
            mxxyMyzz  += OxyyMxzz * (-mxxyMyzz);
            mxxzPyyz  += OxyyPxzz * (-mxxzPyyz);
            mxxzMyyz  += OxyyMxzz * (-mxxzMyyz);
            mxyyPxzz  += OxyyPxzz * (-mxyyPxzz);
            mxyyMxzz  += OxyyMxzz * (-mxyyMxzz);
            //////////////////////////////////////////////////////////////////////////

            mfcba = (mxxyMyzz + mxxyPyzz) * c1o2;
            mfabc = (-mxxyMyzz + mxxyPyzz) * c1o2;
            mfcab = (mxxzMyyz + mxxzPyyz) * c1o2;
            mfacb = (-mxxzMyyz + mxxzPyyz) * c1o2;
            mfbca = (mxyyMxzz + mxyyPxzz) * c1o2;
            mfbac = (-mxyyMxzz + mxyyPxzz) * c1o2;

            //4.
            //////////////////////////////////////////////////////////////////////////
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
            mfcbb = CUMcbb + c1o3*((c3o1*mfcaa + c1o1) * mfabb + c6o1 * mfbba * mfbab) * OOrho; 
            mfbcb = CUMbcb + c1o3*((c3o1*mfaca + c1o1) * mfbab + c6o1 * mfbba * mfabb) * OOrho;
            mfbbc = CUMbbc + c1o3*((c3o1*mfaac + c1o1) * mfbba + c6o1 * mfbab * mfabb) * OOrho;

            mfcca = CUMcca + (((mfcaa * mfaca + c2o1 * mfbba * mfbba)*c9o1 + c3o1 * (mfcaa + mfaca)) * OOrho - (drho * OOrho))*c1o9;
            mfcac = CUMcac + (((mfcaa * mfaac + c2o1 * mfbab * mfbab)*c9o1 + c3o1 * (mfcaa + mfaac)) * OOrho - (drho * OOrho))*c1o9;
            mfacc = CUMacc + (((mfaac * mfaca + c2o1 * mfabb * mfabb)*c9o1 + c3o1 * (mfaac + mfaca)) * OOrho - (drho * OOrho))*c1o9;

            //5.
            mfbcc = CUMbcc + c1o3 *(c3o1*(mfaac * mfbca + mfaca * mfbac + c4o1 * mfabb * mfbbb + c2o1 * (mfbab * mfacb + mfbba * mfabc)) + (mfbca + mfbac)) * OOrho;
            mfcbc = CUMcbc + c1o3 *(c3o1*(mfaac * mfcba + mfcaa * mfabc + c4o1 * mfbab * mfbbb + c2o1 * (mfabb * mfcab + mfbba * mfbac)) + (mfcba + mfabc)) * OOrho;
            mfccb = CUMccb + c1o3 *(c3o1*(mfcaa * mfacb + mfaca * mfcab + c4o1 * mfbba * mfbbb + c2o1 * (mfbab * mfbca + mfabb * mfcba)) +  (mfacb + mfcab)) * OOrho;

            //6.
            mfccc = 
                CUMccc - ((-c4o1 *  mfbbb * mfbbb
                - (mfcaa * mfacc + mfaca * mfcac + mfaac * mfcca)
                - c4o1 * (mfabb * mfcbb + mfbab * mfbcb + mfbba * mfbbc)
                - c2o1 * (mfbca * mfbac + mfcba * mfabc + mfcab * mfacb)) * OOrho
                + (c4o1 * (mfbab * mfbab * mfaca + mfabb * mfabb * mfcaa + mfbba * mfbba * mfaac)
                + c2o1 * (mfcaa * mfaca * mfaac)
                + c16o1 *  mfbba * mfbab * mfabb) * OOrho * OOrho
                - c1o3 * (mfacc + mfcac + mfcca) * OOrho
                - c1o9 * (mfcaa + mfaca + mfaac) * OOrho
                + (c2o1 * (mfbab * mfbab + mfabb * mfabb + mfbba * mfbba)
                + (mfaac * mfaca + mfaac * mfcaa + mfaca * mfcaa) + c1o3 *(mfaac + mfaca + mfcaa)) * OOrho * OOrho * c2o3
                + c1o27*((drho * drho - drho) * OOrho * OOrho ));

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
            // X - Dir
            vf::lbm::backwardInverseChimeraWithK(mfaaa, mfbaa, mfcaa, vvx, vx2, c1o1, c1o1);
            vf::lbm::backwardChimera(			mfaba, mfbba, mfcba, vvx, vx2);
            vf::lbm::backwardInverseChimeraWithK(mfaca, mfbca, mfcca, vvx, vx2, 3.0f, c1o3);
            vf::lbm::backwardChimera(			mfaab, mfbab, mfcab, vvx, vx2);
            vf::lbm::backwardChimera(			mfabb, mfbbb, mfcbb, vvx, vx2);
            vf::lbm::backwardChimera(			mfacb, mfbcb, mfccb, vvx, vx2);
            vf::lbm::backwardInverseChimeraWithK(mfaac, mfbac, mfcac, vvx, vx2, 3.0f, c1o3);
            vf::lbm::backwardChimera(			mfabc, mfbbc, mfcbc, vvx, vx2);
            vf::lbm::backwardInverseChimeraWithK(mfacc, mfbcc, mfccc, vvx, vx2, 9.0f, c1o9);

            ////////////////////////////////////////////////////////////////////////////////////
            //mit 1/6, 2/3, 1/6, 0, 0, 0, 1/18, 2/9, 1/18   Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Y - Dir
            vf::lbm::backwardInverseChimeraWithK(mfaaa, mfaba, mfaca, vvy, vy2, 6.0f , c1o6 );
            vf::lbm::backwardChimera(			mfaab, mfabb, mfacb, vvy, vy2);
            vf::lbm::backwardInverseChimeraWithK(mfaac, mfabc, mfacc, vvy, vy2, 18.0f, c1o18);
            vf::lbm::backwardInverseChimeraWithK(mfbaa, mfbba, mfbca, vvy, vy2, 1.5f , c2o3 );
            vf::lbm::backwardChimera(			mfbab, mfbbb, mfbcb, vvy, vy2);
            vf::lbm::backwardInverseChimeraWithK(mfbac, mfbbc, mfbcc, vvy, vy2, 4.5f , c2o9 );
            vf::lbm::backwardInverseChimeraWithK(mfcaa, mfcba, mfcca, vvy, vy2, 6.0f , c1o6 );
            vf::lbm::backwardChimera(			mfcab, mfcbb, mfccb, vvy, vy2);
            vf::lbm::backwardInverseChimeraWithK(mfcac, mfcbc, mfccc, vvy, vy2, 18.0f, c1o18);

            ////////////////////////////////////////////////////////////////////////////////////
            //mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36 Konditionieren
            ////////////////////////////////////////////////////////////////////////////////////
            // Z - Dir
            vf::lbm::backwardInverseChimeraWithK(mfaaa, mfaab, mfaac, vvz, vz2, 36.0f, c1o36);
            vf::lbm::backwardInverseChimeraWithK(mfaba, mfabb, mfabc, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::backwardInverseChimeraWithK(mfaca, mfacb, mfacc, vvz, vz2, 36.0f, c1o36);
            vf::lbm::backwardInverseChimeraWithK(mfbaa, mfbab, mfbac, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::backwardInverseChimeraWithK(mfbba, mfbbb, mfbbc, vvz, vz2, 2.25f, c4o9 );
            vf::lbm::backwardInverseChimeraWithK(mfbca, mfbcb, mfbcc, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::backwardInverseChimeraWithK(mfcaa, mfcab, mfcac, vvz, vz2, 36.0f, c1o36);
            vf::lbm::backwardInverseChimeraWithK(mfcba, mfcbb, mfcbc, vvz, vz2, 9.0f , c1o9 );
            vf::lbm::backwardInverseChimeraWithK(mfcca, mfccb, mfccc, vvz, vz2, 36.0f, c1o36);

            //////////////////////////////////////////////////////////////////////////////////////
            real drhoPost =
                ((((mfccc + mfaaa) + (mfaca + mfcac)) + ((mfacc + mfcaa) + (mfaac + mfcca))) +
                (((mfbac + mfbca) + (mfbaa + mfbcc)) + ((mfabc + mfcba) + (mfaba + mfcbc)) + ((mfacb + mfcab) + (mfaab + mfccb))) +
                    ((mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc))) + mfbbb;
            mfbbb += drho - drhoPost;
            ////////////////////////////////////////////////////////////////////////////////////
            (D.f[DIR_P00])[k   ] = mfabb;                                                                   
            (D.f[DIR_M00])[kw  ] = mfcbb;                                                                 
            (D.f[DIR_0P0])[k   ] = mfbab;
            (D.f[DIR_0M0])[ks  ] = mfbcb;
            (D.f[DIR_00P])[k   ] = mfbba;
            (D.f[DIR_00M])[kb  ] = mfbbc;
            (D.f[DIR_PP0])[k   ] = mfaab;
            (D.f[DIR_MM0])[ksw ] = mfccb;
            (D.f[DIR_PM0])[ks  ] = mfacb;
            (D.f[DIR_MP0])[kw  ] = mfcab;
            (D.f[DIR_P0P])[k   ] = mfaba;
            (D.f[DIR_M0M])[kbw ] = mfcbc;
            (D.f[DIR_P0M])[kb  ] = mfabc;
            (D.f[DIR_M0P])[kw  ] = mfcba;
            (D.f[DIR_0PP])[k   ] = mfbaa;
            (D.f[DIR_0MM])[kbs ] = mfbcc;
            (D.f[DIR_0PM])[kb  ] = mfbac;
            (D.f[DIR_0MP])[ks  ] = mfbca;
            (D.f[DIR_000])[k   ] = mfbbb;
            (D.f[DIR_PPP])[k   ] = mfaaa;
            (D.f[DIR_PMP])[ks  ] = mfaca;
            (D.f[DIR_PPM])[kb  ] = mfaac;
            (D.f[DIR_PMM])[kbs ] = mfacc;
            (D.f[DIR_MPP])[kw  ] = mfcaa;
            (D.f[DIR_MMP])[ksw ] = mfcca;
            (D.f[DIR_MPM])[kbw ] = mfcac;
            (D.f[DIR_MMM])[kbsw] = mfccc;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////








































