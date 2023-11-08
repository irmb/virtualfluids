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
//! \file EsoTwistD3Q27System.h
//! \ingroup Data
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef ESOTWISTD3Q27SYSTEM_H
#define ESOTWISTD3Q27SYSTEM_H

#include "D3Q27System.h"

//!
struct EsoTwistD3Q27System {
    const static int FSTARTDIR = D3Q27System::FSTARTDIR;
    const static int FENDDIR = D3Q27System::FENDDIR; // gellerstyle: meint alle frichtungen OHNE f0

    const static int STARTF = D3Q27System::STARTF;
    const static int ENDF   = D3Q27System::ENDF;

 //   const static int STARTDIR = D3Q27System::STARTDIR;
    const static int ENDDIR   = D3Q27System::ENDDIR;

    static const int REST = vf::lbm::dir::d000; /*f0 */
    static const int E    = vf::lbm::dir::dP00;    /*f1 */
    static const int W    = vf::lbm::dir::dM00;    /*f2 */
    static const int N    = vf::lbm::dir::d0P0;    /*f3 */
    static const int S    = vf::lbm::dir::d0M0;    /*f4 */
    static const int T    = vf::lbm::dir::d00P;    /*f5 */
    static const int B    = vf::lbm::dir::d00M;    /*f6 */
    static const int NE   = vf::lbm::dir::dPP0;   /*f7 */
    static const int SW   = vf::lbm::dir::dMM0;   /*f8 */
    static const int SE   = vf::lbm::dir::dPM0;   /*f9 */
    static const int NW   = vf::lbm::dir::dMP0;   /*f10*/
    static const int TE   = vf::lbm::dir::dP0P;   /*f11*/
    static const int BW   = vf::lbm::dir::dM0M;   /*f12*/
    static const int BE   = vf::lbm::dir::dP0M;   /*f13*/
    static const int TW   = vf::lbm::dir::dM0P;   /*f14*/
    static const int TN   = vf::lbm::dir::d0PP;   /*f15*/
    static const int BS   = vf::lbm::dir::d0MM;   /*f16*/
    static const int BN   = vf::lbm::dir::d0PM;   /*f17*/
    static const int TS   = vf::lbm::dir::d0MP;   /*f18*/
    static const int TNE  = vf::lbm::dir::dPPP;
    static const int TNW  = vf::lbm::dir::dMPP;
    static const int TSE  = vf::lbm::dir::dPMP;
    static const int TSW  = vf::lbm::dir::dMMP;
    static const int BNE  = vf::lbm::dir::dPPM;
    static const int BNW  = vf::lbm::dir::dMPM;
    static const int BSE  = vf::lbm::dir::dPMM;
    static const int BSW  = vf::lbm::dir::dMMM;

    static const int INV_E   = vf::lbm::dir::dM00;
    static const int INV_W   = vf::lbm::dir::dP00;
    static const int INV_N   = vf::lbm::dir::d0M0;
    static const int INV_S   = vf::lbm::dir::d0P0;
    static const int INV_T   = vf::lbm::dir::d00M;
    static const int INV_B   = vf::lbm::dir::d00P;
    static const int INV_NE  = vf::lbm::dir::dMM0;
    static const int INV_SW  = vf::lbm::dir::dPP0;
    static const int INV_SE  = vf::lbm::dir::dMP0;
    static const int INV_NW  = vf::lbm::dir::dPM0;
    static const int INV_TE  = vf::lbm::dir::dM0M;
    static const int INV_BW  = vf::lbm::dir::dP0P;
    static const int INV_BE  = vf::lbm::dir::dM0P;
    static const int INV_TW  = vf::lbm::dir::dP0M;
    static const int INV_TN  = vf::lbm::dir::d0MM;
    static const int INV_BS  = vf::lbm::dir::d0PP;
    static const int INV_BN  = vf::lbm::dir::d0MP;
    static const int INV_TS  = vf::lbm::dir::d0PM;
    static const int INV_TNE = vf::lbm::dir::dMMM;
    static const int INV_TNW = vf::lbm::dir::dPMM;
    static const int INV_TSE = vf::lbm::dir::dMPM;
    static const int INV_TSW = vf::lbm::dir::dPPM;
    static const int INV_BNE = vf::lbm::dir::dMMP;
    static const int INV_BNW = vf::lbm::dir::dPMP;
    static const int INV_BSE = vf::lbm::dir::dMPP;
    static const int INV_BSW = vf::lbm::dir::dPPP;

    static const unsigned long int etZERO; // 1;/*f0 */
    static const unsigned long int etE;    //  2;    /*f1 */
    static const unsigned long int etW;    //  4;    /*f2 */
    static const unsigned long int etN;    //  8;    /*f3 */
    static const unsigned long int etS;    //  16;   /*f4 */
    static const unsigned long int etT;    //  32;    /*f5 */
    static const unsigned long int etB;    //  64;   /*f6 */
    static const unsigned long int etNE;   // 128;  /*f7 */
    static const unsigned long int etSW;   // 256;  /*f8 */
    static const unsigned long int etSE;   // 512;  /*f9 */
    static const unsigned long int etNW;   // 1024;  /*f10*/
    static const unsigned long int etTE;   // 2048;  /*f11*/
    static const unsigned long int etBW;   // 4096;  /*f12*/
    static const unsigned long int etBE;   // 8192;  /*f13*/
    static const unsigned long int etTW;   // 16384;  /*f14*/
    static const unsigned long int etTN;   // 32768;  /*f15*/
    static const unsigned long int etBS;   // 65536;  /*f16*/
    static const unsigned long int etBN;   // 131072;  /*f17*/
    static const unsigned long int etTS;   // 262144;  /*f18*/
    static const unsigned long int etTNE;  // 524288;
    static const unsigned long int etTNW;  // 1048576;
    static const unsigned long int etTSE;  // 2097152;
    static const unsigned long int etTSW;  // 4194304;
    static const unsigned long int etBNE;  // 8388608;
    static const unsigned long int etBNW;  // 16777216;
    static const unsigned long int etBSE;  // 33554432;
    static const unsigned long int etBSW;  // = 67108864;

    const static int ETX1[ENDF + 1];
    const static int ETX2[ENDF + 1];
    const static int ETX3[ENDF + 1];
    const static int etINVDIR[ENDF + 1];
    const static unsigned long int etDIR[ENDF + 1];
};

#endif
