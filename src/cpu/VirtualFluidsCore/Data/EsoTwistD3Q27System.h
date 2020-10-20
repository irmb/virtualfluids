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

    const static int STARTDIR = D3Q27System::STARTDIR;
    const static int ENDDIR   = D3Q27System::ENDDIR;

    static const int ZERO = D3Q27System::ZERO; /*f0 */
    static const int E    = D3Q27System::E;    /*f1 */
    static const int W    = D3Q27System::W;    /*f2 */
    static const int N    = D3Q27System::N;    /*f3 */
    static const int S    = D3Q27System::S;    /*f4 */
    static const int T    = D3Q27System::T;    /*f5 */
    static const int B    = D3Q27System::B;    /*f6 */
    static const int NE   = D3Q27System::NE;   /*f7 */
    static const int SW   = D3Q27System::SW;   /*f8 */
    static const int SE   = D3Q27System::SE;   /*f9 */
    static const int NW   = D3Q27System::NW;   /*f10*/
    static const int TE   = D3Q27System::TE;   /*f11*/
    static const int BW   = D3Q27System::BW;   /*f12*/
    static const int BE   = D3Q27System::BE;   /*f13*/
    static const int TW   = D3Q27System::TW;   /*f14*/
    static const int TN   = D3Q27System::TN;   /*f15*/
    static const int BS   = D3Q27System::BS;   /*f16*/
    static const int BN   = D3Q27System::BN;   /*f17*/
    static const int TS   = D3Q27System::TS;   /*f18*/
    static const int TNE  = D3Q27System::TNE;
    static const int TNW  = D3Q27System::TNW;
    static const int TSE  = D3Q27System::TSE;
    static const int TSW  = D3Q27System::TSW;
    static const int BNE  = D3Q27System::BNE;
    static const int BNW  = D3Q27System::BNW;
    static const int BSE  = D3Q27System::BSE;
    static const int BSW  = D3Q27System::BSW;

    static const int INV_E   = D3Q27System::W;
    static const int INV_W   = D3Q27System::E;
    static const int INV_N   = D3Q27System::S;
    static const int INV_S   = D3Q27System::N;
    static const int INV_T   = D3Q27System::B;
    static const int INV_B   = D3Q27System::T;
    static const int INV_NE  = D3Q27System::SW;
    static const int INV_SW  = D3Q27System::NE;
    static const int INV_SE  = D3Q27System::NW;
    static const int INV_NW  = D3Q27System::SE;
    static const int INV_TE  = D3Q27System::BW;
    static const int INV_BW  = D3Q27System::TE;
    static const int INV_BE  = D3Q27System::TW;
    static const int INV_TW  = D3Q27System::BE;
    static const int INV_TN  = D3Q27System::BS;
    static const int INV_BS  = D3Q27System::TN;
    static const int INV_BN  = D3Q27System::TS;
    static const int INV_TS  = D3Q27System::BN;
    static const int INV_TNE = D3Q27System::BSW;
    static const int INV_TNW = D3Q27System::BSE;
    static const int INV_TSE = D3Q27System::BNW;
    static const int INV_TSW = D3Q27System::BNE;
    static const int INV_BNE = D3Q27System::TSW;
    static const int INV_BNW = D3Q27System::TSE;
    static const int INV_BSE = D3Q27System::TNW;
    static const int INV_BSW = D3Q27System::TNE;

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
