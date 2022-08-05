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
//! \file EsoTwistD3Q27System.cpp
//! \ingroup Data
//! \author Konstantin Kutscher
//=======================================================================================

#include "EsoTwistD3Q27System.h"

// index                                                              0   1   2   3   4   5  6   7   8    9  10  11  12
// 13  14  15  16  17  18  19  20  21  22  23  24  25  26 f: E,  W,  N,  S,  T,  B, NE, SW, SE, NW, TE, BW, BE, TW, TN,
// BS, BN, TS, TNE TNW TSE TSW BNE BNW BSE BSW REST
const int EsoTwistD3Q27System::ETX1[EsoTwistD3Q27System::ENDF + 1] = { 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1,
                                                                       0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0 };
const int EsoTwistD3Q27System::ETX2[EsoTwistD3Q27System::ENDF + 1] = { 0, 0, 0,  1, 0, 0,  0, 1, 0, -1, 0, 0, 0, 0,
                                                                       0, 1, -1, 0, 0, -1, 0, 1, 0, -1, 0, 1, 0 };
const int EsoTwistD3Q27System::ETX3[EsoTwistD3Q27System::ENDF + 1] = { 0, 0, 0, 0, 0, 1,  0, 0,  0, 0, 0, 1, 0, -1,
                                                                       0, 1, 1, 0, 0, -1, 0, -1, 0, 1, 0, 1, 0 };

const int EsoTwistD3Q27System::etINVDIR[EsoTwistD3Q27System::ENDF + 1] = {
    D3Q27System::INV_P00,   D3Q27System::INV_M00,   D3Q27System::INV_0P0,   D3Q27System::INV_0M0,   D3Q27System::INV_00P,
    D3Q27System::INV_00M,   D3Q27System::INV_PP0,  D3Q27System::INV_MM0,  D3Q27System::INV_PM0,  D3Q27System::INV_MP0,
    D3Q27System::INV_P0P,  D3Q27System::INV_M0M,  D3Q27System::INV_P0M,  D3Q27System::INV_M0P,  D3Q27System::INV_0PP,
    D3Q27System::INV_0MM,  D3Q27System::INV_0PM,  D3Q27System::INV_0MP,  D3Q27System::INV_PPP, D3Q27System::INV_MPP,
    D3Q27System::INV_PMP, D3Q27System::INV_MMP, D3Q27System::INV_PPM, D3Q27System::INV_MPM, D3Q27System::INV_PMM,
    D3Q27System::INV_MMM, D3Q27System::DIR_000
};

const unsigned long int EsoTwistD3Q27System::etDIR[EsoTwistD3Q27System::ENDF + 1] = {
    etE,  etW,  etN,  etS,  etT,   etB,   etNE,  etSW,  etSE,  etNW,  etTE,  etBW,  etBE,  etTW,
    etTN, etBS, etBN, etTS, etTNE, etTNW, etTSE, etTSW, etBNE, etBNW, etBSE, etBSW, etZERO
};

const unsigned long int EsoTwistD3Q27System::etZERO = 1;      /*f0 */
const unsigned long int EsoTwistD3Q27System::etE    = 2;      /*f1 */
const unsigned long int EsoTwistD3Q27System::etW    = 4;      /*f2 */
const unsigned long int EsoTwistD3Q27System::etN    = 8;      /*f3 */
const unsigned long int EsoTwistD3Q27System::etS    = 16;     /*f4 */
const unsigned long int EsoTwistD3Q27System::etT    = 32;     /*f5 */
const unsigned long int EsoTwistD3Q27System::etB    = 64;     /*f6 */
const unsigned long int EsoTwistD3Q27System::etNE   = 128;    /*f7 */
const unsigned long int EsoTwistD3Q27System::etSW   = 256;    /*f8 */
const unsigned long int EsoTwistD3Q27System::etSE   = 512;    /*f9 */
const unsigned long int EsoTwistD3Q27System::etNW   = 1024;   /*f10*/
const unsigned long int EsoTwistD3Q27System::etTE   = 2048;   /*f11*/
const unsigned long int EsoTwistD3Q27System::etBW   = 4096;   /*f12*/
const unsigned long int EsoTwistD3Q27System::etBE   = 8192;   /*f13*/
const unsigned long int EsoTwistD3Q27System::etTW   = 16384;  /*f14*/
const unsigned long int EsoTwistD3Q27System::etTN   = 32768;  /*f15*/
const unsigned long int EsoTwistD3Q27System::etBS   = 65536;  /*f16*/
const unsigned long int EsoTwistD3Q27System::etBN   = 131072; /*f17*/
const unsigned long int EsoTwistD3Q27System::etTS   = 262144; /*f18*/
const unsigned long int EsoTwistD3Q27System::etTNE  = 524288;
const unsigned long int EsoTwistD3Q27System::etTNW  = 1048576;
const unsigned long int EsoTwistD3Q27System::etTSE  = 2097152;
const unsigned long int EsoTwistD3Q27System::etTSW  = 4194304;
const unsigned long int EsoTwistD3Q27System::etBNE  = 8388608;
const unsigned long int EsoTwistD3Q27System::etBNW  = 16777216;
const unsigned long int EsoTwistD3Q27System::etBSE  = 33554432;
const unsigned long int EsoTwistD3Q27System::etBSW  = 67108864;
