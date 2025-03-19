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
//! \addtogroup constants
//! \ingroup lbm
//! \{
//=======================================================================================
#ifndef LBM_D3Q27_H
#define LBM_D3Q27_H

#include <map>

#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>

namespace vf::lbm::dir
{
static constexpr size_t STARTDIR = 0;
static constexpr size_t ENDDIR = 26;
static constexpr size_t NUMBER_Of_DIRECTIONS = ENDDIR + 1;

static constexpr size_t d000 = 0;
static constexpr size_t dP00 = 1;
static constexpr size_t dM00 = 2;
static constexpr size_t d0P0 = 3;
static constexpr size_t d0M0 = 4;
static constexpr size_t d00P = 5;
static constexpr size_t d00M = 6;
static constexpr size_t dPP0 = 7;
static constexpr size_t dMM0 = 8;
static constexpr size_t dPM0 = 9;
static constexpr size_t dMP0 = 10;
static constexpr size_t dP0P = 11;
static constexpr size_t dM0M = 12;
static constexpr size_t dP0M = 13;
static constexpr size_t dM0P = 14;
static constexpr size_t d0PP = 15;
static constexpr size_t d0MM = 16;
static constexpr size_t d0PM = 17;
static constexpr size_t d0MP = 18;
static constexpr size_t dPPP = 19;
static constexpr size_t dMPP = 20;
static constexpr size_t dPMP = 21;
static constexpr size_t dMMP = 22;
static constexpr size_t dPPM = 23;
static constexpr size_t dMPM = 24;
static constexpr size_t dPMM = 25;
static constexpr size_t dMMM = 26;

#define NAME_OF(direction) #direction
static const std::map<size_t, std::string> directionNames {
    { d000, NAME_OF(d000) }, { dP00, NAME_OF(dP00) }, { dM00, NAME_OF(dM00) }, { d0P0, NAME_OF(d0P0) },
    { d0M0, NAME_OF(d0M0) }, { d00P, NAME_OF(d00P) }, { d00M, NAME_OF(d00M) }, { dPP0, NAME_OF(dPP0) },
    { dMM0, NAME_OF(dMM0) }, { dPM0, NAME_OF(dPM0) }, { dMP0, NAME_OF(dMP0) }, { dP0P, NAME_OF(dP0P) },
    { dM0M, NAME_OF(dM0M) }, { dP0M, NAME_OF(dP0M) }, { dM0P, NAME_OF(dM0P) }, { d0PP, NAME_OF(d0PP) },
    { d0MM, NAME_OF(d0MM) }, { d0PM, NAME_OF(d0PM) }, { d0MP, NAME_OF(d0MP) }, { dPPP, NAME_OF(dPPP) },
    { dMPP, NAME_OF(dMPP) }, { dPMP, NAME_OF(dPMP) }, { dMMP, NAME_OF(dMMP) }, { dPPM, NAME_OF(dPPM) },
    { dMPM, NAME_OF(dMPM) }, { dPMM, NAME_OF(dPMM) }, { dMMM, NAME_OF(dMMM) },
};

static constexpr size_t i000 = d000;
static constexpr size_t iP00 = dM00;
static constexpr size_t iM00 = dP00;
static constexpr size_t i0P0 = d0M0;
static constexpr size_t i0M0 = d0P0;
static constexpr size_t i00P = d00M;
static constexpr size_t i00M = d00P;
static constexpr size_t iPP0 = dMM0;
static constexpr size_t iMM0 = dPP0;
static constexpr size_t iPM0 = dMP0;
static constexpr size_t iMP0 = dPM0;
static constexpr size_t iP0P = dM0M;
static constexpr size_t iM0M = dP0P;
static constexpr size_t iP0M = dM0P;
static constexpr size_t iM0P = dP0M;
static constexpr size_t i0PP = d0MM;
static constexpr size_t i0MM = d0PP;
static constexpr size_t i0PM = d0MP;
static constexpr size_t i0MP = d0PM;
static constexpr size_t iPPP = dMMM;
static constexpr size_t iMPP = dPMM;
static constexpr size_t iPMP = dMPM;
static constexpr size_t iMMP = dPPM;
static constexpr size_t iPPM = dMMP;
static constexpr size_t iMPM = dPMP;
static constexpr size_t iPMM = dMPP;
static constexpr size_t iMMM = dPPP;

static constexpr size_t eP00 = 0;
static constexpr size_t eM00 = 0;
static constexpr size_t e0P0 = 1;
static constexpr size_t e0M0 = 1;
static constexpr size_t e00P = 2;
static constexpr size_t e00M = 2;
static constexpr size_t ePP0 = 3;
static constexpr size_t eMM0 = 3;
static constexpr size_t ePM0 = 4;
static constexpr size_t eMP0 = 4;
static constexpr size_t eP0P = 5;
static constexpr size_t eM0M = 5;
static constexpr size_t eP0M = 6;
static constexpr size_t eM0P = 6;
static constexpr size_t e0PP = 7;
static constexpr size_t e0MM = 7;
static constexpr size_t e0PM = 8;
static constexpr size_t e0MP = 8;
static constexpr size_t ePPP = 9;
static constexpr size_t eMMM = 9;
static constexpr size_t eMPP = 10;
static constexpr size_t ePMM = 10;
static constexpr size_t ePMP = 11;
static constexpr size_t eMPM = 11;
static constexpr size_t eMMP = 12;
static constexpr size_t ePPM = 12;

static constexpr unsigned long int et000 = 1;
static constexpr unsigned long int etP00 = 2;
static constexpr unsigned long int etM00 = 4;
static constexpr unsigned long int et0P0 = 8;
static constexpr unsigned long int et0M0 = 16;
static constexpr unsigned long int et00P = 32;
static constexpr unsigned long int et00M = 64;
static constexpr unsigned long int etPP0 = 128;
static constexpr unsigned long int etMM0 = 256;
static constexpr unsigned long int etPM0 = 512;
static constexpr unsigned long int etMP0 = 1024;
static constexpr unsigned long int etP0P = 2048;
static constexpr unsigned long int etM0M = 4096;
static constexpr unsigned long int etP0M = 8192;
static constexpr unsigned long int etM0P = 16384;
static constexpr unsigned long int et0PP = 32768;
static constexpr unsigned long int et0MM = 65536;
static constexpr unsigned long int et0PM = 131072;
static constexpr unsigned long int et0MP = 262144;
static constexpr unsigned long int etPPP = 524288;
static constexpr unsigned long int etMPP = 1048576;
static constexpr unsigned long int etPMP = 2097152;
static constexpr unsigned long int etMMP = 4194304;
static constexpr unsigned long int etPPM = 8388608;
static constexpr unsigned long int etMPM = 16777216;
static constexpr unsigned long int etPMM = 33554432;
static constexpr unsigned long int etMMM = 67108864;

template<size_t direction>
constexpr size_t inverseDir()
{
    switch(direction)
    {
        case d000: return i000;
        case dP00: return iP00;
        case dM00: return iM00;
        case d0P0: return i0P0;
        case d0M0: return i0M0;
        case d00P: return i00P;
        case d00M: return i00M;
        case dPP0: return iPP0;
        case dMM0: return iMM0;
        case dPM0: return iPM0;
        case dMP0: return iMP0;
        case dP0P: return iP0P;
        case dM0M: return iM0M;
        case dP0M: return iP0M;
        case dM0P: return iM0P;
        case d0PP: return i0PP;
        case d0MM: return i0MM;
        case d0PM: return i0PM;
        case d0MP: return i0MP;
        case dPPP: return iPPP;
        case dMMM: return iMMM;
        case dPPM: return iPPM;
        case dMMP: return iMMP;
        case dPMP: return iPMP;
        case dMPM: return iMPM;
        case dPMM: return iPMM;
        case dMPP: return iMPP;
    }
    return d000;
}


template<size_t direction>
constexpr real getVelocity(real velocityX, real velocityY, real velocityZ)
{
    switch(direction)
    {
        case d000: return basics::constant::c0o1;
        case dP00: return  velocityX;
        case dM00: return -velocityX;
        case d0P0: return velocityY;
        case d0M0: return -velocityY;
        case d00P: return velocityZ;
        case d00M: return -velocityZ;
        case dPP0: return velocityX + velocityY;
        case dMM0: return -velocityX - velocityY;
        case dPM0: return velocityX - velocityY;
        case dMP0: return -velocityX + velocityY;
        case dP0P: return velocityX + velocityZ;
        case dM0M: return -velocityX - velocityZ;
        case dP0M: return velocityX - velocityZ;
        case dM0P: return -velocityX + velocityZ;
        case d0PP: return velocityY + velocityZ;
        case d0MM: return -velocityY - velocityZ;
        case d0PM: return velocityY - velocityZ;
        case d0MP: return -velocityY + velocityZ;
        case dPPP: return velocityX + velocityY + velocityZ;
        case dMMM: return -velocityX - velocityY - velocityZ;
        case dPPM: return velocityX + velocityY - velocityZ;
        case dMMP: return -velocityX - velocityY + velocityZ;
        case dPMP: return velocityX - velocityY + velocityZ;
        case dMPM: return -velocityX + velocityY - velocityZ;
        case dPMM: return velocityX - velocityY - velocityZ;
        case dMPP: return -velocityX + velocityY + velocityZ;
    }
    return basics::constant::c0o1;
}

template<size_t direction>
constexpr real getWeight()
{
    switch(direction)
    {
        case d000: return basics::constant::c8o27;
        case dP00: return basics::constant::c2o27;
        case dM00: return basics::constant::c2o27;
        case d0P0: return basics::constant::c2o27;
        case d0M0: return basics::constant::c2o27;
        case d00P: return basics::constant::c2o27;
        case d00M: return basics::constant::c2o27;
        case dPP0: return basics::constant::c1o54;
        case dMM0: return basics::constant::c1o54;
        case dPM0: return basics::constant::c1o54;
        case dMP0: return basics::constant::c1o54;
        case dP0P: return basics::constant::c1o54;
        case dM0M: return basics::constant::c1o54;
        case dP0M: return basics::constant::c1o54;
        case dM0P: return basics::constant::c1o54;
        case d0PP: return basics::constant::c1o54;
        case d0MM: return basics::constant::c1o54;
        case d0PM: return basics::constant::c1o54;
        case d0MP: return basics::constant::c1o54;
        case dPPP: return basics::constant::c1o216;
        case dMMM: return basics::constant::c1o216;
        case dPPM: return basics::constant::c1o216;
        case dMMP: return basics::constant::c1o216;
        case dPMP: return basics::constant::c1o216;
        case dMPM: return basics::constant::c1o216;
        case dPMM: return basics::constant::c1o216;
        case dMPP: return basics::constant::c1o216;
    }
    return basics::constant::c8o27;
}
} // namespace vf::lbm::dir
#endif

//! \}
