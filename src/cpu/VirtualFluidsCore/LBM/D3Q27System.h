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
//! \file D3Q27System.h
//! \ingroup LBM
//! \author Konstantin Kutscher, Sebastian Geller, Soeren Freudiger
//=======================================================================================

#ifndef D3Q27SYSTEM_H
#define D3Q27SYSTEM_H

#include <cmath>
#include <iostream>
#include <string>

#include "LBMSystem.h"
#include "UbException.h"
#include "UbMath.h"


//! \brief namespace for global system-functions
namespace D3Q27System
{
//////////////////////////////////////////////////////////////////////////
// DIRECTION STUFF
static const int FSTARTDIR = 1;
static const int FENDDIR   = 26; // D3Q27

static const int STARTF = 0;
static const int ENDF   = 26; // D3Q27

static const int STARTDIR = 1; //0
static const int ENDDIR   = 27;//26 // all geometric directions

extern const int DX1[ENDDIR + 1];
extern const int DX2[ENDDIR + 1];
extern const int DX3[ENDDIR + 1];
extern const double WEIGTH[ENDDIR + 1];

extern const double cNorm[3][ENDDIR];

static const int MINLEVEL = 0;
static const int MAXLEVEL = 25;

extern const int EX1[ENDDIR + 1];
extern const int EX2[ENDDIR + 1];
extern const int EX3[ENDDIR + 1];

//static const int E    = 0;
//static const int W    = 1;
//static const int N    = 2;
//static const int S    = 3;
//static const int T    = 4;
//static const int B    = 5;
//static const int NE   = 6;
//static const int SW   = 7;
//static const int SE   = 8;
//static const int NW   = 9;
//static const int TE   = 10;
//static const int BW   = 11;
//static const int BE   = 12;
//static const int TW   = 13;
//static const int TN   = 14;
//static const int BS   = 15;
//static const int BN   = 16;
//static const int TS   = 17;
//static const int TNE  = 18;
//static const int TNW  = 19;
//static const int TSE  = 20;
//static const int TSW  = 21;
//static const int BNE  = 22;
//static const int BNW  = 23;
//static const int BSE  = 24;
//static const int BSW  = 25;
//static const int REST = 26;

//static constexpr int REST = 0;
//static constexpr int E = 1;
//static constexpr int W = 2;
//static constexpr int N = 3;
//static constexpr int S = 4;
//static constexpr int T = 5;
//static constexpr int B = 6;
//static constexpr int NE = 7;
//static constexpr int SW = 8;
//static constexpr int SE = 9;
//static constexpr int NW = 10;
//static constexpr int TE = 11;
//static constexpr int BW = 12;
//static constexpr int BE = 13;
//static constexpr int TW = 14;
//static constexpr int TN = 15;
//static constexpr int BS = 16;
//static constexpr int BN = 17;
//static constexpr int TS = 18;
//static constexpr int TNE = 19;
//static constexpr int TNW = 20;
//static constexpr int TSE = 21;
//static constexpr int TSW = 22;
//static constexpr int BNE = 23;
//static constexpr int BNW = 24;
//static constexpr int BSE = 25;
//static constexpr int BSW = 26;

static constexpr int DIR_000 = 0;
static constexpr int E = 1;
static constexpr int W = 2;
static constexpr int N = 3;
static constexpr int S = 4;
static constexpr int T = 5;
static constexpr int B = 6;
static constexpr int NE = 7;
static constexpr int SW = 8;
static constexpr int SE = 9;
static constexpr int NW = 10;
static constexpr int TE = 11;
static constexpr int BW = 12;
static constexpr int BE = 13;
static constexpr int TW = 14;
static constexpr int TN = 15;
static constexpr int BS = 16;
static constexpr int BN = 17;
static constexpr int TS = 18;
static constexpr int TNE = 19;
static constexpr int TNW = 20;
static constexpr int TSE = 21;
static constexpr int TSW = 22;
static constexpr int BNE = 23;
static constexpr int BNW = 24;
static constexpr int BSE = 25;
static constexpr int BSW = 26;

static const int INV_E   = W;
static const int INV_W   = E;
static const int INV_N   = S;
static const int INV_S   = N;
static const int INV_T   = B;
static const int INV_B   = T;
static const int INV_NE  = SW;
static const int INV_SW  = NE;
static const int INV_SE  = NW;
static const int INV_NW  = SE;
static const int INV_TE  = BW;
static const int INV_BW  = TE;
static const int INV_BE  = TW;
static const int INV_TW  = BE;
static const int INV_TN  = BS;
static const int INV_BS  = TN;
static const int INV_BN  = TS;
static const int INV_TS  = BN;
static const int INV_TNE = BSW;
static const int INV_TNW = BSE;
static const int INV_TSE = BNW;
static const int INV_TSW = BNE;
static const int INV_BNE = TSW;
static const int INV_BNW = TSE;
static const int INV_BSE = TNW;
static const int INV_BSW = TNE;

extern const int INVDIR[ENDDIR + 1];

static const int ET_E   = 0;
static const int ET_W   = 0;
static const int ET_N   = 1;
static const int ET_S   = 1;
static const int ET_T   = 2;
static const int ET_B   = 2;
static const int ET_NE  = 3;
static const int ET_SW  = 3;
static const int ET_SE  = 4;
static const int ET_NW  = 4;
static const int ET_TE  = 5;
static const int ET_BW  = 5;
static const int ET_BE  = 6;
static const int ET_TW  = 6;
static const int ET_TN  = 7;
static const int ET_BS  = 7;
static const int ET_BN  = 8;
static const int ET_TS  = 8;
static const int ET_TNE = 9;
static const int ET_BSW = 9;
static const int ET_TNW = 10;
static const int ET_BSE = 10;
static const int ET_TSE = 11;
static const int ET_BNW = 11;
static const int ET_TSW = 12;
static const int ET_BNE = 12;

//////////////////////////////////////////////////////////////////////////
inline std::string getDirectionString(int direction)
{
    switch (direction) {
        case E:
            return "E";
        case W:
            return "W";
        case N:
            return "N";
        case S:
            return "S";
        case T:
            return "T";
        case B:
            return "B";
        case NE:
            return "NE";
        case NW:
            return "NW";
        case SE:
            return "SE";
        case SW:
            return "SW";
        case TE:
            return "TE";
        case TW:
            return "TW";
        case BE:
            return "BE";
        case BW:
            return "BW";
        case TN:
            return "TN";
        case TS:
            return "TS";
        case BN:
            return "BN";
        case BS:
            return "BS";
        case TNE:
            return "TNE";
        case TNW:
            return "TNW";
        case TSE:
            return "TSE";
        case TSW:
            return "TSW";
        case BNE:
            return "BNE";
        case BNW:
            return "BNW";
        case BSE:
            return "BSE";
        case BSW:
            return "BSW";
        default:
            return "Cell3DSystem::getDrectionString(...) - unknown dir";
    }
}
//////////////////////////////////////////////////////////////////////////
static inline void setNeighborCoordinatesForDirection(int &x1, int &x2, int &x3, const int &direction)
{
    switch (direction) {
        case D3Q27System::E:
            x1++;
            break;
        case D3Q27System::N:
            x2++;
            break;
        case D3Q27System::T:
            x3++;
            break;
        case D3Q27System::W:
            x1--;
            break;
        case D3Q27System::S:
            x2--;
            break;
        case D3Q27System::B:
            x3--;
            break;
        case D3Q27System::NE:
            x1++;
            x2++;
            break;
        case D3Q27System::NW:
            x1--;
            x2++;
            break;
        case D3Q27System::SW:
            x1--;
            x2--;
            break;
        case D3Q27System::SE:
            x1++;
            x2--;
            break;
        case D3Q27System::TE:
            x1++;
            x3++;
            break;
        case D3Q27System::BW:
            x1--;
            x3--;
            break;
        case D3Q27System::BE:
            x1++;
            x3--;
            break;
        case D3Q27System::TW:
            x1--;
            x3++;
            break;
        case D3Q27System::TN:
            x2++;
            x3++;
            break;
        case D3Q27System::BS:
            x2--;
            x3--;
            break;
        case D3Q27System::BN:
            x2++;
            x3--;
            break;
        case D3Q27System::TS:
            x2--;
            x3++;
            break;
        case D3Q27System::TNE:
            x1++;
            x2++;
            x3++;
            break;
        case D3Q27System::TNW:
            x1--;
            x2++;
            x3++;
            break;
        case D3Q27System::TSE:
            x1++;
            x2--;
            x3++;
            break;
        case D3Q27System::TSW:
            x1--;
            x2--;
            x3++;
            break;
        case D3Q27System::BNE:
            x1++;
            x2++;
            x3--;
            break;
        case D3Q27System::BNW:
            x1--;
            x2++;
            x3--;
            break;
        case D3Q27System::BSE:
            x1++;
            x2--;
            x3--;
            break;
        case D3Q27System::BSW:
            x1--;
            x2--;
            x3--;
            break;
        default:
            throw UbException(UB_EXARGS, "no direction ...");
    }
}

//////////////////////////////////////////////////////////////////////////
// MACROSCOPIC VALUES
/*=====================================================================*/
LBMReal getDensity(const LBMReal *const &f /*[27]*/);
/*=====================================================================*/
static LBMReal getPressure(const LBMReal *const &f /*[27]*/) { return REAL_CAST(UbMath::c1o3) * getDensity(f); }
/*=====================================================================*/
LBMReal getIncompVelocityX1(const LBMReal *const &f /*[27]*/);
/*=====================================================================*/
LBMReal getIncompVelocityX2(const LBMReal *const &f /*[27]*/);
/*=====================================================================*/
LBMReal getIncompVelocityX3(const LBMReal *const &f /*[27]*/);


/*=====================================================================*/
static void calcDensity(const LBMReal *const &f /*[27]*/, LBMReal &rho)
{
    rho = ((f[TNE] + f[BSW]) + (f[TSE] + f[BNW])) + ((f[BSE] + f[TNW]) + (f[TSW] + f[BNE])) +
          (((f[NE] + f[SW]) + (f[SE] + f[NW])) + ((f[TE] + f[BW]) + (f[BE] + f[TW])) +
           ((f[BN] + f[TS]) + (f[TN] + f[BS]))) +
          ((f[E] + f[W]) + (f[N] + f[S]) + (f[T] + f[B])) + f[DIR_000];
}
/*=====================================================================*/
static void calcIncompVelocityX1(const LBMReal *const &f /*[27]*/, LBMReal &vx1)
{
    vx1 = ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[BSE] - f[TNW]) + (f[BNE] - f[TSW]))) +
           (((f[BE] - f[TW]) + (f[TE] - f[BW])) + ((f[SE] - f[NW]) + (f[NE] - f[SW]))) + (f[E] - f[W]));
}
/*=====================================================================*/
static void calcIncompVelocityX2(const LBMReal *const &f /*[27]*/, LBMReal &vx2)
{
    vx2 = ((((f[TNE] - f[BSW]) + (f[BNW] - f[TSE])) + ((f[TNW] - f[BSE]) + (f[BNE] - f[TSW]))) +
           (((f[BN] - f[TS]) + (f[TN] - f[BS])) + ((f[NW] - f[SE]) + (f[NE] - f[SW]))) + (f[N] - f[S]));
}
/*=====================================================================*/
static void calcIncompVelocityX3(const LBMReal *const &f /*[27]*/, LBMReal &vx3)
{
    vx3 = ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[TNW] - f[BSE]) + (f[TSW] - f[BNE]))) +
           (((f[TS] - f[BN]) + (f[TN] - f[BS])) + ((f[TW] - f[BE]) + (f[TE] - f[BW]))) + (f[T] - f[B]));
}
/*=====================================================================*/
static LBMReal getCompVelocityX1(const LBMReal *const &f /*[27]*/)
{
    return ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[BSE] - f[TNW]) + (f[BNE] - f[TSW]))) +
            (((f[BE] - f[TW]) + (f[TE] - f[BW])) + ((f[SE] - f[NW]) + (f[NE] - f[SW]))) + (f[E] - f[W])) /
           getDensity(f);
}
/*=====================================================================*/
static LBMReal getCompVelocityX2(const LBMReal *const &f /*[27]*/)
{
    return ((((f[TNE] - f[BSW]) + (f[BNW] - f[TSE])) + ((f[TNW] - f[BSE]) + (f[BNE] - f[TSW]))) +
            (((f[BN] - f[TS]) + (f[TN] - f[BS])) + ((f[NW] - f[SE]) + (f[NE] - f[SW]))) + (f[N] - f[S])) /
           getDensity(f);
}
/*=====================================================================*/
static LBMReal getCompVelocityX3(const LBMReal *const &f /*[27]*/)
{
    return ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[TNW] - f[BSE]) + (f[TSW] - f[BNE]))) +
            (((f[TS] - f[BN]) + (f[TN] - f[BS])) + ((f[TW] - f[BE]) + (f[TE] - f[BW]))) + (f[T] - f[B])) /
           getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX1(const LBMReal *const &f /*[27]*/, LBMReal &vx1)
{
    vx1 = ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[BSE] - f[TNW]) + (f[BNE] - f[TSW]))) +
           (((f[BE] - f[TW]) + (f[TE] - f[BW])) + ((f[SE] - f[NW]) + (f[NE] - f[SW]))) + (f[E] - f[W])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX2(const LBMReal *const &f /*[27]*/, LBMReal &vx2)
{
    vx2 = ((((f[TNE] - f[BSW]) + (f[BNW] - f[TSE])) + ((f[TNW] - f[BSE]) + (f[BNE] - f[TSW]))) +
           (((f[BN] - f[TS]) + (f[TN] - f[BS])) + ((f[NW] - f[SE]) + (f[NE] - f[SW]))) + (f[N] - f[S])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX3(const LBMReal *const &f /*[27]*/, LBMReal &vx3)
{
    vx3 = ((((f[TNE] - f[BSW]) + (f[TSE] - f[BNW])) + ((f[TNW] - f[BSE]) + (f[TSW] - f[BNE]))) +
           (((f[TS] - f[BN]) + (f[TN] - f[BS])) + ((f[TW] - f[BE]) + (f[TE] - f[BW]))) + (f[T] - f[B])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcIncompMacroscopicValues(const LBMReal *const &f /*[27]*/, LBMReal &rho, LBMReal &vx1, LBMReal &vx2,
                                        LBMReal &vx3)
{
    D3Q27System::calcDensity(f, rho);
    D3Q27System::calcIncompVelocityX1(f, vx1);
    D3Q27System::calcIncompVelocityX2(f, vx2);
    D3Q27System::calcIncompVelocityX3(f, vx3);
}

/*=====================================================================*/
static void calcCompMacroscopicValues(const LBMReal *const &f /*[27]*/, LBMReal &drho, LBMReal &vx1, LBMReal &vx2,
                                      LBMReal &vx3)
{
    D3Q27System::calcDensity(f, drho);
    D3Q27System::calcIncompVelocityX1(f, vx1);
    D3Q27System::calcIncompVelocityX2(f, vx2);
    D3Q27System::calcIncompVelocityX3(f, vx3);
    LBMReal rho = drho + UbMath::one;
    vx1 /= rho;
    vx2 /= rho;
    vx3 /= rho;
}
//////////////////////////////////////////////////////////////////////////
static LBMReal getCompFeqForDirection(const int &direction, const LBMReal &drho, const LBMReal &vx1, const LBMReal &vx2,
                                      const LBMReal &vx3)
{
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);
    LBMReal rho   = drho + UbMath::one;
    switch (direction) {
        case DIR_000:
            return REAL_CAST(UbMath::c8o27 * (drho + rho * (-cu_sq)));
        case E:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq)));
        case W:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq)));
        case N:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq)));
        case S:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq)));
        case T:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq)));
        case B:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq)));
        case NE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq)));
        case SW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq)));
        case SE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq)));
        case NW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq)));
        case TE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq)));
        case BW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq)));
        case BE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq)));
        case TW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq)));
        case TN:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq)));
        case BS:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq)));
        case BN:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq)));
        case TS:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq)));
        case TNE:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 + vx2 + vx3) +
                                            UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq)));
        case BSW:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (-vx1 - vx2 - vx3) +
                                            UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq)));
        case BNE:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 + vx2 - vx3) +
                                            UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq)));
        case TSW:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (-vx1 - vx2 + vx3) +
                                            UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq)));
        case TSE:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 - vx2 + vx3) +
                                            UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq)));
        case BNW:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (-vx1 + vx2 - vx3) +
                                            UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq)));
        case BSE:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 - vx2 - vx3) +
                                            UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq)));
        case TNW:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (-vx1 + vx2 + vx3) +
                                            UbMath::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq)));
        default:
            throw UbException(UB_EXARGS, "unknown dir");
    }
}
//////////////////////////////////////////////////////////////////////////
static void calcCompFeq(LBMReal *const &feq /*[27]*/, const LBMReal &drho, const LBMReal &vx1, const LBMReal &vx2,
                        const LBMReal &vx3)
{
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);
    LBMReal rho   = drho + UbMath::one;

    feq[DIR_000] = UbMath::c8o27 * (drho + rho * (-cu_sq));
    feq[E]    = UbMath::c2o27 * (drho + rho * (3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq));
    feq[W]    = UbMath::c2o27 * (drho + rho * (3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[N]    = UbMath::c2o27 * (drho + rho * (3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq));
    feq[S]    = UbMath::c2o27 * (drho + rho * (3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[T]    = UbMath::c2o27 * (drho + rho * (3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq));
    feq[B]    = UbMath::c2o27 * (drho + rho * (3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[NE]   = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[SW]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[SE]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[NW]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[TE]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[BW]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[BE]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[TW]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[TN]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[BS]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[BN]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[TS]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[TNE] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 + vx2 + vx3) + UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[BSW] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 - vx2 - vx3) + UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[BNE] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 + vx2 - vx3) + UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[TSW] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 - vx2 + vx3) + UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[TSE] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 - vx2 + vx3) + UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[BNW] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 + vx2 - vx3) + UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[BSE] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 - vx2 - vx3) + UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[TNW] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 + vx2 + vx3) + UbMath::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static LBMReal getIncompFeqForDirection(const int &direction, const LBMReal &drho, const LBMReal &vx1,
                                        const LBMReal &vx2, const LBMReal &vx3)
{
    LBMReal cu_sq = 1.5f * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    switch (direction) {
        case DIR_000:
            return REAL_CAST(UbMath::c8o27 * (drho - cu_sq));
        case E:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq));
        case W:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq));
        case N:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq));
        case S:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq));
        case T:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq));
        case B:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq));
        case NE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
        case SW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
        case SE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
        case NW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
        case TE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
        case BW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
        case BE:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
        case TW:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
        case TN:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
        case BS:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
        case BN:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
        case TS:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
        case TNE:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 + vx2 + vx3) +
                                               UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
        case BSW:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (-vx1 - vx2 - vx3) +
                                               UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
        case BNE:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 + vx2 - vx3) +
                                               UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
        case TSW:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (-vx1 - vx2 + vx3) +
                                               UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
        case TSE:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 - vx2 + vx3) +
                                               UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
        case BNW:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (-vx1 + vx2 - vx3) +
                                               UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
        case BSE:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 - vx2 - vx3) +
                                               UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
        case TNW:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (-vx1 + vx2 + vx3) +
                                               UbMath::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
        default:
            throw UbException(UB_EXARGS, "unknown dir");
    }
}
//////////////////////////////////////////////////////////////////////////
static void calcIncompFeq(LBMReal *const &feq /*[27]*/, const LBMReal &drho, const LBMReal &vx1, const LBMReal &vx2,
                          const LBMReal &vx3)
{
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[DIR_000] = UbMath::c8o27 * (drho - cu_sq);
    feq[E]    = UbMath::c2o27 * (drho + 3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq);
    feq[W]    = UbMath::c2o27 * (drho + 3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq);
    feq[N]    = UbMath::c2o27 * (drho + 3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq);
    feq[S]    = UbMath::c2o27 * (drho + 3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq);
    feq[T]    = UbMath::c2o27 * (drho + 3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq);
    feq[B]    = UbMath::c2o27 * (drho + 3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq);
    feq[NE]   = UbMath::c1o54 * (drho + 3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
    feq[SW]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
    feq[SE]   = UbMath::c1o54 * (drho + 3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
    feq[NW]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
    feq[TE]   = UbMath::c1o54 * (drho + 3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
    feq[BW]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
    feq[BE]   = UbMath::c1o54 * (drho + 3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
    feq[TW]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
    feq[TN]   = UbMath::c1o54 * (drho + 3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
    feq[BS]   = UbMath::c1o54 * (drho + 3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
    feq[BN]   = UbMath::c1o54 * (drho + 3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
    feq[TS]   = UbMath::c1o54 * (drho + 3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
    feq[TNE]  = UbMath::c1o216 *
               (drho + 3.0 * (vx1 + vx2 + vx3) + UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
    feq[BSW] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 - vx2 - vx3) + UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
    feq[BNE] = UbMath::c1o216 *
               (drho + 3.0 * (vx1 + vx2 - vx3) + UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
    feq[TSW] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 - vx2 + vx3) + UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
    feq[TSE] = UbMath::c1o216 *
               (drho + 3.0 * (vx1 - vx2 + vx3) + UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
    feq[BNW] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 + vx2 - vx3) + UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
    feq[BSE] = UbMath::c1o216 *
               (drho + 3.0 * (vx1 - vx2 - vx3) + UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
    feq[TNW] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 + vx2 + vx3) + UbMath::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
}
//////////////////////////////////////////////////////////////////////////
static inline float getBoundaryVelocityForDirection(const int &direction, const float &bcVelocityX1,
                                                    const float &bcVelocityX2, const float &bcVelocityX3)
{
    switch (direction) {
        case E:
            return (float)(UbMath::c4o9 * (+bcVelocityX1));
        case W:
            return (float)(UbMath::c4o9 * (-bcVelocityX1));
        case N:
            return (float)(UbMath::c4o9 * (+bcVelocityX2));
        case S:
            return (float)(UbMath::c4o9 * (-bcVelocityX2));
        case T:
            return (float)(UbMath::c4o9 * (+bcVelocityX3));
        case B:
            return (float)(UbMath::c4o9 * (-bcVelocityX3));
        case NE:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 + bcVelocityX2));
        case SW:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 - bcVelocityX2));
        case SE:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 - bcVelocityX2));
        case NW:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 + bcVelocityX2));
        case TE:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 + bcVelocityX3));
        case BW:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 - bcVelocityX3));
        case BE:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 - bcVelocityX3));
        case TW:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 + bcVelocityX3));
        case TN:
            return (float)(UbMath::c1o9 * (+bcVelocityX2 + bcVelocityX3));
        case BS:
            return (float)(UbMath::c1o9 * (-bcVelocityX2 - bcVelocityX3));
        case BN:
            return (float)(UbMath::c1o9 * (+bcVelocityX2 - bcVelocityX3));
        case TS:
            return (float)(UbMath::c1o9 * (-bcVelocityX2 + bcVelocityX3));
        case TNE:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
        case BSW:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
        case BNE:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
        case TSW:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
        case TSE:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
        case BNW:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
        case BSE:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
        case TNW:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
        default:
            throw UbException(UB_EXARGS, "unknown direction");
    }
}
/*=====================================================================*/
static const int &getInvertDirection(const int &direction)
{
#ifdef _DEBUG
    if (direction < STARTDIR || direction > ENDDIR)
        throw UbException(UB_EXARGS, "unknown direction");
#endif
    return INVDIR[direction];
}
/*=====================================================================*/
static void getLBMDirections(std::vector<int> &dirs, bool onlyLBdirs = false)
{
    std::vector<int> D3Q27Dirs;
    if (onlyLBdirs) /*FSTARTDIR->FENDDIR*/
    {
        dirs.resize(FENDDIR + 1);
        for (int dir = FSTARTDIR; dir <= FENDDIR; ++dir)
            dirs[dir] = dir;
    } else /*STARTDIR->ENDDIR*/
    {
        dirs.resize(ENDDIR + 1);
        for (int dir = STARTDIR; dir <= ENDDIR; ++dir)
            dirs[dir] = dir;
    }
}
//////////////////////////////////////////////////////////////////////////
static std::vector<int> getEX(const int &exn)
{
    std::vector<int> ex;
    ex.resize(ENDDIR + 1);
    switch (exn) {
        case 1:
            for (int dir = STARTDIR; dir < ENDDIR; ++dir)
                ex[dir] = DX1[dir];
            break;
        case 2:
            for (int dir = STARTDIR; dir < ENDDIR; ++dir)
                ex[dir] = DX2[dir];
            break;
        case 3:
            for (int dir = STARTDIR; dir < ENDDIR; ++dir)
                ex[dir] = DX3[dir];
            break;
    }
    return ex;
}
//////////////////////////////////////////////////////////////////////////
static inline void calcDistanceToNeighbors(std::vector<double> &distNeigh, const double &deltaX1)
{
    // distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);

    distNeigh[E] = distNeigh[W] = distNeigh[N] = deltaX1;
    distNeigh[S] = distNeigh[T] = distNeigh[B] = deltaX1;
    distNeigh[NE] = distNeigh[NW] = distNeigh[SW] = distNeigh[SE] = UbMath::sqrt2 * deltaX1;
    distNeigh[TE] = distNeigh[TN] = distNeigh[TW] = distNeigh[TS] = UbMath::sqrt2 * deltaX1;
    distNeigh[BE] = distNeigh[BN] = distNeigh[BW] = distNeigh[BS] = UbMath::sqrt2 * deltaX1;
    distNeigh[TNE] = distNeigh[TNW] = distNeigh[TSE] = distNeigh[TSW] = UbMath::sqrt3 * deltaX1;
    distNeigh[BNE] = distNeigh[BNW] = distNeigh[BSE] = distNeigh[BSW] = UbMath::sqrt3 * deltaX1;
}
//////////////////////////////////////////////////////////////////////////
static inline void calcDistanceToNeighbors(std::vector<double> &distNeigh, const double &deltaX1, const double &deltaX2,
                                           const double &deltaX3)
{
    // distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);
    distNeigh[E] = distNeigh[W] = deltaX1;
    distNeigh[N] = distNeigh[S] = deltaX2;
    distNeigh[T] = distNeigh[B] = deltaX3;
    distNeigh[NE] = distNeigh[NW] = distNeigh[SW] = distNeigh[SE] = sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2);
    distNeigh[TE] = distNeigh[TN] = distNeigh[TW] = distNeigh[TS] = sqrt(deltaX1 * deltaX1 + deltaX3 * deltaX3);
    distNeigh[BE] = distNeigh[BN] = distNeigh[BW] = distNeigh[BS] = sqrt(deltaX2 * deltaX2 + deltaX3 * deltaX3);
    distNeigh[TNE] = distNeigh[TNW] = distNeigh[TSE] = distNeigh[TSW] =
        sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3);
    distNeigh[BNE] = distNeigh[BNW] = distNeigh[BSE] = distNeigh[BSW] =
        sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3);
}
//////////////////////////////////////////////////////////////////////////
static inline void initRayVectors(double *const &rayX1, double *const &rayX2, double *const &rayX3)
{
    int fdir;
    double c1oS2 = UbMath::one_over_sqrt2;
    double c1oS3 = UbMath::one_over_sqrt3;
    fdir         = E;
    rayX1[fdir]  = 1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = W;
    rayX1[fdir]  = -1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = N;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 1.0;
    rayX3[fdir]  = 0.0;
    fdir         = S;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -1.0;
    rayX3[fdir]  = 0.0;
    fdir         = T;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 1.0;
    fdir         = B;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -1.0;
    fdir         = NE;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = SW;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = SE;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = NW;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = TE;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = BW;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = BE;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = TW;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = TN;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = c1oS2;
    fdir         = BS;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = BN;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = TS;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = c1oS2;
    fdir         = TNE;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = TNW;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = TSE;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = TSW;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = BNE;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = BNW;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = BSE;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = BSW;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = -c1oS3;
}
//////////////////////////////////////////////////////////////////////////
static inline LBMReal calcPress(const LBMReal *const f, LBMReal rho, LBMReal vx1, LBMReal vx2, LBMReal vx3)
{
    LBMReal op = 1.0;
    return ((f[E] + f[W] + f[N] + f[S] + f[T] + f[B] +
             2. * (f[NE] + f[SW] + f[SE] + f[NW] + f[TE] + f[BW] + f[BE] + f[TW] + f[TN] + f[BS] + f[BN] + f[TS]) +
             3. * (f[TNE] + f[TSW] + f[TSE] + f[TNW] + f[BNE] + f[BSW] + f[BSE] + f[BNW]) -
             (vx1 * vx1 + vx2 * vx2 + vx3 * vx3)) *
                (1 - 0.5 * op) +
            op * 0.5 * (rho)) *
           UbMath::c1o3;
}
//////////////////////////////////////////////////////////////////////////
static inline LBMReal getShearRate(const LBMReal *const f, LBMReal collFactorF)
{
    LBMReal mfcbb = f[E];
    LBMReal mfbcb = f[N];
    LBMReal mfbbc = f[T];
    LBMReal mfccb = f[NE];
    LBMReal mfacb = f[NW];
    LBMReal mfcbc = f[TE];
    LBMReal mfabc = f[TW];
    LBMReal mfbcc = f[TN];
    LBMReal mfbac = f[TS];
    LBMReal mfccc = f[TNE];
    LBMReal mfacc = f[TNW];
    LBMReal mfcac = f[TSE];
    LBMReal mfaac = f[TSW];

    LBMReal mfabb = f[W];
    LBMReal mfbab = f[S];
    LBMReal mfbba = f[B];
    LBMReal mfaab = f[SW];
    LBMReal mfcab = f[SE];
    LBMReal mfaba = f[BW];
    LBMReal mfcba = f[BE];
    LBMReal mfbaa = f[BS];
    LBMReal mfbca = f[BN];
    LBMReal mfaaa = f[BSW];
    LBMReal mfcaa = f[BSE];
    LBMReal mfaca = f[BNW];
    LBMReal mfcca = f[BNE];

    LBMReal mfbbb = f[DIR_000];

    LBMReal m0, m1, m2;

    LBMReal rho = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) + (mfaab + mfacb + mfcab + mfccb) +
                  (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) +
                  (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

    LBMReal vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                   (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) + (mfcbb - mfabb));
    LBMReal vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                   (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) + (mfbcb - mfbab));
    LBMReal vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                   (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) + (mfbbc - mfbba));

    LBMReal oMdrho;

    oMdrho = mfccc + mfaaa;
    m0     = mfaca + mfcac;
    m1     = mfacc + mfcaa;
    m2     = mfaac + mfcca;
    oMdrho += m0;
    m1 += m2;
    oMdrho += m1;
    m0 = mfbac + mfbca;
    m1 = mfbaa + mfbcc;
    m0 += m1;
    m1 = mfabc + mfcba;
    m2 = mfaba + mfcbc;
    m1 += m2;
    m0 += m1;
    m1 = mfacb + mfcab;
    m2 = mfaab + mfccb;
    m1 += m2;
    m0 += m1;
    oMdrho += m0;
    m0 = mfabb + mfcbb;
    m1 = mfbab + mfbcb;
    m2 = mfbba + mfbbc;
    m0 += m1 + m2;
    m0 += mfbbb; // hat gefehlt
    oMdrho = 1. - (oMdrho + m0);

    LBMReal vx2;
    LBMReal vy2;
    LBMReal vz2;
    vx2 = vvx * vvx;
    vy2 = vvy * vvy;
    vz2 = vvz * vvz;
    ////////////////////////////////////////////////////////////////////////////////////
    // Hin
    ////////////////////////////////////////////////////////////////////////////////////
    // mit 1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36  Konditionieren
    ////////////////////////////////////////////////////////////////////////////////////
    // Z - Dir
    m2    = mfaaa + mfaac;
    m1    = mfaac - mfaaa;
    m0    = m2 + mfaab;
    mfaaa = m0;
    m0 += UbMath::c1o36 * oMdrho;
    mfaab = m1 - m0 * vvz;
    mfaac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaba + mfabc;
    m1    = mfabc - mfaba;
    m0    = m2 + mfabb;
    mfaba = m0;
    m0 += UbMath::c1o9 * oMdrho;
    mfabb = m1 - m0 * vvz;
    mfabc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaca + mfacc;
    m1    = mfacc - mfaca;
    m0    = m2 + mfacb;
    mfaca = m0;
    m0 += UbMath::c1o36 * oMdrho;
    mfacb = m1 - m0 * vvz;
    mfacc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbaa + mfbac;
    m1    = mfbac - mfbaa;
    m0    = m2 + mfbab;
    mfbaa = m0;
    m0 += UbMath::c1o9 * oMdrho;
    mfbab = m1 - m0 * vvz;
    mfbac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbba + mfbbc;
    m1    = mfbbc - mfbba;
    m0    = m2 + mfbbb;
    mfbba = m0;
    m0 += UbMath::c4o9 * oMdrho;
    mfbbb = m1 - m0 * vvz;
    mfbbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbca + mfbcc;
    m1    = mfbcc - mfbca;
    m0    = m2 + mfbcb;
    mfbca = m0;
    m0 += UbMath::c1o9 * oMdrho;
    mfbcb = m1 - m0 * vvz;
    mfbcc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcaa + mfcac;
    m1    = mfcac - mfcaa;
    m0    = m2 + mfcab;
    mfcaa = m0;
    m0 += UbMath::c1o36 * oMdrho;
    mfcab = m1 - m0 * vvz;
    mfcac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcba + mfcbc;
    m1    = mfcbc - mfcba;
    m0    = m2 + mfcbb;
    mfcba = m0;
    m0 += UbMath::c1o9 * oMdrho;
    mfcbb = m1 - m0 * vvz;
    mfcbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcca + mfccc;
    m1    = mfccc - mfcca;
    m0    = m2 + mfccb;
    mfcca = m0;
    m0 += UbMath::c1o36 * oMdrho;
    mfccb = m1 - m0 * vvz;
    mfccc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    // mit  1/6, 0, 1/18, 2/3, 0, 2/9, 1/6, 0, 1/18 Konditionieren
    ////////////////////////////////////////////////////////////////////////////////////
    // Y - Dir
    m2    = mfaaa + mfaca;
    m1    = mfaca - mfaaa;
    m0    = m2 + mfaba;
    mfaaa = m0;
    m0 += UbMath::c1o6 * oMdrho;
    mfaba = m1 - m0 * vvy;
    mfaca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaab + mfacb;
    m1    = mfacb - mfaab;
    m0    = m2 + mfabb;
    mfaab = m0;
    mfabb = m1 - m0 * vvy;
    mfacb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaac + mfacc;
    m1    = mfacc - mfaac;
    m0    = m2 + mfabc;
    mfaac = m0;
    m0 += UbMath::c1o18 * oMdrho;
    mfabc = m1 - m0 * vvy;
    mfacc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbaa + mfbca;
    m1    = mfbca - mfbaa;
    m0    = m2 + mfbba;
    mfbaa = m0;
    m0 += UbMath::c2o3 * oMdrho;
    mfbba = m1 - m0 * vvy;
    mfbca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbab + mfbcb;
    m1    = mfbcb - mfbab;
    m0    = m2 + mfbbb;
    mfbab = m0;
    mfbbb = m1 - m0 * vvy;
    mfbcb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbac + mfbcc;
    m1    = mfbcc - mfbac;
    m0    = m2 + mfbbc;
    mfbac = m0;
    m0 += UbMath::c2o9 * oMdrho;
    mfbbc = m1 - m0 * vvy;
    mfbcc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcaa + mfcca;
    m1    = mfcca - mfcaa;
    m0    = m2 + mfcba;
    mfcaa = m0;
    m0 += UbMath::c1o6 * oMdrho;
    mfcba = m1 - m0 * vvy;
    mfcca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcab + mfccb;
    m1    = mfccb - mfcab;
    m0    = m2 + mfcbb;
    mfcab = m0;
    mfcbb = m1 - m0 * vvy;
    mfccb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcac + mfccc;
    m1    = mfccc - mfcac;
    m0    = m2 + mfcbc;
    mfcac = m0;
    m0 += UbMath::c1o18 * oMdrho;
    mfcbc = m1 - m0 * vvy;
    mfccc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    m2    = mfaaa + mfcaa;
    m1    = mfcaa - mfaaa;
    m0    = m2 + mfbaa;
    mfaaa = m0;
    m0 += 1. * oMdrho;
    mfbaa = m1 - m0 * vvx;
    mfcaa = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaba + mfcba;
    m1    = mfcba - mfaba;
    m0    = m2 + mfbba;
    mfaba = m0;
    mfbba = m1 - m0 * vvx;
    mfcba = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaca + mfcca;
    m1    = mfcca - mfaca;
    m0    = m2 + mfbca;
    mfaca = m0;
    m0 += UbMath::c1o3 * oMdrho;
    mfbca = m1 - m0 * vvx;
    mfcca = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaab + mfcab;
    m1    = mfcab - mfaab;
    m0    = m2 + mfbab;
    mfaab = m0;
    mfbab = m1 - m0 * vvx;
    mfcab = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfabb + mfcbb;
    m1    = mfcbb - mfabb;
    m0    = m2 + mfbbb;
    mfabb = m0;
    mfbbb = m1 - m0 * vvx;
    mfcbb = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfacb + mfccb;
    m1    = mfccb - mfacb;
    m0    = m2 + mfbcb;
    mfacb = m0;
    mfbcb = m1 - m0 * vvx;
    mfccb = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaac + mfcac;
    m1    = mfcac - mfaac;
    m0    = m2 + mfbac;
    mfaac = m0;
    m0 += UbMath::c1o3 * oMdrho;
    mfbac = m1 - m0 * vvx;
    mfcac = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfabc + mfcbc;
    m1    = mfcbc - mfabc;
    m0    = m2 + mfbbc;
    mfabc = m0;
    mfbbc = m1 - m0 * vvx;
    mfcbc = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfacc + mfccc;
    m1    = mfccc - mfacc;
    m0    = m2 + mfbcc;
    mfacc = m0;
    m0 += UbMath::c1o9 * oMdrho;
    mfbcc = m1 - m0 * vvx;
    mfccc = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    // Cumulants
    ////////////////////////////////////////////////////////////////////////////////////
    LBMReal OxxPyyPzz = 1.; // omega2 or bulk viscosity

    LBMReal mxxPyyPzz = mfcaa + mfaca + mfaac;
    LBMReal mxxMyy    = mfcaa - mfaca;
    LBMReal mxxMzz    = mfcaa - mfaac;

    LBMReal dxux = -UbMath::c1o2 * collFactorF * (mxxMyy + mxxMzz) + UbMath::c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
    LBMReal dyuy = dxux + collFactorF * UbMath::c3o2 * mxxMyy;
    LBMReal dzuz = dxux + collFactorF * UbMath::c3o2 * mxxMzz;

    LBMReal Dxy = -UbMath::three * collFactorF * mfbba;
    LBMReal Dxz = -UbMath::three * collFactorF * mfbab;
    LBMReal Dyz = -UbMath::three * collFactorF * mfabb;

    return sqrt(UbMath::c2 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz) + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz) /
           (rho + UbMath::one);
}
//Multiphase stuff
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseFeq(LBMReal *const &feq /*[27]*/, const LBMReal &rho, const LBMReal &p1, const LBMReal &vx1,
                              const LBMReal &vx2, const LBMReal &vx3)
{
    using namespace UbMath;
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[DIR_000] = c8o27 * (p1 + rho * c1o3 * (-cu_sq));
    feq[E]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq));
    feq[W]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[N]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq));
    feq[S]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[T]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq));
    feq[B]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[NE]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[SW]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[SE]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[NW]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[TE]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[BW]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[BE]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[TW]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[TN]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[BS]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[BN]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[TS]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[TNE] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[BSW] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[BNE] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[TSW] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[TSE] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[BNW] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[BSE] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[TNW] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseFeqVB(LBMReal *const &feq /*[27]*/, const LBMReal &p1, const LBMReal &vx1, const LBMReal &vx2,
                                const LBMReal &vx3)
{
    using namespace UbMath;
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[DIR_000] = p1 + c8o27 * (-cu_sq);
    feq[E]    = c2o27 * ((3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq));
    feq[W]    = c2o27 * ((3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[N]    = c2o27 * ((3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq));
    feq[S]    = c2o27 * ((3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[T]    = c2o27 * ((3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq));
    feq[B]    = c2o27 * ((3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[NE]   = c1o54 * ((3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[SW]   = c1o54 * ((3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[SE]   = c1o54 * ((3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[NW]   = c1o54 * ((3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[TE]   = c1o54 * ((3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[BW]   = c1o54 * ((3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[BE]   = c1o54 * ((3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[TW]   = c1o54 * ((3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[TN]   = c1o54 * ((3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[BS]   = c1o54 * ((3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[BN]   = c1o54 * ((3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[TS]   = c1o54 * ((3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[TNE]  = c1o216 * ((3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[BSW]  = c1o216 * ((3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[BNE]  = c1o216 * ((3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[TSW]  = c1o216 * ((3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[TSE]  = c1o216 * ((3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[BNW]  = c1o216 * ((3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[BSE]  = c1o216 * ((3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[TNW]  = c1o216 * ((3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseHeq(LBMReal *const &heq /*[27]*/, const LBMReal &phi, const LBMReal &vx1, const LBMReal &vx2,
                              const LBMReal &vx3)
{
    using namespace UbMath;
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    heq[DIR_000] = c8o27 * phi * (1.0 - cu_sq);
    heq[E]    = c2o27 * phi * (1.0 + 3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq);
    heq[W]    = c2o27 * phi * (1.0 + 3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq);
    heq[N]    = c2o27 * phi * (1.0 + 3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq);
    heq[S]    = c2o27 * phi * (1.0 + 3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq);
    heq[T]    = c2o27 * phi * (1.0 + 3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq);
    heq[B]    = c2o27 * phi * (1.0 + 3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq);
    heq[NE]   = c1o54 * phi * (1.0 + 3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
    heq[SW]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
    heq[SE]   = c1o54 * phi * (1.0 + 3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
    heq[NW]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
    heq[TE]   = c1o54 * phi * (1.0 + 3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
    heq[BW]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
    heq[BE]   = c1o54 * phi * (1.0 + 3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
    heq[TW]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
    heq[TN]   = c1o54 * phi * (1.0 + 3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
    heq[BS]   = c1o54 * phi * (1.0 + 3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
    heq[BN]   = c1o54 * phi * (1.0 + 3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
    heq[TS]   = c1o54 * phi * (1.0 + 3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
    heq[TNE]  = c1o216 * phi * (1.0 + 3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
    heq[BSW] = c1o216 * phi * (1.0 + 3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
    heq[BNE] = c1o216 * phi * (1.0 + 3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
    heq[TSW] = c1o216 * phi * (1.0 + 3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
    heq[TSE] = c1o216 * phi * (1.0 + 3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
    heq[BNW] = c1o216 * phi * (1.0 + 3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
    heq[BSE] = c1o216 * phi * (1.0 + 3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
    heq[TNW] = c1o216 * phi * (1.0 + 3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
}

} // namespace D3Q27System
#endif
