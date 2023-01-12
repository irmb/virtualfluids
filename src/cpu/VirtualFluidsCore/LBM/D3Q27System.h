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
#include <string>
#include <iostream>

#include "lbm/constants/D3Q27.h"
#include "LBMSystem.h"
#include "UbException.h"
#include "UbMath.h"

//using namespace vf::lbm::dir;

//! \brief namespace for global system-functions
namespace D3Q27System
{
//////////////////////////////////////////////////////////////////////////
// DIRECTION STUFF
static const int FSTARTDIR = 1;
static const int FENDDIR   = 26; // D3Q27

static const int STARTF = 0;
static const int ENDF   = 26; // D3Q27

//static const int STARTDIR = 1; //0
static const int ENDDIR   = 26;//26 // all geometric directions

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

//static constexpr int DIR_000 = 0;
//static constexpr int DIR_P00 = 1;
//static constexpr int DIR_M00 = 2;
//static constexpr int DIR_0P0 = 3;
//static constexpr int DIR_0M0 = 4;
//static constexpr int DIR_00P = 5;
//static constexpr int DIR_00M = 6;
//static constexpr int DIR_PP0 = 7;
//static constexpr int DIR_MM0 = 8;
//static constexpr int DIR_PM0 = 9;
//static constexpr int DIR_MP0 = 10;
//static constexpr int DIR_P0P = 11;
//static constexpr int DIR_M0M = 12;
//static constexpr int DIR_P0M = 13;
//static constexpr int DIR_M0P = 14;
//static constexpr int DIR_0PP = 15;
//static constexpr int DIR_0MM = 16;
//static constexpr int DIR_0PM = 17;
//static constexpr int DIR_0MP = 18;
//static constexpr int DIR_PPP = 19;
//static constexpr int DIR_MPP = 20;
//static constexpr int DIR_PMP = 21;
//static constexpr int DIR_MMP = 22;
//static constexpr int DIR_PPM = 23;
//static constexpr int DIR_MPM = 24;
//static constexpr int DIR_PMM = 25;
//static constexpr int DIR_MMM = 26;

//static constexpr int INV_P00 = DIR_M00;
//static constexpr int INV_M00 = DIR_P00;
//static constexpr int INV_0P0 = DIR_0M0;
//static constexpr int INV_0M0 = DIR_0P0;
//static constexpr int INV_00P = DIR_00M;
//static constexpr int INV_00M = DIR_00P;
//static constexpr int INV_PP0 = DIR_MM0;
//static constexpr int INV_MM0 = DIR_PP0;
//static constexpr int INV_PM0 = DIR_MP0;
//static constexpr int INV_MP0 = DIR_PM0;
//static constexpr int INV_P0P = DIR_M0M;
//static constexpr int INV_M0M = DIR_P0P;
//static constexpr int INV_P0M = DIR_M0P;
//static constexpr int INV_M0P = DIR_P0M;
//static constexpr int INV_0PP = DIR_0MM;
//static constexpr int INV_0MM = DIR_0PP;
//static constexpr int INV_0PM = DIR_0MP;
//static constexpr int INV_0MP = DIR_0PM;
//static constexpr int INV_PPP = DIR_MMM;
//static constexpr int INV_MPP = DIR_PMM;
//static constexpr int INV_PMP = DIR_MPM;
//static constexpr int INV_MMP = DIR_PPM;
//static constexpr int INV_PPM = DIR_MMP;
//static constexpr int INV_MPM = DIR_PMP;
//static constexpr int INV_PMM = DIR_MPP;
//static constexpr int INV_MMM = DIR_PPP;

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
    using namespace vf::lbm::dir;

    switch (direction) {
        case DIR_P00:
            return "E";
        case DIR_M00:
            return "W";
        case DIR_0P0:
            return "N";
        case DIR_0M0:
            return "S";
        case DIR_00P:
            return "T";
        case DIR_00M:
            return "B";
        case DIR_PP0:
            return "NE";
        case DIR_MP0:
            return "NW";
        case DIR_PM0:
            return "SE";
        case DIR_MM0:
            return "SW";
        case DIR_P0P:
            return "TE";
        case DIR_M0P:
            return "TW";
        case DIR_P0M:
            return "BE";
        case DIR_M0M:
            return "BW";
        case DIR_0PP:
            return "TN";
        case DIR_0MP:
            return "TS";
        case DIR_0PM:
            return "BN";
        case DIR_0MM:
            return "BS";
        case DIR_PPP:
            return "TNE";
        case DIR_MPP:
            return "TNW";
        case DIR_PMP:
            return "TSE";
        case DIR_MMP:
            return "TSW";
        case DIR_PPM:
            return "BNE";
        case DIR_MPM:
            return "BNW";
        case DIR_PMM:
            return "BSE";
        case DIR_MMM:
            return "BSW";
        default:
            return "Cell3DSystem::getDrectionString(...) - unknown dir";
    }
}
//////////////////////////////////////////////////////////////////////////
static inline void setNeighborCoordinatesForDirection(int &x1, int &x2, int &x3, const int &direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case DIR_P00:
            x1++;
            break;
        case DIR_0P0:
            x2++;
            break;
        case DIR_00P:
            x3++;
            break;
        case DIR_M00:
            x1--;
            break;
        case DIR_0M0:
            x2--;
            break;
        case DIR_00M:
            x3--;
            break;
        case DIR_PP0:
            x1++;
            x2++;
            break;
        case DIR_MP0:
            x1--;
            x2++;
            break;
        case DIR_MM0:
            x1--;
            x2--;
            break;
        case DIR_PM0:
            x1++;
            x2--;
            break;
        case DIR_P0P:
            x1++;
            x3++;
            break;
        case DIR_M0M:
            x1--;
            x3--;
            break;
        case DIR_P0M:
            x1++;
            x3--;
            break;
        case DIR_M0P:
            x1--;
            x3++;
            break;
        case DIR_0PP:
            x2++;
            x3++;
            break;
        case DIR_0MM:
            x2--;
            x3--;
            break;
        case DIR_0PM:
            x2++;
            x3--;
            break;
        case DIR_0MP:
            x2--;
            x3++;
            break;
        case DIR_PPP:
            x1++;
            x2++;
            x3++;
            break;
        case DIR_MPP:
            x1--;
            x2++;
            x3++;
            break;
        case DIR_PMP:
            x1++;
            x2--;
            x3++;
            break;
        case DIR_MMP:
            x1--;
            x2--;
            x3++;
            break;
        case DIR_PPM:
            x1++;
            x2++;
            x3--;
            break;
        case DIR_MPM:
            x1--;
            x2++;
            x3--;
            break;
        case DIR_PMM:
            x1++;
            x2--;
            x3--;
            break;
        case DIR_MMM:
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
    using namespace vf::lbm::dir;

    rho = ((f[DIR_PPP] + f[DIR_MMM]) + (f[DIR_PMP] + f[DIR_MPM])) + ((f[DIR_PMM] + f[DIR_MPP]) + (f[DIR_MMP] + f[DIR_PPM])) +
          (((f[DIR_PP0] + f[DIR_MM0]) + (f[DIR_PM0] + f[DIR_MP0])) + ((f[DIR_P0P] + f[DIR_M0M]) + (f[DIR_P0M] + f[DIR_M0P])) +
           ((f[DIR_0PM] + f[DIR_0MP]) + (f[DIR_0PP] + f[DIR_0MM]))) +
          ((f[DIR_P00] + f[DIR_M00]) + (f[DIR_0P0] + f[DIR_0M0]) + (f[DIR_00P] + f[DIR_00M])) + f[DIR_000];
}
/*=====================================================================*/
static void calcIncompVelocityX1(const LBMReal *const &f /*[27]*/, LBMReal &vx1)
{
    using namespace vf::lbm::dir;

    vx1 = ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_PMM] - f[DIR_MPP]) + (f[DIR_PPM] - f[DIR_MMP]))) +
           (((f[DIR_P0M] - f[DIR_M0P]) + (f[DIR_P0P] - f[DIR_M0M])) + ((f[DIR_PM0] - f[DIR_MP0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_P00] - f[DIR_M00]));
}
/*=====================================================================*/
static void calcIncompVelocityX2(const LBMReal *const &f /*[27]*/, LBMReal &vx2)
{
    using namespace vf::lbm::dir;

    vx2 = ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_MPM] - f[DIR_PMP])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_PPM] - f[DIR_MMP]))) +
           (((f[DIR_0PM] - f[DIR_0MP]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_MP0] - f[DIR_PM0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_0P0] - f[DIR_0M0]));
}
/*=====================================================================*/
static void calcIncompVelocityX3(const LBMReal *const &f /*[27]*/, LBMReal &vx3)
{
    using namespace vf::lbm::dir;

    vx3 = ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_MMP] - f[DIR_PPM]))) +
           (((f[DIR_0MP] - f[DIR_0PM]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_M0P] - f[DIR_P0M]) + (f[DIR_P0P] - f[DIR_M0M]))) + (f[DIR_00P] - f[DIR_00M]));
}
/*=====================================================================*/
static LBMReal getCompVelocityX1(const LBMReal *const &f /*[27]*/)
{
    using namespace vf::lbm::dir;

    return ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_PMM] - f[DIR_MPP]) + (f[DIR_PPM] - f[DIR_MMP]))) +
            (((f[DIR_P0M] - f[DIR_M0P]) + (f[DIR_P0P] - f[DIR_M0M])) + ((f[DIR_PM0] - f[DIR_MP0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_P00] - f[DIR_M00])) /
           getDensity(f);
}
/*=====================================================================*/
static LBMReal getCompVelocityX2(const LBMReal *const &f /*[27]*/)
{
    using namespace vf::lbm::dir;

    return ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_MPM] - f[DIR_PMP])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_PPM] - f[DIR_MMP]))) +
            (((f[DIR_0PM] - f[DIR_0MP]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_MP0] - f[DIR_PM0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_0P0] - f[DIR_0M0])) /
           getDensity(f);
}
/*=====================================================================*/
static LBMReal getCompVelocityX3(const LBMReal *const &f /*[27]*/)
{
    using namespace vf::lbm::dir;

    return ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_MMP] - f[DIR_PPM]))) +
            (((f[DIR_0MP] - f[DIR_0PM]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_M0P] - f[DIR_P0M]) + (f[DIR_P0P] - f[DIR_M0M]))) + (f[DIR_00P] - f[DIR_00M])) /
           getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX1(const LBMReal *const &f /*[27]*/, LBMReal &vx1)
{
    using namespace vf::lbm::dir;

    vx1 = ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_PMM] - f[DIR_MPP]) + (f[DIR_PPM] - f[DIR_MMP]))) +
           (((f[DIR_P0M] - f[DIR_M0P]) + (f[DIR_P0P] - f[DIR_M0M])) + ((f[DIR_PM0] - f[DIR_MP0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_P00] - f[DIR_M00])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX2(const LBMReal *const &f /*[27]*/, LBMReal &vx2)
{
    using namespace vf::lbm::dir;

    vx2 = ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_MPM] - f[DIR_PMP])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_PPM] - f[DIR_MMP]))) +
           (((f[DIR_0PM] - f[DIR_0MP]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_MP0] - f[DIR_PM0]) + (f[DIR_PP0] - f[DIR_MM0]))) + (f[DIR_0P0] - f[DIR_0M0])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX3(const LBMReal *const &f /*[27]*/, LBMReal &vx3)
{
    using namespace vf::lbm::dir;

    vx3 = ((((f[DIR_PPP] - f[DIR_MMM]) + (f[DIR_PMP] - f[DIR_MPM])) + ((f[DIR_MPP] - f[DIR_PMM]) + (f[DIR_MMP] - f[DIR_PPM]))) +
           (((f[DIR_0MP] - f[DIR_0PM]) + (f[DIR_0PP] - f[DIR_0MM])) + ((f[DIR_M0P] - f[DIR_P0M]) + (f[DIR_P0P] - f[DIR_M0M]))) + (f[DIR_00P] - f[DIR_00M])) /
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
    using namespace vf::lbm::dir;

    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);
    LBMReal rho   = drho + UbMath::one;
    switch (direction) {
        case DIR_000:
            return REAL_CAST(UbMath::c8o27 * (drho + rho * (-cu_sq)));
        case DIR_P00:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq)));
        case DIR_M00:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq)));
        case DIR_0P0:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq)));
        case DIR_0M0:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq)));
        case DIR_00P:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq)));
        case DIR_00M:
            return REAL_CAST(UbMath::c2o27 * (drho + rho * (3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq)));
        case DIR_PP0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq)));
        case DIR_MM0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq)));
        case DIR_PM0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq)));
        case DIR_MP0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq)));
        case DIR_P0P:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq)));
        case DIR_M0M:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq)));
        case DIR_P0M:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq)));
        case DIR_M0P:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq)));
        case DIR_0PP:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq)));
        case DIR_0MM:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq)));
        case DIR_0PM:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq)));
        case DIR_0MP:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + rho * (3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq)));
        case DIR_PPP:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 + vx2 + vx3) +
                                            UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq)));
        case DIR_MMM:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (-vx1 - vx2 - vx3) +
                                            UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq)));
        case DIR_PPM:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 + vx2 - vx3) +
                                            UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq)));
        case DIR_MMP:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (-vx1 - vx2 + vx3) +
                                            UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq)));
        case DIR_PMP:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 - vx2 + vx3) +
                                            UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq)));
        case DIR_MPM:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (-vx1 + vx2 - vx3) +
                                            UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq)));
        case DIR_PMM:
            return REAL_CAST(UbMath::c1o216 *
                             (drho + rho * (3.0 * (vx1 - vx2 - vx3) +
                                            UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq)));
        case DIR_MPP:
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
    using namespace vf::lbm::dir;

    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);
    LBMReal rho   = drho + UbMath::one;

    feq[DIR_000] = UbMath::c8o27 * (drho + rho * (-cu_sq));
    feq[DIR_P00]    = UbMath::c2o27 * (drho + rho * (3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq));
    feq[DIR_M00]    = UbMath::c2o27 * (drho + rho * (3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[DIR_0P0]    = UbMath::c2o27 * (drho + rho * (3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq));
    feq[DIR_0M0]    = UbMath::c2o27 * (drho + rho * (3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[DIR_00P]    = UbMath::c2o27 * (drho + rho * (3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq));
    feq[DIR_00M]    = UbMath::c2o27 * (drho + rho * (3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[DIR_PP0]   = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[DIR_MM0]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[DIR_PM0]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[DIR_MP0]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[DIR_P0P]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[DIR_M0M]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[DIR_P0M]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[DIR_M0P]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[DIR_0PP]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[DIR_0MM]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[DIR_0PM]  = UbMath::c1o54 * (drho + rho * (3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[DIR_0MP]  = UbMath::c1o54 * (drho + rho * (3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[DIR_PPP] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 + vx2 + vx3) + UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[DIR_MMM] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 - vx2 - vx3) + UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[DIR_PPM] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 + vx2 - vx3) + UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[DIR_MMP] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 - vx2 + vx3) + UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[DIR_PMP] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 - vx2 + vx3) + UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[DIR_MPM] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 + vx2 - vx3) + UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[DIR_PMM] = UbMath::c1o216 *
               (drho + rho * (3.0 * (vx1 - vx2 - vx3) + UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[DIR_MPP] =
        UbMath::c1o216 *
        (drho + rho * (3.0 * (-vx1 + vx2 + vx3) + UbMath::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static LBMReal getIncompFeqForDirection(const int &direction, const LBMReal &drho, const LBMReal &vx1,
                                        const LBMReal &vx2, const LBMReal &vx3)
{
    using namespace vf::lbm::dir;

    LBMReal cu_sq = 1.5f * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    switch (direction) {
        case DIR_000:
            return REAL_CAST(UbMath::c8o27 * (drho - cu_sq));
        case DIR_P00:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq));
        case DIR_M00:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq));
        case DIR_0P0:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq));
        case DIR_0M0:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq));
        case DIR_00P:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq));
        case DIR_00M:
            return REAL_CAST(UbMath::c2o27 * (drho + 3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq));
        case DIR_PP0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
        case DIR_MM0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
        case DIR_PM0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
        case DIR_MP0:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
        case DIR_P0P:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
        case DIR_M0M:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
        case DIR_P0M:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
        case DIR_M0P:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
        case DIR_0PP:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
        case DIR_0MM:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
        case DIR_0PM:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
        case DIR_0MP:
            return REAL_CAST(UbMath::c1o54 *
                             (drho + 3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
        case DIR_PPP:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 + vx2 + vx3) +
                                               UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
        case DIR_MMM:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (-vx1 - vx2 - vx3) +
                                               UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
        case DIR_PPM:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 + vx2 - vx3) +
                                               UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
        case DIR_MMP:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (-vx1 - vx2 + vx3) +
                                               UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
        case DIR_PMP:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 - vx2 + vx3) +
                                               UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
        case DIR_MPM:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (-vx1 + vx2 - vx3) +
                                               UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
        case DIR_PMM:
            return REAL_CAST(UbMath::c1o216 * (drho + 3.0 * (vx1 - vx2 - vx3) +
                                               UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
        case DIR_MPP:
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
    using namespace vf::lbm::dir;

    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[DIR_000] = UbMath::c8o27 * (drho - cu_sq);
    feq[DIR_P00]    = UbMath::c2o27 * (drho + 3.0 * (vx1) + UbMath::c9o2 * (vx1) * (vx1)-cu_sq);
    feq[DIR_M00]    = UbMath::c2o27 * (drho + 3.0 * (-vx1) + UbMath::c9o2 * (-vx1) * (-vx1) - cu_sq);
    feq[DIR_0P0]    = UbMath::c2o27 * (drho + 3.0 * (vx2) + UbMath::c9o2 * (vx2) * (vx2)-cu_sq);
    feq[DIR_0M0]    = UbMath::c2o27 * (drho + 3.0 * (-vx2) + UbMath::c9o2 * (-vx2) * (-vx2) - cu_sq);
    feq[DIR_00P]    = UbMath::c2o27 * (drho + 3.0 * (vx3) + UbMath::c9o2 * (vx3) * (vx3)-cu_sq);
    feq[DIR_00M]    = UbMath::c2o27 * (drho + 3.0 * (-vx3) + UbMath::c9o2 * (-vx3) * (-vx3) - cu_sq);
    feq[DIR_PP0]   = UbMath::c1o54 * (drho + 3.0 * (vx1 + vx2) + UbMath::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
    feq[DIR_MM0]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 - vx2) + UbMath::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
    feq[DIR_PM0]   = UbMath::c1o54 * (drho + 3.0 * (vx1 - vx2) + UbMath::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
    feq[DIR_MP0]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 + vx2) + UbMath::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
    feq[DIR_P0P]   = UbMath::c1o54 * (drho + 3.0 * (vx1 + vx3) + UbMath::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
    feq[DIR_M0M]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 - vx3) + UbMath::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
    feq[DIR_P0M]   = UbMath::c1o54 * (drho + 3.0 * (vx1 - vx3) + UbMath::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
    feq[DIR_M0P]   = UbMath::c1o54 * (drho + 3.0 * (-vx1 + vx3) + UbMath::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
    feq[DIR_0PP]   = UbMath::c1o54 * (drho + 3.0 * (vx2 + vx3) + UbMath::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
    feq[DIR_0MM]   = UbMath::c1o54 * (drho + 3.0 * (-vx2 - vx3) + UbMath::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
    feq[DIR_0PM]   = UbMath::c1o54 * (drho + 3.0 * (vx2 - vx3) + UbMath::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
    feq[DIR_0MP]   = UbMath::c1o54 * (drho + 3.0 * (-vx2 + vx3) + UbMath::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
    feq[DIR_PPP]  = UbMath::c1o216 *
               (drho + 3.0 * (vx1 + vx2 + vx3) + UbMath::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
    feq[DIR_MMM] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 - vx2 - vx3) + UbMath::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
    feq[DIR_PPM] = UbMath::c1o216 *
               (drho + 3.0 * (vx1 + vx2 - vx3) + UbMath::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
    feq[DIR_MMP] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 - vx2 + vx3) + UbMath::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
    feq[DIR_PMP] = UbMath::c1o216 *
               (drho + 3.0 * (vx1 - vx2 + vx3) + UbMath::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
    feq[DIR_MPM] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 + vx2 - vx3) + UbMath::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
    feq[DIR_PMM] = UbMath::c1o216 *
               (drho + 3.0 * (vx1 - vx2 - vx3) + UbMath::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
    feq[DIR_MPP] = UbMath::c1o216 *
               (drho + 3.0 * (-vx1 + vx2 + vx3) + UbMath::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
}
//////////////////////////////////////////////////////////////////////////
static inline float getBoundaryVelocityForDirection(const int &direction, const float &bcVelocityX1,
                                                    const float &bcVelocityX2, const float &bcVelocityX3)
{
    using namespace vf::lbm::dir;
 
    switch (direction) {
        case DIR_P00:
            return (float)(UbMath::c4o9 * (+bcVelocityX1));
        case DIR_M00:
            return (float)(UbMath::c4o9 * (-bcVelocityX1));
        case DIR_0P0:
            return (float)(UbMath::c4o9 * (+bcVelocityX2));
        case DIR_0M0:
            return (float)(UbMath::c4o9 * (-bcVelocityX2));
        case DIR_00P:
            return (float)(UbMath::c4o9 * (+bcVelocityX3));
        case DIR_00M:
            return (float)(UbMath::c4o9 * (-bcVelocityX3));
        case DIR_PP0:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 + bcVelocityX2));
        case DIR_MM0:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 - bcVelocityX2));
        case DIR_PM0:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 - bcVelocityX2));
        case DIR_MP0:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 + bcVelocityX2));
        case DIR_P0P:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 + bcVelocityX3));
        case DIR_M0M:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 - bcVelocityX3));
        case DIR_P0M:
            return (float)(UbMath::c1o9 * (+bcVelocityX1 - bcVelocityX3));
        case DIR_M0P:
            return (float)(UbMath::c1o9 * (-bcVelocityX1 + bcVelocityX3));
        case DIR_0PP:
            return (float)(UbMath::c1o9 * (+bcVelocityX2 + bcVelocityX3));
        case DIR_0MM:
            return (float)(UbMath::c1o9 * (-bcVelocityX2 - bcVelocityX3));
        case DIR_0PM:
            return (float)(UbMath::c1o9 * (+bcVelocityX2 - bcVelocityX3));
        case DIR_0MP:
            return (float)(UbMath::c1o9 * (-bcVelocityX2 + bcVelocityX3));
        case DIR_PPP:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
        case DIR_MMM:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
        case DIR_PPM:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
        case DIR_MMP:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
        case DIR_PMP:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
        case DIR_MPM:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
        case DIR_PMM:
            return (float)(UbMath::c1o36 * (+bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
        case DIR_MPP:
            return (float)(UbMath::c1o36 * (-bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
        default:
            throw UbException(UB_EXARGS, "unknown direction");
    }
}
/*=====================================================================*/
static const int &getInvertDirection(const int &direction)
{
#ifdef _DEBUG
 //   if (direction < STARTDIR || direction > ENDDIR)
     if (direction < FSTARTDIR || direction > FENDDIR)
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
        for (int dir = STARTF; dir <= ENDF; ++dir)
            dirs[dir] = dir;
    }
}
//////////////////////////////////////////////////////////////////////////
static std::vector<int> getDX(const int &exn)
{
    std::vector<int> ex;
    ex.resize(ENDDIR + 1);
    switch (exn) {
        case 1:
            for (int dir = FSTARTDIR; dir <= FENDDIR; ++dir)
                ex[dir] = DX1[dir];
            break;
        case 2:
            for (int dir = FSTARTDIR; dir <= FENDDIR; ++dir)
                ex[dir] = DX2[dir];
            break;
        case 3:
            for (int dir = FSTARTDIR; dir <= FENDDIR; ++dir)
                ex[dir] = DX3[dir];
            break;
    }
    return ex;
}
//////////////////////////////////////////////////////////////////////////
static inline void calcDistanceToNeighbors(std::vector<double> &distNeigh, const double &deltaX1)
{
    using namespace vf::lbm::dir;

    // distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);

    distNeigh[DIR_P00] = distNeigh[DIR_M00] = distNeigh[DIR_0P0] = deltaX1;
    distNeigh[DIR_0M0] = distNeigh[DIR_00P] = distNeigh[DIR_00M] = deltaX1;
    distNeigh[DIR_PP0] = distNeigh[DIR_MP0] = distNeigh[DIR_MM0] = distNeigh[DIR_PM0] = UbMath::sqrt2 * deltaX1;
    distNeigh[DIR_P0P] = distNeigh[DIR_0PP] = distNeigh[DIR_M0P] = distNeigh[DIR_0MP] = UbMath::sqrt2 * deltaX1;
    distNeigh[DIR_P0M] = distNeigh[DIR_0PM] = distNeigh[DIR_M0M] = distNeigh[DIR_0MM] = UbMath::sqrt2 * deltaX1;
    distNeigh[DIR_PPP] = distNeigh[DIR_MPP] = distNeigh[DIR_PMP] = distNeigh[DIR_MMP] = UbMath::sqrt3 * deltaX1;
    distNeigh[DIR_PPM] = distNeigh[DIR_MPM] = distNeigh[DIR_PMM] = distNeigh[DIR_MMM] = UbMath::sqrt3 * deltaX1;
}
//////////////////////////////////////////////////////////////////////////
static inline void calcDistanceToNeighbors(std::vector<double> &distNeigh, const double &deltaX1, const double &deltaX2,
                                           const double &deltaX3)
{
    using namespace vf::lbm::dir;

    // distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);
    distNeigh[DIR_P00] = distNeigh[DIR_M00] = deltaX1;
    distNeigh[DIR_0P0] = distNeigh[DIR_0M0] = deltaX2;
    distNeigh[DIR_00P] = distNeigh[DIR_00M] = deltaX3;
    distNeigh[DIR_PP0] = distNeigh[DIR_MP0] = distNeigh[DIR_MM0] = distNeigh[DIR_PM0] = sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2);
    distNeigh[DIR_P0P] = distNeigh[DIR_0PP] = distNeigh[DIR_M0P] = distNeigh[DIR_0MP] = sqrt(deltaX1 * deltaX1 + deltaX3 * deltaX3);
    distNeigh[DIR_P0M] = distNeigh[DIR_0PM] = distNeigh[DIR_M0M] = distNeigh[DIR_0MM] = sqrt(deltaX2 * deltaX2 + deltaX3 * deltaX3);
    distNeigh[DIR_PPP] = distNeigh[DIR_MPP] = distNeigh[DIR_PMP] = distNeigh[DIR_MMP] =
        sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3);
    distNeigh[DIR_PPM] = distNeigh[DIR_MPM] = distNeigh[DIR_PMM] = distNeigh[DIR_MMM] =
        sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3);
}
//////////////////////////////////////////////////////////////////////////
static inline void initRayVectors(double *const &rayX1, double *const &rayX2, double *const &rayX3)
{
    using namespace vf::lbm::dir;

    int fdir;
    double c1oS2 = UbMath::one_over_sqrt2;
    double c1oS3 = UbMath::one_over_sqrt3;
    fdir         = DIR_P00;
    rayX1[fdir]  = 1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_M00;
    rayX1[fdir]  = -1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_0P0;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 1.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_0M0;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -1.0;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_00P;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 1.0;
    fdir         = DIR_00M;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -1.0;
    fdir         = DIR_PP0;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_MM0;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_PM0;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_MP0;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = DIR_P0P;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = DIR_M0M;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_P0M;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_M0P;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = DIR_0PP;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = c1oS2;
    fdir         = DIR_0MM;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_0PM;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = DIR_0MP;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = c1oS2;
    fdir         = DIR_PPP;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = DIR_MPP;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = DIR_PMP;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = DIR_MMP;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = DIR_PPM;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = DIR_MPM;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = DIR_PMM;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = DIR_MMM;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = -c1oS3;
}
//////////////////////////////////////////////////////////////////////////
static inline LBMReal calcPress(const LBMReal *const f, LBMReal rho, LBMReal vx1, LBMReal vx2, LBMReal vx3)
{
    using namespace vf::lbm::dir;

    LBMReal op = 1.0;
    return ((f[DIR_P00] + f[DIR_M00] + f[DIR_0P0] + f[DIR_0M0] + f[DIR_00P] + f[DIR_00M] +
             2. * (f[DIR_PP0] + f[DIR_MM0] + f[DIR_PM0] + f[DIR_MP0] + f[DIR_P0P] + f[DIR_M0M] + f[DIR_P0M] + f[DIR_M0P] + f[DIR_0PP] + f[DIR_0MM] + f[DIR_0PM] + f[DIR_0MP]) +
             3. * (f[DIR_PPP] + f[DIR_MMP] + f[DIR_PMP] + f[DIR_MPP] + f[DIR_PPM] + f[DIR_MMM] + f[DIR_PMM] + f[DIR_MPM]) -
             (vx1 * vx1 + vx2 * vx2 + vx3 * vx3)) *
                (1 - 0.5 * op) +
            op * 0.5 * (rho)) *
           UbMath::c1o3;
}
//////////////////////////////////////////////////////////////////////////
static inline LBMReal getShearRate(const LBMReal *const f, LBMReal collFactorF)
{
    using namespace vf::lbm::dir;

    LBMReal mfcbb = f[DIR_P00];
    LBMReal mfbcb = f[DIR_0P0];
    LBMReal mfbbc = f[DIR_00P];
    LBMReal mfccb = f[DIR_PP0];
    LBMReal mfacb = f[DIR_MP0];
    LBMReal mfcbc = f[DIR_P0P];
    LBMReal mfabc = f[DIR_M0P];
    LBMReal mfbcc = f[DIR_0PP];
    LBMReal mfbac = f[DIR_0MP];
    LBMReal mfccc = f[DIR_PPP];
    LBMReal mfacc = f[DIR_MPP];
    LBMReal mfcac = f[DIR_PMP];
    LBMReal mfaac = f[DIR_MMP];

    LBMReal mfabb = f[DIR_M00];
    LBMReal mfbab = f[DIR_0M0];
    LBMReal mfbba = f[DIR_00M];
    LBMReal mfaab = f[DIR_MM0];
    LBMReal mfcab = f[DIR_PM0];
    LBMReal mfaba = f[DIR_M0M];
    LBMReal mfcba = f[DIR_P0M];
    LBMReal mfbaa = f[DIR_0MM];
    LBMReal mfbca = f[DIR_0PM];
    LBMReal mfaaa = f[DIR_MMM];
    LBMReal mfcaa = f[DIR_PMM];
    LBMReal mfaca = f[DIR_MPM];
    LBMReal mfcca = f[DIR_PPM];

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
    using namespace vf::lbm::dir;

    using namespace UbMath;
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[DIR_000] = c8o27 * (p1 + rho * c1o3 * (-cu_sq));
    feq[DIR_P00]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq));
    feq[DIR_M00]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[DIR_0P0]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq));
    feq[DIR_0M0]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[DIR_00P]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq));
    feq[DIR_00M]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[DIR_PP0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[DIR_MM0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[DIR_PM0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[DIR_MP0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[DIR_P0P]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[DIR_M0M]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[DIR_P0M]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[DIR_M0P]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[DIR_0PP]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[DIR_0MM]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[DIR_0PM]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[DIR_0MP]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[DIR_PPP] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[DIR_MMM] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[DIR_PPM] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[DIR_MMP] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[DIR_PMP] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[DIR_MPM] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[DIR_PMM] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[DIR_MPP] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseFeqVB(LBMReal *const &feq /*[27]*/, const LBMReal &p1, const LBMReal &vx1, const LBMReal &vx2,
                                const LBMReal &vx3)
{
    using namespace vf::lbm::dir;

    using namespace UbMath;
    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[DIR_000] = p1 + c8o27 * (-cu_sq);
    feq[DIR_P00]    = c2o27 * ((3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq));
    feq[DIR_M00]    = c2o27 * ((3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[DIR_0P0]    = c2o27 * ((3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq));
    feq[DIR_0M0]    = c2o27 * ((3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[DIR_00P]    = c2o27 * ((3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq));
    feq[DIR_00M]    = c2o27 * ((3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[DIR_PP0]   = c1o54 * ((3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[DIR_MM0]   = c1o54 * ((3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[DIR_PM0]   = c1o54 * ((3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[DIR_MP0]   = c1o54 * ((3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[DIR_P0P]   = c1o54 * ((3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[DIR_M0M]   = c1o54 * ((3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[DIR_P0M]   = c1o54 * ((3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[DIR_M0P]   = c1o54 * ((3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[DIR_0PP]   = c1o54 * ((3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[DIR_0MM]   = c1o54 * ((3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[DIR_0PM]   = c1o54 * ((3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[DIR_0MP]   = c1o54 * ((3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[DIR_PPP]  = c1o216 * ((3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[DIR_MMM]  = c1o216 * ((3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[DIR_PPM]  = c1o216 * ((3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[DIR_MMP]  = c1o216 * ((3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[DIR_PMP]  = c1o216 * ((3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[DIR_MPM]  = c1o216 * ((3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[DIR_PMM]  = c1o216 * ((3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[DIR_MPP]  = c1o216 * ((3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseHeq(LBMReal *const &heq /*[27]*/, const LBMReal &phi, const LBMReal &vx1, const LBMReal &vx2,
                              const LBMReal &vx3)
{
    using namespace vf::lbm::dir;
    using namespace UbMath;

    LBMReal cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    heq[DIR_000] = c8o27 * phi * (1.0 - cu_sq);
    heq[DIR_P00]    = c2o27 * phi * (1.0 + 3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq);
    heq[DIR_M00]    = c2o27 * phi * (1.0 + 3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq);
    heq[DIR_0P0]    = c2o27 * phi * (1.0 + 3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq);
    heq[DIR_0M0]    = c2o27 * phi * (1.0 + 3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq);
    heq[DIR_00P]    = c2o27 * phi * (1.0 + 3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq);
    heq[DIR_00M]    = c2o27 * phi * (1.0 + 3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq);
    heq[DIR_PP0]   = c1o54 * phi * (1.0 + 3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
    heq[DIR_MM0]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
    heq[DIR_PM0]   = c1o54 * phi * (1.0 + 3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
    heq[DIR_MP0]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
    heq[DIR_P0P]   = c1o54 * phi * (1.0 + 3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
    heq[DIR_M0M]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
    heq[DIR_P0M]   = c1o54 * phi * (1.0 + 3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
    heq[DIR_M0P]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
    heq[DIR_0PP]   = c1o54 * phi * (1.0 + 3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
    heq[DIR_0MM]   = c1o54 * phi * (1.0 + 3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
    heq[DIR_0PM]   = c1o54 * phi * (1.0 + 3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
    heq[DIR_0MP]   = c1o54 * phi * (1.0 + 3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
    heq[DIR_PPP]  = c1o216 * phi * (1.0 + 3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
    heq[DIR_MMM] = c1o216 * phi * (1.0 + 3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
    heq[DIR_PPM] = c1o216 * phi * (1.0 + 3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
    heq[DIR_MMP] = c1o216 * phi * (1.0 + 3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
    heq[DIR_PMP] = c1o216 * phi * (1.0 + 3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
    heq[DIR_MPM] = c1o216 * phi * (1.0 + 3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
    heq[DIR_PMM] = c1o216 * phi * (1.0 + 3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
    heq[DIR_MPP] = c1o216 * phi * (1.0 + 3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
}

} // namespace D3Q27System
#endif
