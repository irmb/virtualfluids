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
#include <array>

#include "lbm/constants/D3Q27.h"
#include "LBMSystem.h"
#include "UbException.h"
#include "UbMath.h"
#include "basics/constants/NumericConstants.h"

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
extern const real WEIGTH[ENDDIR + 1];

extern const real cNorm[3][ENDDIR];

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

//static constexpr int d000 = 0;
//static constexpr int dP00 = 1;
//static constexpr int dM00 = 2;
//static constexpr int d0P0 = 3;
//static constexpr int d0M0 = 4;
//static constexpr int d00P = 5;
//static constexpr int d00M = 6;
//static constexpr int dPP0 = 7;
//static constexpr int dMM0 = 8;
//static constexpr int dPM0 = 9;
//static constexpr int dMP0 = 10;
//static constexpr int dP0P = 11;
//static constexpr int dM0M = 12;
//static constexpr int dP0M = 13;
//static constexpr int dM0P = 14;
//static constexpr int d0PP = 15;
//static constexpr int d0MM = 16;
//static constexpr int d0PM = 17;
//static constexpr int d0MP = 18;
//static constexpr int dPPP = 19;
//static constexpr int dMPP = 20;
//static constexpr int dPMP = 21;
//static constexpr int dMMP = 22;
//static constexpr int dPPM = 23;
//static constexpr int dMPM = 24;
//static constexpr int dPMM = 25;
//static constexpr int dMMM = 26;

//static constexpr int iP00 = dM00;
//static constexpr int INV_M00 = dP00;
//static constexpr int INV_0P0 = d0M0;
//static constexpr int INV_0M0 = d0P0;
//static constexpr int INV_00P = d00M;
//static constexpr int INV_00M = d00P;
//static constexpr int INV_PP0 = dMM0;
//static constexpr int INV_MM0 = dPP0;
//static constexpr int INV_PM0 = dMP0;
//static constexpr int INV_MP0 = dPM0;
//static constexpr int INV_P0P = dM0M;
//static constexpr int INV_M0M = dP0P;
//static constexpr int INV_P0M = dM0P;
//static constexpr int INV_M0P = dP0M;
//static constexpr int INV_0PP = d0MM;
//static constexpr int INV_0MM = d0PP;
//static constexpr int INV_0PM = d0MP;
//static constexpr int INV_0MP = d0PM;
//static constexpr int INV_PPP = dMMM;
//static constexpr int INV_MPP = dPMM;
//static constexpr int INV_PMP = dMPM;
//static constexpr int INV_MMP = dPPM;
//static constexpr int INV_PPM = dMMP;
//static constexpr int INV_MPM = dPMP;
//static constexpr int INV_PMM = dMPP;
//static constexpr int INV_MMM = dPPP;

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

static const int ET_P00 = 0;
static const int ET_M00 = 0;
static const int ET_0P0 = 1;
static const int ET_0M0 = 1;
static const int ET_00P = 2;
static const int ET_00M = 2;
static const int ET_PP0 = 3;
static const int ET_MM0 = 3;
static const int ET_PM0 = 4;
static const int ET_MP0 = 4;
static const int ET_P0P = 5;
static const int ET_M0M = 5;
static const int ET_P0M = 6;
static const int ET_M0P = 6;
static const int ET_0PP = 7;
static const int ET_0MM = 7;
static const int ET_0PM = 8;
static const int ET_0MP = 8;
static const int ET_PPP = 9;
static const int ET_MMM = 9;
static const int ET_MPP = 10;
static const int ET_PMM = 10;
static const int ET_PMP = 11;
static const int ET_MPM = 11;
static const int ET_MMP = 12;
static const int ET_PPM = 12;

//////////////////////////////////////////////////////////////////////////
inline std::string getDirectionString(int direction)
{
    using namespace vf::lbm::dir;

    switch (direction) {
        case dP00:
            return "E";
        case dM00:
            return "W";
        case d0P0:
            return "N";
        case d0M0:
            return "S";
        case d00P:
            return "T";
        case d00M:
            return "B";
        case dPP0:
            return "NE";
        case dMP0:
            return "NW";
        case dPM0:
            return "SE";
        case dMM0:
            return "SW";
        case dP0P:
            return "TE";
        case dM0P:
            return "TW";
        case dP0M:
            return "BE";
        case dM0M:
            return "BW";
        case d0PP:
            return "TN";
        case d0MP:
            return "TS";
        case d0PM:
            return "BN";
        case d0MM:
            return "BS";
        case dPPP:
            return "TNE";
        case dMPP:
            return "TNW";
        case dPMP:
            return "TSE";
        case dMMP:
            return "TSW";
        case dPPM:
            return "BNE";
        case dMPM:
            return "BNW";
        case dPMM:
            return "BSE";
        case dMMM:
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
        case dP00:
            x1++;
            break;
        case d0P0:
            x2++;
            break;
        case d00P:
            x3++;
            break;
        case dM00:
            x1--;
            break;
        case d0M0:
            x2--;
            break;
        case d00M:
            x3--;
            break;
        case dPP0:
            x1++;
            x2++;
            break;
        case dMP0:
            x1--;
            x2++;
            break;
        case dMM0:
            x1--;
            x2--;
            break;
        case dPM0:
            x1++;
            x2--;
            break;
        case dP0P:
            x1++;
            x3++;
            break;
        case dM0M:
            x1--;
            x3--;
            break;
        case dP0M:
            x1++;
            x3--;
            break;
        case dM0P:
            x1--;
            x3++;
            break;
        case d0PP:
            x2++;
            x3++;
            break;
        case d0MM:
            x2--;
            x3--;
            break;
        case d0PM:
            x2++;
            x3--;
            break;
        case d0MP:
            x2--;
            x3++;
            break;
        case dPPP:
            x1++;
            x2++;
            x3++;
            break;
        case dMPP:
            x1--;
            x2++;
            x3++;
            break;
        case dPMP:
            x1++;
            x2--;
            x3++;
            break;
        case dMMP:
            x1--;
            x2--;
            x3++;
            break;
        case dPPM:
            x1++;
            x2++;
            x3--;
            break;
        case dMPM:
            x1--;
            x2++;
            x3--;
            break;
        case dPMM:
            x1++;
            x2--;
            x3--;
            break;
        case dMMM:
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
real getDensity(const real *const &f /*[27]*/);
/*=====================================================================*/
static real getPressure(const real *const &f /*[27]*/) { return REAL_CAST(vf::basics::constant::c1o3) * getDensity(f); }
/*=====================================================================*/
real getIncompVelocityX1(const real *const &f /*[27]*/);
/*=====================================================================*/
real getIncompVelocityX2(const real *const &f /*[27]*/);
/*=====================================================================*/
real getIncompVelocityX3(const real *const &f /*[27]*/);


/*=====================================================================*/
static void calcDensity(const real *const &f /*[27]*/, real &rho)
{
    using namespace vf::lbm::dir;

    rho = ((f[dPPP] + f[dMMM]) + (f[dPMP] + f[dMPM])) + ((f[dPMM] + f[dMPP]) + (f[dMMP] + f[dPPM])) +
          (((f[dPP0] + f[dMM0]) + (f[dPM0] + f[dMP0])) + ((f[dP0P] + f[dM0M]) + (f[dP0M] + f[dM0P])) +
           ((f[d0PM] + f[d0MP]) + (f[d0PP] + f[d0MM]))) +
          ((f[dP00] + f[dM00]) + (f[d0P0] + f[d0M0]) + (f[d00P] + f[d00M])) + f[d000];
}
/*=====================================================================*/
static void calcIncompVelocityX1(const real *const &f /*[27]*/, real &vx1)
{
    using namespace vf::lbm::dir;

    vx1 = ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dPMM] - f[dMPP]) + (f[dPPM] - f[dMMP]))) +
           (((f[dP0M] - f[dM0P]) + (f[dP0P] - f[dM0M])) + ((f[dPM0] - f[dMP0]) + (f[dPP0] - f[dMM0]))) + (f[dP00] - f[dM00]));
}
/*=====================================================================*/
static void calcIncompVelocityX2(const real *const &f /*[27]*/, real &vx2)
{
    using namespace vf::lbm::dir;

    vx2 = ((((f[dPPP] - f[dMMM]) + (f[dMPM] - f[dPMP])) + ((f[dMPP] - f[dPMM]) + (f[dPPM] - f[dMMP]))) +
           (((f[d0PM] - f[d0MP]) + (f[d0PP] - f[d0MM])) + ((f[dMP0] - f[dPM0]) + (f[dPP0] - f[dMM0]))) + (f[d0P0] - f[d0M0]));
}
/*=====================================================================*/
static void calcIncompVelocityX3(const real *const &f /*[27]*/, real &vx3)
{
    using namespace vf::lbm::dir;

    vx3 = ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dMPP] - f[dPMM]) + (f[dMMP] - f[dPPM]))) +
           (((f[d0MP] - f[d0PM]) + (f[d0PP] - f[d0MM])) + ((f[dM0P] - f[dP0M]) + (f[dP0P] - f[dM0M]))) + (f[d00P] - f[d00M]));
}
/*=====================================================================*/
static real getCompVelocityX1(const real *const &f /*[27]*/)
{
    using namespace vf::lbm::dir;

    return ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dPMM] - f[dMPP]) + (f[dPPM] - f[dMMP]))) +
            (((f[dP0M] - f[dM0P]) + (f[dP0P] - f[dM0M])) + ((f[dPM0] - f[dMP0]) + (f[dPP0] - f[dMM0]))) + (f[dP00] - f[dM00])) /
           getDensity(f);
}
/*=====================================================================*/
static real getCompVelocityX2(const real *const &f /*[27]*/)
{
    using namespace vf::lbm::dir;

    return ((((f[dPPP] - f[dMMM]) + (f[dMPM] - f[dPMP])) + ((f[dMPP] - f[dPMM]) + (f[dPPM] - f[dMMP]))) +
            (((f[d0PM] - f[d0MP]) + (f[d0PP] - f[d0MM])) + ((f[dMP0] - f[dPM0]) + (f[dPP0] - f[dMM0]))) + (f[d0P0] - f[d0M0])) /
           getDensity(f);
}
/*=====================================================================*/
static real getCompVelocityX3(const real *const &f /*[27]*/)
{
    using namespace vf::lbm::dir;

    return ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dMPP] - f[dPMM]) + (f[dMMP] - f[dPPM]))) +
            (((f[d0MP] - f[d0PM]) + (f[d0PP] - f[d0MM])) + ((f[dM0P] - f[dP0M]) + (f[dP0P] - f[dM0M]))) + (f[d00P] - f[d00M])) /
           getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX1(const real *const &f /*[27]*/, real &vx1)
{
    using namespace vf::lbm::dir;

    vx1 = ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dPMM] - f[dMPP]) + (f[dPPM] - f[dMMP]))) +
           (((f[dP0M] - f[dM0P]) + (f[dP0P] - f[dM0M])) + ((f[dPM0] - f[dMP0]) + (f[dPP0] - f[dMM0]))) + (f[dP00] - f[dM00])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX2(const real *const &f /*[27]*/, real &vx2)
{
    using namespace vf::lbm::dir;

    vx2 = ((((f[dPPP] - f[dMMM]) + (f[dMPM] - f[dPMP])) + ((f[dMPP] - f[dPMM]) + (f[dPPM] - f[dMMP]))) +
           (((f[d0PM] - f[d0MP]) + (f[d0PP] - f[d0MM])) + ((f[dMP0] - f[dPM0]) + (f[dPP0] - f[dMM0]))) + (f[d0P0] - f[d0M0])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcCompVelocityX3(const real *const &f /*[27]*/, real &vx3)
{
    using namespace vf::lbm::dir;

    vx3 = ((((f[dPPP] - f[dMMM]) + (f[dPMP] - f[dMPM])) + ((f[dMPP] - f[dPMM]) + (f[dMMP] - f[dPPM]))) +
           (((f[d0MP] - f[d0PM]) + (f[d0PP] - f[d0MM])) + ((f[dM0P] - f[dP0M]) + (f[dP0P] - f[dM0M]))) + (f[d00P] - f[d00M])) /
          getDensity(f);
}
/*=====================================================================*/
static void calcIncompMacroscopicValues(const real *const &f /*[27]*/, real &rho, real &vx1, real &vx2,
                                        real &vx3)
{
    D3Q27System::calcDensity(f, rho);
    D3Q27System::calcIncompVelocityX1(f, vx1);
    D3Q27System::calcIncompVelocityX2(f, vx2);
    D3Q27System::calcIncompVelocityX3(f, vx3);
}

/*=====================================================================*/
static void calcCompMacroscopicValues(const real *const &f /*[27]*/, real &drho, real &vx1, real &vx2,
                                      real &vx3)
{
    D3Q27System::calcDensity(f, drho);
    D3Q27System::calcIncompVelocityX1(f, vx1);
    D3Q27System::calcIncompVelocityX2(f, vx2);
    D3Q27System::calcIncompVelocityX3(f, vx3);
    //real rho = drho + vf::basics::constant::one;
    real rho = drho + vf::basics::constant::c1o1;
    vx1 /= rho;
    vx2 /= rho;
    vx3 /= rho;
}
//////////////////////////////////////////////////////////////////////////
static real getCompFeqForDirection(const int &direction, const real &drho, const real &vx1, const real &vx2,
                                      const real &vx3)
{
    using namespace vf::lbm::dir;

    real cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);
    real rho   = drho + vf::basics::constant::c1o1;
    switch (direction) {
        case d000:
            return REAL_CAST(vf::basics::constant::c8o27 * (drho + rho * (-cu_sq)));
        case dP00:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + rho * (3.0 * (vx1) +vf::basics::constant::c9o2 * (vx1) * (vx1)-cu_sq)));
        case dM00:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + rho * (3.0 * (-vx1) + vf::basics::constant::c9o2 * (-vx1) * (-vx1) - cu_sq)));
        case d0P0:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + rho * (3.0 * (vx2) +vf::basics::constant::c9o2 * (vx2) * (vx2)-cu_sq)));
        case d0M0:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + rho * (3.0 * (-vx2) + vf::basics::constant::c9o2 * (-vx2) * (-vx2) - cu_sq)));
        case d00P:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + rho * (3.0 * (vx3) + vf::basics::constant::c9o2 * (vx3) * (vx3)-cu_sq)));
        case d00M:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + rho * (3.0 * (-vx3) + vf::basics::constant::c9o2 * (-vx3) * (-vx3) - cu_sq)));
        case dPP0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (vx1 + vx2) + vf::basics::constant::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq)));
        case dMM0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (-vx1 - vx2) + vf::basics::constant::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq)));
        case dPM0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (vx1 - vx2) + vf::basics::constant::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq)));
        case dMP0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (-vx1 + vx2) + vf::basics::constant::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq)));
        case dP0P:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (vx1 + vx3) + vf::basics::constant::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq)));
        case dM0M:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (-vx1 - vx3) + vf::basics::constant::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq)));
        case dP0M:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (vx1 - vx3) + vf::basics::constant::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq)));
        case dM0P:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (-vx1 + vx3) + vf::basics::constant::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq)));
        case d0PP:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (vx2 + vx3) + vf::basics::constant::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq)));
        case d0MM:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (-vx2 - vx3) + vf::basics::constant::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq)));
        case d0PM:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (vx2 - vx3) + vf::basics::constant::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq)));
        case d0MP:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + rho * (3.0 * (-vx2 + vx3) + vf::basics::constant::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq)));
        case dPPP:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (vx1 + vx2 + vx3) +
                                 vf::basics::constant::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq)));
        case dMMM:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (-vx1 - vx2 - vx3) +
                                 vf::basics::constant::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq)));
        case dPPM:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (vx1 + vx2 - vx3) +
                                 vf::basics::constant::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq)));
        case dMMP:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (-vx1 - vx2 + vx3) +
                                            vf::basics::constant::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq)));
        case dPMP:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (vx1 - vx2 + vx3) +
                                 vf::basics::constant::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq)));
        case dMPM:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (-vx1 + vx2 - vx3) +
                                 vf::basics::constant::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq)));
        case dPMM:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (vx1 - vx2 - vx3) +
                                 vf::basics::constant::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq)));
        case dMPP:
            return REAL_CAST(vf::basics::constant::c1o216 *
                             (drho + rho * (3.0 * (-vx1 + vx2 + vx3) +
                                 vf::basics::constant::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq)));
        default:
            throw UbException(UB_EXARGS, "unknown dir");
    }
}
//////////////////////////////////////////////////////////////////////////
static void calcCompFeq(real *const &feq /*[27]*/, const real &drho, const real &vx1, const real &vx2,
                        const real &vx3)
{
    using namespace vf::lbm::dir;

    real cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);
    real rho   = drho + vf::basics::constant::c1o1;

    feq[d000] = vf::basics::constant::c8o27 * (drho + rho * (-cu_sq));
    feq[dP00]    = vf::basics::constant::c2o27 * (drho + rho * (3.0 * (vx1) + vf::basics::constant::c9o2 * (vx1) * (vx1)-cu_sq));
    feq[dM00]    = vf::basics::constant::c2o27 * (drho + rho * (3.0 * (-vx1) + vf::basics::constant::c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[d0P0]    = vf::basics::constant::c2o27 * (drho + rho * (3.0 * (vx2) + vf::basics::constant::c9o2 * (vx2) * (vx2)-cu_sq));
    feq[d0M0]    = vf::basics::constant::c2o27 * (drho + rho * (3.0 * (-vx2) + vf::basics::constant::c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[d00P]    = vf::basics::constant::c2o27 * (drho + rho * (3.0 * (vx3) + vf::basics::constant::c9o2 * (vx3) * (vx3)-cu_sq));
    feq[d00M]    = vf::basics::constant::c2o27 * (drho + rho * (3.0 * (-vx3) + vf::basics::constant::c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[dPP0]   = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (vx1 + vx2) + vf::basics::constant::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[dMM0]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (-vx1 - vx2) + vf::basics::constant::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[dPM0]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (vx1 - vx2) + vf::basics::constant::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[dMP0]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (-vx1 + vx2) + vf::basics::constant::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[dP0P]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (vx1 + vx3) + vf::basics::constant::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[dM0M]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (-vx1 - vx3) + vf::basics::constant::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[dP0M]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (vx1 - vx3) + vf::basics::constant::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[dM0P]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (-vx1 + vx3) + vf::basics::constant::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[d0PP]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (vx2 + vx3) + vf::basics::constant::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[d0MM]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (-vx2 - vx3) + vf::basics::constant::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[d0PM]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (vx2 - vx3) + vf::basics::constant::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[d0MP]  = vf::basics::constant::c1o54 * (drho + rho * (3.0 * (-vx2 + vx3) + vf::basics::constant::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[dPPP] = vf::basics::constant::c1o216 *
               (drho + rho * (3.0 * (vx1 + vx2 + vx3) + vf::basics::constant::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[dMMM] =
        vf::basics::constant::c1o216 *
        (drho + rho * (3.0 * (-vx1 - vx2 - vx3) + vf::basics::constant::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[dPPM] = vf::basics::constant::c1o216 *
               (drho + rho * (3.0 * (vx1 + vx2 - vx3) + vf::basics::constant::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[dMMP] =
        vf::basics::constant::c1o216 *
        (drho + rho * (3.0 * (-vx1 - vx2 + vx3) + vf::basics::constant::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[dPMP] = vf::basics::constant::c1o216 *
               (drho + rho * (3.0 * (vx1 - vx2 + vx3) + vf::basics::constant::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[dMPM] =
        vf::basics::constant::c1o216 *
        (drho + rho * (3.0 * (-vx1 + vx2 - vx3) + vf::basics::constant::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[dPMM] = vf::basics::constant::c1o216 *
               (drho + rho * (3.0 * (vx1 - vx2 - vx3) + vf::basics::constant::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[dMPP] =
        vf::basics::constant::c1o216 *
        (drho + rho * (3.0 * (-vx1 + vx2 + vx3) + vf::basics::constant::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static real getIncompFeqForDirection(const int &direction, const real &drho, const real &vx1,
                                        const real &vx2, const real &vx3)
{
    using namespace vf::lbm::dir;

    real cu_sq = 1.5f * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    switch (direction) {
        case d000:
            return REAL_CAST(vf::basics::constant::c8o27 * (drho - cu_sq));
        case dP00:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + 3.0 * (vx1) + vf::basics::constant::c9o2 * (vx1) * (vx1)-cu_sq));
        case dM00:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + 3.0 * (-vx1) + vf::basics::constant::c9o2 * (-vx1) * (-vx1) - cu_sq));
        case d0P0:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + 3.0 * (vx2) + vf::basics::constant::c9o2 * (vx2) * (vx2)-cu_sq));
        case d0M0:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + 3.0 * (-vx2) + vf::basics::constant::c9o2 * (-vx2) * (-vx2) - cu_sq));
        case d00P:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + 3.0 * (vx3) + vf::basics::constant::c9o2 * (vx3) * (vx3)-cu_sq));
        case d00M:
            return REAL_CAST(vf::basics::constant::c2o27 * (drho + 3.0 * (-vx3) + vf::basics::constant::c9o2 * (-vx3) * (-vx3) - cu_sq));
        case dPP0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (vx1 + vx2) + vf::basics::constant::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
        case dMM0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (-vx1 - vx2) + vf::basics::constant::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
        case dPM0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (vx1 - vx2) + vf::basics::constant::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
        case dMP0:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (-vx1 + vx2) + vf::basics::constant::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
        case dP0P:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (vx1 + vx3) + vf::basics::constant::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
        case dM0M:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (-vx1 - vx3) + vf::basics::constant::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
        case dP0M:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (vx1 - vx3) + vf::basics::constant::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
        case dM0P:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (-vx1 + vx3) + vf::basics::constant::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
        case d0PP:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (vx2 + vx3) + vf::basics::constant::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
        case d0MM:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (-vx2 - vx3) + vf::basics::constant::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
        case d0PM:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (vx2 - vx3) + vf::basics::constant::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
        case d0MP:
            return REAL_CAST(vf::basics::constant::c1o54 *
                             (drho + 3.0 * (-vx2 + vx3) + vf::basics::constant::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
        case dPPP:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (vx1 + vx2 + vx3) +
                                               vf::basics::constant::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
        case dMMM:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (-vx1 - vx2 - vx3) +
                                               vf::basics::constant::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
        case dPPM:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (vx1 + vx2 - vx3) +
                                               vf::basics::constant::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
        case dMMP:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (-vx1 - vx2 + vx3) +
                                               vf::basics::constant::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
        case dPMP:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (vx1 - vx2 + vx3) +
                                               vf::basics::constant::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
        case dMPM:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (-vx1 + vx2 - vx3) +
                                               vf::basics::constant::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
        case dPMM:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (vx1 - vx2 - vx3) +
                                               vf::basics::constant::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
        case dMPP:
            return REAL_CAST(vf::basics::constant::c1o216 * (drho + 3.0 * (-vx1 + vx2 + vx3) +
                                               vf::basics::constant::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
        default:
            throw UbException(UB_EXARGS, "unknown dir");
    }
}
//////////////////////////////////////////////////////////////////////////
static void calcIncompFeq(real *const &feq /*[27]*/, const real &drho, const real &vx1, const real &vx2,
                          const real &vx3)
{
    using namespace vf::lbm::dir;

    real cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[d000] = vf::basics::constant::c8o27 * (drho - cu_sq);
    feq[dP00]    = vf::basics::constant::c2o27 * (drho + 3.0 * (vx1) + vf::basics::constant::c9o2 * (vx1) * (vx1)-cu_sq);
    feq[dM00]    = vf::basics::constant::c2o27 * (drho + 3.0 * (-vx1) + vf::basics::constant::c9o2 * (-vx1) * (-vx1) - cu_sq);
    feq[d0P0]    = vf::basics::constant::c2o27 * (drho + 3.0 * (vx2) + vf::basics::constant::c9o2 * (vx2) * (vx2)-cu_sq);
    feq[d0M0]    = vf::basics::constant::c2o27 * (drho + 3.0 * (-vx2) + vf::basics::constant::c9o2 * (-vx2) * (-vx2) - cu_sq);
    feq[d00P]    = vf::basics::constant::c2o27 * (drho + 3.0 * (vx3) + vf::basics::constant::c9o2 * (vx3) * (vx3)-cu_sq);
    feq[d00M]    = vf::basics::constant::c2o27 * (drho + 3.0 * (-vx3) + vf::basics::constant::c9o2 * (-vx3) * (-vx3) - cu_sq);
    feq[dPP0]   = vf::basics::constant::c1o54 * (drho + 3.0 * (vx1 + vx2) + vf::basics::constant::c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
    feq[dMM0]   = vf::basics::constant::c1o54 * (drho + 3.0 * (-vx1 - vx2) + vf::basics::constant::c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
    feq[dPM0]   = vf::basics::constant::c1o54 * (drho + 3.0 * (vx1 - vx2) + vf::basics::constant::c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
    feq[dMP0]   = vf::basics::constant::c1o54 * (drho + 3.0 * (-vx1 + vx2) + vf::basics::constant::c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
    feq[dP0P]   = vf::basics::constant::c1o54 * (drho + 3.0 * (vx1 + vx3) + vf::basics::constant::c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
    feq[dM0M]   = vf::basics::constant::c1o54 * (drho + 3.0 * (-vx1 - vx3) + vf::basics::constant::c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
    feq[dP0M]   = vf::basics::constant::c1o54 * (drho + 3.0 * (vx1 - vx3) + vf::basics::constant::c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
    feq[dM0P]   = vf::basics::constant::c1o54 * (drho + 3.0 * (-vx1 + vx3) + vf::basics::constant::c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
    feq[d0PP]   = vf::basics::constant::c1o54 * (drho + 3.0 * (vx2 + vx3) + vf::basics::constant::c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
    feq[d0MM]   = vf::basics::constant::c1o54 * (drho + 3.0 * (-vx2 - vx3) + vf::basics::constant::c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
    feq[d0PM]   = vf::basics::constant::c1o54 * (drho + 3.0 * (vx2 - vx3) + vf::basics::constant::c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
    feq[d0MP]   = vf::basics::constant::c1o54 * (drho + 3.0 * (-vx2 + vx3) + vf::basics::constant::c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
    feq[dPPP]  = vf::basics::constant::c1o216 *
               (drho + 3.0 * (vx1 + vx2 + vx3) + vf::basics::constant::c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
    feq[dMMM] = vf::basics::constant::c1o216 *
               (drho + 3.0 * (-vx1 - vx2 - vx3) + vf::basics::constant::c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
    feq[dPPM] = vf::basics::constant::c1o216 *
               (drho + 3.0 * (vx1 + vx2 - vx3) + vf::basics::constant::c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
    feq[dMMP] = vf::basics::constant::c1o216 *
               (drho + 3.0 * (-vx1 - vx2 + vx3) + vf::basics::constant::c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
    feq[dPMP] = vf::basics::constant::c1o216 *
               (drho + 3.0 * (vx1 - vx2 + vx3) + vf::basics::constant::c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
    feq[dMPM] = vf::basics::constant::c1o216 *
               (drho + 3.0 * (-vx1 + vx2 - vx3) + vf::basics::constant::c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
    feq[dPMM] = vf::basics::constant::c1o216 *
               (drho + 3.0 * (vx1 - vx2 - vx3) + vf::basics::constant::c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
    feq[dMPP] = vf::basics::constant::c1o216 *
               (drho + 3.0 * (-vx1 + vx2 + vx3) + vf::basics::constant::c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
}
//////////////////////////////////////////////////////////////////////////
static inline real getBoundaryVelocityForDirection(const int &direction, const real &bcVelocityX1,
                                                    const real &bcVelocityX2, const real &bcVelocityX3)
{
    using namespace vf::lbm::dir;
 
    switch (direction) {
        case dP00:
            return (real)(vf::basics::constant::c4o9 * (+bcVelocityX1));
        case dM00:
            return (real)(vf::basics::constant::c4o9 * (-bcVelocityX1));
        case d0P0:
            return (real)(vf::basics::constant::c4o9 * (+bcVelocityX2));
        case d0M0:
            return (real)(vf::basics::constant::c4o9 * (-bcVelocityX2));
        case d00P:
            return (real)(vf::basics::constant::c4o9 * (+bcVelocityX3));
        case d00M:
            return (real)(vf::basics::constant::c4o9 * (-bcVelocityX3));
        case dPP0:
            return (real)(vf::basics::constant::c1o9 * (+bcVelocityX1 + bcVelocityX2));
        case dMM0:
            return (real)(vf::basics::constant::c1o9 * (-bcVelocityX1 - bcVelocityX2));
        case dPM0:
            return (real)(vf::basics::constant::c1o9 * (+bcVelocityX1 - bcVelocityX2));
        case dMP0:
            return (real)(vf::basics::constant::c1o9 * (-bcVelocityX1 + bcVelocityX2));
        case dP0P:
            return (real)(vf::basics::constant::c1o9 * (+bcVelocityX1 + bcVelocityX3));
        case dM0M:
            return (real)(vf::basics::constant::c1o9 * (-bcVelocityX1 - bcVelocityX3));
        case dP0M:
            return (real)(vf::basics::constant::c1o9 * (+bcVelocityX1 - bcVelocityX3));
        case dM0P:
            return (real)(vf::basics::constant::c1o9 * (-bcVelocityX1 + bcVelocityX3));
        case d0PP:
            return (real)(vf::basics::constant::c1o9 * (+bcVelocityX2 + bcVelocityX3));
        case d0MM:
            return (real)(vf::basics::constant::c1o9 * (-bcVelocityX2 - bcVelocityX3));
        case d0PM:
            return (real)(vf::basics::constant::c1o9 * (+bcVelocityX2 - bcVelocityX3));
        case d0MP:
            return (real)(vf::basics::constant::c1o9 * (-bcVelocityX2 + bcVelocityX3));
        case dPPP:
            return (real)(vf::basics::constant::c1o36 * (+bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
        case dMMM:
            return (real)(vf::basics::constant::c1o36 * (-bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
        case dPPM:
            return (real)(vf::basics::constant::c1o36 * (+bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
        case dMMP:
            return (real)(vf::basics::constant::c1o36 * (-bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
        case dPMP:
            return (real)(vf::basics::constant::c1o36 * (+bcVelocityX1 - bcVelocityX2 + bcVelocityX3));
        case dMPM:
            return (real)(vf::basics::constant::c1o36 * (-bcVelocityX1 + bcVelocityX2 - bcVelocityX3));
        case dPMM:
            return (real)(vf::basics::constant::c1o36 * (+bcVelocityX1 - bcVelocityX2 - bcVelocityX3));
        case dMPP:
            return (real)(vf::basics::constant::c1o36 * (-bcVelocityX1 + bcVelocityX2 + bcVelocityX3));
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
static inline void calcDistanceToNeighbors(std::vector<real> &distNeigh, const real &deltaX1)
{
    using namespace vf::lbm::dir;

    // distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);

    distNeigh[dP00] = distNeigh[dM00] = distNeigh[d0P0] = deltaX1;
    distNeigh[d0M0] = distNeigh[d00P] = distNeigh[d00M] = deltaX1;
    distNeigh[dPP0] = distNeigh[dMP0] = distNeigh[dMM0] = distNeigh[dPM0] = vf::basics::constant::sqrt2 * deltaX1;
    distNeigh[dP0P] = distNeigh[d0PP] = distNeigh[dM0P] = distNeigh[d0MP] = vf::basics::constant::sqrt2 * deltaX1;
    distNeigh[dP0M] = distNeigh[d0PM] = distNeigh[dM0M] = distNeigh[d0MM] = vf::basics::constant::sqrt2 * deltaX1;
    distNeigh[dPPP] = distNeigh[dMPP] = distNeigh[dPMP] = distNeigh[dMMP] = vf::basics::constant::sqrt3 * deltaX1;
    distNeigh[dPPM] = distNeigh[dMPM] = distNeigh[dPMM] = distNeigh[dMMM] = vf::basics::constant::sqrt3 * deltaX1;
}
//////////////////////////////////////////////////////////////////////////
static inline void calcDistanceToNeighbors(std::vector<real> &distNeigh, const real &deltaX1, const real &deltaX2,
                                           const real &deltaX3)
{
    using namespace vf::lbm::dir;

    // distNeigh.resize(FENDDIR+1, UbMath::sqrt2*deltaX1);
    distNeigh[dP00] = distNeigh[dM00] = deltaX1;
    distNeigh[d0P0] = distNeigh[d0M0] = deltaX2;
    distNeigh[d00P] = distNeigh[d00M] = deltaX3;
    distNeigh[dPP0] = distNeigh[dMP0] = distNeigh[dMM0] = distNeigh[dPM0] = sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2);
    distNeigh[dP0P] = distNeigh[d0PP] = distNeigh[dM0P] = distNeigh[d0MP] = sqrt(deltaX1 * deltaX1 + deltaX3 * deltaX3);
    distNeigh[dP0M] = distNeigh[d0PM] = distNeigh[dM0M] = distNeigh[d0MM] = sqrt(deltaX2 * deltaX2 + deltaX3 * deltaX3);
    distNeigh[dPPP] = distNeigh[dMPP] = distNeigh[dPMP] = distNeigh[dMMP] =
        sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3);
    distNeigh[dPPM] = distNeigh[dMPM] = distNeigh[dPMM] = distNeigh[dMMM] =
        sqrt(deltaX1 * deltaX1 + deltaX2 * deltaX2 + deltaX3 * deltaX3);
}
//////////////////////////////////////////////////////////////////////////
static inline void initRayVectors(real *const &rayX1, real *const &rayX2, real *const &rayX3)
{
    using namespace vf::lbm::dir;

    int fdir;
    real c1oS2 = vf::basics::constant::one_over_sqrt2;
    real c1oS3 = vf::basics::constant::one_over_sqrt3;
    fdir         = dP00;
    rayX1[fdir]  = 1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = dM00;
    rayX1[fdir]  = -1.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 0.0;
    fdir         = d0P0;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 1.0;
    rayX3[fdir]  = 0.0;
    fdir         = d0M0;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -1.0;
    rayX3[fdir]  = 0.0;
    fdir         = d00P;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = 1.0;
    fdir         = d00M;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -1.0;
    fdir         = dPP0;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = dMM0;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = dPM0;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = dMP0;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = 0.0;
    fdir         = dP0P;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = dM0M;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = dP0M;
    rayX1[fdir]  = c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = -c1oS2;
    fdir         = dM0P;
    rayX1[fdir]  = -c1oS2;
    rayX2[fdir]  = 0.0;
    rayX3[fdir]  = c1oS2;
    fdir         = d0PP;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = c1oS2;
    fdir         = d0MM;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = d0PM;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = c1oS2;
    rayX3[fdir]  = -c1oS2;
    fdir         = d0MP;
    rayX1[fdir]  = 0.0;
    rayX2[fdir]  = -c1oS2;
    rayX3[fdir]  = c1oS2;
    fdir         = dPPP;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = dMPP;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = dPMP;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = dMMP;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = c1oS3;
    fdir         = dPPM;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = dMPM;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = dPMM;
    rayX1[fdir]  = c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = -c1oS3;
    fdir         = dMMM;
    rayX1[fdir]  = -c1oS3;
    rayX2[fdir]  = -c1oS3;
    rayX3[fdir]  = -c1oS3;
}
//////////////////////////////////////////////////////////////////////////
static inline real calcPress(const real *const f, real rho, real vx1, real vx2, real vx3)
{
    using namespace vf::lbm::dir;

    real op = 1.0;
    return ((f[dP00] + f[dM00] + f[d0P0] + f[d0M0] + f[d00P] + f[d00M] +
             2. * (f[dPP0] + f[dMM0] + f[dPM0] + f[dMP0] + f[dP0P] + f[dM0M] + f[dP0M] + f[dM0P] + f[d0PP] + f[d0MM] + f[d0PM] + f[d0MP]) +
             3. * (f[dPPP] + f[dMMP] + f[dPMP] + f[dMPP] + f[dPPM] + f[dMMM] + f[dPMM] + f[dMPM]) -
             (vx1 * vx1 + vx2 * vx2 + vx3 * vx3)) *
                (1 - 0.5 * op) +
            op * 0.5 * (rho)) *
           vf::basics::constant::c1o3;
}
//////////////////////////////////////////////////////////////////////////
static inline real getShearRate(const real *const f, real collFactorF)
{
    using namespace vf::lbm::dir;

    real mfcbb = f[dP00];
    real mfbcb = f[d0P0];
    real mfbbc = f[d00P];
    real mfccb = f[dPP0];
    real mfacb = f[dMP0];
    real mfcbc = f[dP0P];
    real mfabc = f[dM0P];
    real mfbcc = f[d0PP];
    real mfbac = f[d0MP];
    real mfccc = f[dPPP];
    real mfacc = f[dMPP];
    real mfcac = f[dPMP];
    real mfaac = f[dMMP];

    real mfabb = f[dM00];
    real mfbab = f[d0M0];
    real mfbba = f[d00M];
    real mfaab = f[dMM0];
    real mfcab = f[dPM0];
    real mfaba = f[dM0M];
    real mfcba = f[dP0M];
    real mfbaa = f[d0MM];
    real mfbca = f[d0PM];
    real mfaaa = f[dMMM];
    real mfcaa = f[dPMM];
    real mfaca = f[dMPM];
    real mfcca = f[dPPM];

    real mfbbb = f[d000];

    real m0, m1, m2;

    real rho = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) + (mfaab + mfacb + mfcab + mfccb) +
                  (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) +
                  (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

    real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) +
                   (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) + (mfcbb - mfabb));
    real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) +
                   (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) + (mfbcb - mfbab));
    real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) +
                   (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) + (mfbbc - mfbba));

    real oMdrho;

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

    real vx2;
    real vy2;
    real vz2;
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
    m0 += vf::basics::constant::c1o36 * oMdrho;
    mfaab = m1 - m0 * vvz;
    mfaac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaba + mfabc;
    m1    = mfabc - mfaba;
    m0    = m2 + mfabb;
    mfaba = m0;
    m0 += vf::basics::constant::c1o9 * oMdrho;
    mfabb = m1 - m0 * vvz;
    mfabc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfaca + mfacc;
    m1    = mfacc - mfaca;
    m0    = m2 + mfacb;
    mfaca = m0;
    m0 += vf::basics::constant::c1o36 * oMdrho;
    mfacb = m1 - m0 * vvz;
    mfacc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbaa + mfbac;
    m1    = mfbac - mfbaa;
    m0    = m2 + mfbab;
    mfbaa = m0;
    m0 += vf::basics::constant::c1o9 * oMdrho;
    mfbab = m1 - m0 * vvz;
    mfbac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbba + mfbbc;
    m1    = mfbbc - mfbba;
    m0    = m2 + mfbbb;
    mfbba = m0;
    m0 += vf::basics::constant::c4o9 * oMdrho;
    mfbbb = m1 - m0 * vvz;
    mfbbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbca + mfbcc;
    m1    = mfbcc - mfbca;
    m0    = m2 + mfbcb;
    mfbca = m0;
    m0 += vf::basics::constant::c1o9 * oMdrho;
    mfbcb = m1 - m0 * vvz;
    mfbcc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcaa + mfcac;
    m1    = mfcac - mfcaa;
    m0    = m2 + mfcab;
    mfcaa = m0;
    m0 += vf::basics::constant::c1o36 * oMdrho;
    mfcab = m1 - m0 * vvz;
    mfcac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcba + mfcbc;
    m1    = mfcbc - mfcba;
    m0    = m2 + mfcbb;
    mfcba = m0;
    m0 += vf::basics::constant::c1o9 * oMdrho;
    mfcbb = m1 - m0 * vvz;
    mfcbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcca + mfccc;
    m1    = mfccc - mfcca;
    m0    = m2 + mfccb;
    mfcca = m0;
    m0 += vf::basics::constant::c1o36 * oMdrho;
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
    m0 += vf::basics::constant::c1o6 * oMdrho;
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
    m0 += vf::basics::constant::c1o18 * oMdrho;
    mfabc = m1 - m0 * vvy;
    mfacc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfbaa + mfbca;
    m1    = mfbca - mfbaa;
    m0    = m2 + mfbba;
    mfbaa = m0;
    m0 += vf::basics::constant::c2o3 * oMdrho;
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
    m0 += vf::basics::constant::c2o9 * oMdrho;
    mfbbc = m1 - m0 * vvy;
    mfbcc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2    = mfcaa + mfcca;
    m1    = mfcca - mfcaa;
    m0    = m2 + mfcba;
    mfcaa = m0;
    m0 += vf::basics::constant::c1o6 * oMdrho;
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
    m0 += vf::basics::constant::c1o18 * oMdrho;
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
    m0 += vf::basics::constant::c1o3 * oMdrho;
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
    m0 += vf::basics::constant::c1o3 * oMdrho;
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
    m0 += vf::basics::constant::c1o9 * oMdrho;
    mfbcc = m1 - m0 * vvx;
    mfccc = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    // Cumulants
    ////////////////////////////////////////////////////////////////////////////////////
    real OxxPyyPzz = 1.; // omega2 or bulk viscosity

    real mxxPyyPzz = mfcaa + mfaca + mfaac;
    real mxxMyy    = mfcaa - mfaca;
    real mxxMzz    = mfcaa - mfaac;

    real dxux = -vf::basics::constant::c1o2 * collFactorF * (mxxMyy + mxxMzz) + vf::basics::constant::c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
    real dyuy = dxux + collFactorF * vf::basics::constant::c3o2 * mxxMyy;
    real dzuz = dxux + collFactorF * vf::basics::constant::c3o2 * mxxMzz;

    real Dxy = -vf::basics::constant::c3o1 * collFactorF * mfbba;
    real Dxz = -vf::basics::constant::c3o1 * collFactorF * mfbab;
    real Dyz = -vf::basics::constant::c3o1 * collFactorF * mfabb;

    return sqrt(vf::basics::constant::c2o1 * (dxux * dxux + dyuy * dyuy + dzuz * dzuz) + Dxy * Dxy + Dxz * Dxz + Dyz * Dyz) /
           (rho + vf::basics::constant::c1o1);
}

static inline std::array<real,6> getSecondMoments(const real *const f, real collFactorF)
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    real mfcbb = f[dP00];
    real mfbcb = f[d0P0];
    real mfbbc = f[d00P];
    real mfccb = f[dPP0];
    real mfacb = f[dMP0];
    real mfcbc = f[dP0P];
    real mfabc = f[dM0P];
    real mfbcc = f[d0PP];
    real mfbac = f[d0MP];
    real mfccc = f[dPPP];
    real mfacc = f[dMPP];
    real mfcac = f[dPMP];
    real mfaac = f[dMMP];

    real mfabb = f[dM00];
    real mfbab = f[d0M0];
    real mfbba = f[d00M];
    real mfaab = f[dMM0];
    real mfcab = f[dPM0];
    real mfaba = f[dM0M];
    real mfcba = f[dP0M];
    real mfbaa = f[d0MM];
    real mfbca = f[d0PM];
    real mfaaa = f[dMMM];
    real mfcaa = f[dPMM];
    real mfaca = f[dMPM];
    real mfcca = f[dPPM];

    real mfbbb = f[d000];

    real m0, m1, m2;

    //real rho = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) + (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

    real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) + (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) + (mfcbb - mfabb));
    real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) + (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) + (mfbcb - mfbab));
    real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) + (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) + (mfbbc - mfbba));

    real oMdrho;

    oMdrho = mfccc + mfaaa;
    m0 = mfaca + mfcac;
    m1 = mfacc + mfcaa;
    m2 = mfaac + mfcca;
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

    real vx2;
    real vy2;
    real vz2;
    vx2 = vvx * vvx;
    vy2 = vvy * vvy;
    vz2 = vvz * vvz;
    ////////////////////////////////////////////////////////////////////////////////////
    // Hin
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
    mfaac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaba + mfabc;
    m1 = mfabc - mfaba;
    m0 = m2 + mfabb;
    mfaba = m0;
    m0 += c1o9 * oMdrho;
    mfabb = m1 - m0 * vvz;
    mfabc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaca + mfacc;
    m1 = mfacc - mfaca;
    m0 = m2 + mfacb;
    mfaca = m0;
    m0 += c1o36 * oMdrho;
    mfacb = m1 - m0 * vvz;
    mfacc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbaa + mfbac;
    m1 = mfbac - mfbaa;
    m0 = m2 + mfbab;
    mfbaa = m0;
    m0 += c1o9 * oMdrho;
    mfbab = m1 - m0 * vvz;
    mfbac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbba + mfbbc;
    m1 = mfbbc - mfbba;
    m0 = m2 + mfbbb;
    mfbba = m0;
    m0 += c4o9 * oMdrho;
    mfbbb = m1 - m0 * vvz;
    mfbbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbca + mfbcc;
    m1 = mfbcc - mfbca;
    m0 = m2 + mfbcb;
    mfbca = m0;
    m0 += c1o9 * oMdrho;
    mfbcb = m1 - m0 * vvz;
    mfbcc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcaa + mfcac;
    m1 = mfcac - mfcaa;
    m0 = m2 + mfcab;
    mfcaa = m0;
    m0 += c1o36 * oMdrho;
    mfcab = m1 - m0 * vvz;
    mfcac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcba + mfcbc;
    m1 = mfcbc - mfcba;
    m0 = m2 + mfcbb;
    mfcba = m0;
    m0 += c1o9 * oMdrho;
    mfcbb = m1 - m0 * vvz;
    mfcbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcca + mfccc;
    m1 = mfccc - mfcca;
    m0 = m2 + mfccb;
    mfcca = m0;
    m0 += c1o36 * oMdrho;
    mfccb = m1 - m0 * vvz;
    mfccc = m2 - 2. * m1 * vvz + vz2 * m0;
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
    mfaca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaab + mfacb;
    m1 = mfacb - mfaab;
    m0 = m2 + mfabb;
    mfaab = m0;
    mfabb = m1 - m0 * vvy;
    mfacb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaac + mfacc;
    m1 = mfacc - mfaac;
    m0 = m2 + mfabc;
    mfaac = m0;
    m0 += c1o18 * oMdrho;
    mfabc = m1 - m0 * vvy;
    mfacc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbaa + mfbca;
    m1 = mfbca - mfbaa;
    m0 = m2 + mfbba;
    mfbaa = m0;
    m0 += c2o3 * oMdrho;
    mfbba = m1 - m0 * vvy;
    mfbca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbab + mfbcb;
    m1 = mfbcb - mfbab;
    m0 = m2 + mfbbb;
    mfbab = m0;
    mfbbb = m1 - m0 * vvy;
    mfbcb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbac + mfbcc;
    m1 = mfbcc - mfbac;
    m0 = m2 + mfbbc;
    mfbac = m0;
    m0 += c2o9 * oMdrho;
    mfbbc = m1 - m0 * vvy;
    mfbcc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcaa + mfcca;
    m1 = mfcca - mfcaa;
    m0 = m2 + mfcba;
    mfcaa = m0;
    m0 += c1o6 * oMdrho;
    mfcba = m1 - m0 * vvy;
    mfcca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcab + mfccb;
    m1 = mfccb - mfcab;
    m0 = m2 + mfcbb;
    mfcab = m0;
    mfcbb = m1 - m0 * vvy;
    mfccb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcac + mfccc;
    m1 = mfccc - mfcac;
    m0 = m2 + mfcbc;
    mfcac = m0;
    m0 += c1o18 * oMdrho;
    mfcbc = m1 - m0 * vvy;
    mfccc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    m2 = mfaaa + mfcaa;
    m1 = mfcaa - mfaaa;
    m0 = m2 + mfbaa;
    mfaaa = m0;
    m0 += 1. * oMdrho;
    mfbaa = m1 - m0 * vvx;
    mfcaa = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaba + mfcba;
    m1 = mfcba - mfaba;
    m0 = m2 + mfbba;
    mfaba = m0;
    mfbba = m1 - m0 * vvx;
    mfcba = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaca + mfcca;
    m1 = mfcca - mfaca;
    m0 = m2 + mfbca;
    mfaca = m0;
    m0 += c1o3 * oMdrho;
    mfbca = m1 - m0 * vvx;
    mfcca = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaab + mfcab;
    m1 = mfcab - mfaab;
    m0 = m2 + mfbab;
    mfaab = m0;
    mfbab = m1 - m0 * vvx;
    mfcab = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfabb + mfcbb;
    m1 = mfcbb - mfabb;
    m0 = m2 + mfbbb;
    mfabb = m0;
    mfbbb = m1 - m0 * vvx;
    mfcbb = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfacb + mfccb;
    m1 = mfccb - mfacb;
    m0 = m2 + mfbcb;
    mfacb = m0;
    mfbcb = m1 - m0 * vvx;
    mfccb = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaac + mfcac;
    m1 = mfcac - mfaac;
    m0 = m2 + mfbac;
    mfaac = m0;
    m0 += c1o3 * oMdrho;
    mfbac = m1 - m0 * vvx;
    mfcac = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfabc + mfcbc;
    m1 = mfcbc - mfabc;
    m0 = m2 + mfbbc;
    mfabc = m0;
    mfbbc = m1 - m0 * vvx;
    mfcbc = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfacc + mfccc;
    m1 = mfccc - mfacc;
    m0 = m2 + mfbcc;
    mfacc = m0;
    m0 += c1o9 * oMdrho;
    mfbcc = m1 - m0 * vvx;
    mfccc = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    // Cumulants
    ////////////////////////////////////////////////////////////////////////////////////
    real OxxPyyPzz = 1.; // omega2 or bulk viscosity

    real mxxPyyPzz = mfcaa + mfaca + mfaac;
    real mxxMyy = mfcaa - mfaca;
    real mxxMzz = mfcaa - mfaac;

   // average pre and post collision
    std::array<real, 6> moments = {
    (mxxPyyPzz-mfaaa) * (c1o1 - c1o2 * OxxPyyPzz),
    (mxxMyy) * (c1o1 - c1o2 * collFactorF),
    (mxxMzz) * (c1o1 - c1o2 * collFactorF),
    (mfbba)  * (c1o1 - c1o2 * collFactorF),
    (mfbab)  * (c1o1 - c1o2 * collFactorF),
    (mfabb)  * (c1o1 - c1o2 * collFactorF)
    };

    return moments;
}
static inline std::array<real, 6> getStressTensor(const real *const f, real collFactorF)
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    real mfcbb = f[dP00];
    real mfbcb = f[d0P0];
    real mfbbc = f[d00P];
    real mfccb = f[dPP0];
    real mfacb = f[dMP0];
    real mfcbc = f[dP0P];
    real mfabc = f[dM0P];
    real mfbcc = f[d0PP];
    real mfbac = f[d0MP];
    real mfccc = f[dPPP];
    real mfacc = f[dMPP];
    real mfcac = f[dPMP];
    real mfaac = f[dMMP];

    real mfabb = f[dM00];
    real mfbab = f[d0M0];
    real mfbba = f[d00M];
    real mfaab = f[dMM0];
    real mfcab = f[dPM0];
    real mfaba = f[dM0M];
    real mfcba = f[dP0M];
    real mfbaa = f[d0MM];
    real mfbca = f[d0PM];
    real mfaaa = f[dMMM];
    real mfcaa = f[dPMM];
    real mfaca = f[dMPM];
    real mfcca = f[dPPM];

    real mfbbb = f[d000];

    real m0, m1, m2;

    //real rho = (mfaaa + mfaac + mfaca + mfcaa + mfacc + mfcac + mfccc + mfcca) + (mfaab + mfacb + mfcab + mfccb) + (mfaba + mfabc + mfcba + mfcbc) + (mfbaa + mfbac + mfbca + mfbcc) + (mfabb + mfcbb) + (mfbab + mfbcb) + (mfbba + mfbbc) + mfbbb;

    real vvx = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfcaa - mfacc) + (mfcca - mfaac))) + (((mfcba - mfabc) + (mfcbc - mfaba)) + ((mfcab - mfacb) + (mfccb - mfaab))) + (mfcbb - mfabb));
    real vvy = ((((mfccc - mfaaa) + (mfaca - mfcac)) + ((mfacc - mfcaa) + (mfcca - mfaac))) + (((mfbca - mfbac) + (mfbcc - mfbaa)) + ((mfacb - mfcab) + (mfccb - mfaab))) + (mfbcb - mfbab));
    real vvz = ((((mfccc - mfaaa) + (mfcac - mfaca)) + ((mfacc - mfcaa) + (mfaac - mfcca))) + (((mfbac - mfbca) + (mfbcc - mfbaa)) + ((mfabc - mfcba) + (mfcbc - mfaba))) + (mfbbc - mfbba));

    real oMdrho;

    oMdrho = mfccc + mfaaa;
    m0 = mfaca + mfcac;
    m1 = mfacc + mfcaa;
    m2 = mfaac + mfcca;
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

    real vx2;
    real vy2;
    real vz2;
    vx2 = vvx * vvx;
    vy2 = vvy * vvy;
    vz2 = vvz * vvz;
    ////////////////////////////////////////////////////////////////////////////////////
    // Hin
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
    mfaac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaba + mfabc;
    m1 = mfabc - mfaba;
    m0 = m2 + mfabb;
    mfaba = m0;
    m0 += c1o9 * oMdrho;
    mfabb = m1 - m0 * vvz;
    mfabc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaca + mfacc;
    m1 = mfacc - mfaca;
    m0 = m2 + mfacb;
    mfaca = m0;
    m0 += c1o36 * oMdrho;
    mfacb = m1 - m0 * vvz;
    mfacc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbaa + mfbac;
    m1 = mfbac - mfbaa;
    m0 = m2 + mfbab;
    mfbaa = m0;
    m0 += c1o9 * oMdrho;
    mfbab = m1 - m0 * vvz;
    mfbac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbba + mfbbc;
    m1 = mfbbc - mfbba;
    m0 = m2 + mfbbb;
    mfbba = m0;
    m0 += c4o9 * oMdrho;
    mfbbb = m1 - m0 * vvz;
    mfbbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbca + mfbcc;
    m1 = mfbcc - mfbca;
    m0 = m2 + mfbcb;
    mfbca = m0;
    m0 += c1o9 * oMdrho;
    mfbcb = m1 - m0 * vvz;
    mfbcc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcaa + mfcac;
    m1 = mfcac - mfcaa;
    m0 = m2 + mfcab;
    mfcaa = m0;
    m0 += c1o36 * oMdrho;
    mfcab = m1 - m0 * vvz;
    mfcac = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcba + mfcbc;
    m1 = mfcbc - mfcba;
    m0 = m2 + mfcbb;
    mfcba = m0;
    m0 += c1o9 * oMdrho;
    mfcbb = m1 - m0 * vvz;
    mfcbc = m2 - 2. * m1 * vvz + vz2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcca + mfccc;
    m1 = mfccc - mfcca;
    m0 = m2 + mfccb;
    mfcca = m0;
    m0 += c1o36 * oMdrho;
    mfccb = m1 - m0 * vvz;
    mfccc = m2 - 2. * m1 * vvz + vz2 * m0;
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
    mfaca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaab + mfacb;
    m1 = mfacb - mfaab;
    m0 = m2 + mfabb;
    mfaab = m0;
    mfabb = m1 - m0 * vvy;
    mfacb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaac + mfacc;
    m1 = mfacc - mfaac;
    m0 = m2 + mfabc;
    mfaac = m0;
    m0 += c1o18 * oMdrho;
    mfabc = m1 - m0 * vvy;
    mfacc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbaa + mfbca;
    m1 = mfbca - mfbaa;
    m0 = m2 + mfbba;
    mfbaa = m0;
    m0 += c2o3 * oMdrho;
    mfbba = m1 - m0 * vvy;
    mfbca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbab + mfbcb;
    m1 = mfbcb - mfbab;
    m0 = m2 + mfbbb;
    mfbab = m0;
    mfbbb = m1 - m0 * vvy;
    mfbcb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfbac + mfbcc;
    m1 = mfbcc - mfbac;
    m0 = m2 + mfbbc;
    mfbac = m0;
    m0 += c2o9 * oMdrho;
    mfbbc = m1 - m0 * vvy;
    mfbcc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcaa + mfcca;
    m1 = mfcca - mfcaa;
    m0 = m2 + mfcba;
    mfcaa = m0;
    m0 += c1o6 * oMdrho;
    mfcba = m1 - m0 * vvy;
    mfcca = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcab + mfccb;
    m1 = mfccb - mfcab;
    m0 = m2 + mfcbb;
    mfcab = m0;
    mfcbb = m1 - m0 * vvy;
    mfccb = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfcac + mfccc;
    m1 = mfccc - mfcac;
    m0 = m2 + mfcbc;
    mfcac = m0;
    m0 += c1o18 * oMdrho;
    mfcbc = m1 - m0 * vvy;
    mfccc = m2 - 2. * m1 * vvy + vy2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    // mit     1, 0, 1/3, 0, 0, 0, 1/3, 0, 1/9            Konditionieren
    ////////////////////////////////////////////////////////////////////////////////////
    // X - Dir
    m2 = mfaaa + mfcaa;
    m1 = mfcaa - mfaaa;
    m0 = m2 + mfbaa;
    mfaaa = m0;
    m0 += 1. * oMdrho;
    mfbaa = m1 - m0 * vvx;
    mfcaa = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaba + mfcba;
    m1 = mfcba - mfaba;
    m0 = m2 + mfbba;
    mfaba = m0;
    mfbba = m1 - m0 * vvx;
    mfcba = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaca + mfcca;
    m1 = mfcca - mfaca;
    m0 = m2 + mfbca;
    mfaca = m0;
    m0 += c1o3 * oMdrho;
    mfbca = m1 - m0 * vvx;
    mfcca = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaab + mfcab;
    m1 = mfcab - mfaab;
    m0 = m2 + mfbab;
    mfaab = m0;
    mfbab = m1 - m0 * vvx;
    mfcab = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfabb + mfcbb;
    m1 = mfcbb - mfabb;
    m0 = m2 + mfbbb;
    mfabb = m0;
    mfbbb = m1 - m0 * vvx;
    mfcbb = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfacb + mfccb;
    m1 = mfccb - mfacb;
    m0 = m2 + mfbcb;
    mfacb = m0;
    mfbcb = m1 - m0 * vvx;
    mfccb = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfaac + mfcac;
    m1 = mfcac - mfaac;
    m0 = m2 + mfbac;
    mfaac = m0;
    m0 += c1o3 * oMdrho;
    mfbac = m1 - m0 * vvx;
    mfcac = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfabc + mfcbc;
    m1 = mfcbc - mfabc;
    m0 = m2 + mfbbc;
    mfabc = m0;
    mfbbc = m1 - m0 * vvx;
    mfcbc = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    m2 = mfacc + mfccc;
    m1 = mfccc - mfacc;
    m0 = m2 + mfbcc;
    mfacc = m0;
    m0 += c1o9 * oMdrho;
    mfbcc = m1 - m0 * vvx;
    mfccc = m2 - 2. * m1 * vvx + vx2 * m0;
    ////////////////////////////////////////////////////////////////////////////////////
    // Cumulants
    ////////////////////////////////////////////////////////////////////////////////////
    real OxxPyyPzz = 1.; // omega2 or bulk viscosity

    real mxxPyyPzz = mfcaa + mfaca + mfaac;
    real mxxMyy = mfcaa - mfaca;
    real mxxMzz = mfcaa - mfaac;

    real dxux = -c1o2 * collFactorF * (mxxMyy + mxxMzz) + c1o2 * OxxPyyPzz * (mfaaa - mxxPyyPzz);
    real dyuy = dxux + collFactorF * c3o2 * mxxMyy;
    real dzuz = dxux + collFactorF * c3o2 * mxxMzz;

    real Dxy = -c3o1 * collFactorF * mfbba;
    real Dxz = -c3o1 * collFactorF * mfbab;
    real Dyz = -c3o1 * collFactorF * mfabb;
    real nu = c1o3 * (c1o1 / collFactorF - c1o2);

    // average pre and post collision
    std::array<real, 6> moments = { -c1o3 * mfaaa + c2o1*nu*dxux,
                                    -c1o3 * mfaaa + c2o1*nu*dyuy,
                                    -c1o3 * mfaaa + c2o1*nu*dzuz,
                                     nu*Dxy,nu*Dxz,nu*Dyz};

    return moments;
}
//Multiphase stuff
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseFeq(real *const &feq /*[27]*/, const real &rho, const real &p1, const real &vx1,
                              const real &vx2, const real &vx3)
{
    using namespace vf::lbm::dir;

    using namespace vf::basics::constant;
    real cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[d000] = c8o27 * (p1 + rho * c1o3 * (-cu_sq));
    feq[dP00]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq));
    feq[dM00]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[d0P0]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq));
    feq[d0M0]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[d00P]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq));
    feq[d00M]    = c2o27 * (p1 + rho * c1o3 * (3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[dPP0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[dMM0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[dPM0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[dMP0]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[dP0P]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[dM0M]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[dP0M]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[dM0P]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[d0PP]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[d0MM]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[d0PM]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[d0MP]   = c1o54 * (p1 + rho * c1o3 * (3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[dPPP] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[dMMM] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[dPPM] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[dMMP] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[dPMP] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[dMPM] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[dPMM] =
        c1o216 * (p1 + rho * c1o3 * (3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[dMPP] = c1o216 *
               (p1 + rho * c1o3 * (3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseFeqVB(real *const &feq /*[27]*/, const real &p1, const real &vx1, const real &vx2,
                                const real &vx3)
{
    using namespace vf::lbm::dir;

    using namespace vf::basics::constant;
    real cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    feq[d000] = p1 + c8o27 * (-cu_sq);
    feq[dP00]    = c2o27 * ((3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq));
    feq[dM00]    = c2o27 * ((3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq));
    feq[d0P0]    = c2o27 * ((3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq));
    feq[d0M0]    = c2o27 * ((3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq));
    feq[d00P]    = c2o27 * ((3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq));
    feq[d00M]    = c2o27 * ((3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq));
    feq[dPP0]   = c1o54 * ((3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq));
    feq[dMM0]   = c1o54 * ((3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq));
    feq[dPM0]   = c1o54 * ((3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq));
    feq[dMP0]   = c1o54 * ((3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq));
    feq[dP0P]   = c1o54 * ((3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq));
    feq[dM0M]   = c1o54 * ((3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq));
    feq[dP0M]   = c1o54 * ((3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq));
    feq[dM0P]   = c1o54 * ((3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq));
    feq[d0PP]   = c1o54 * ((3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq));
    feq[d0MM]   = c1o54 * ((3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq));
    feq[d0PM]   = c1o54 * ((3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq));
    feq[d0MP]   = c1o54 * ((3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq));
    feq[dPPP]  = c1o216 * ((3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq));
    feq[dMMM]  = c1o216 * ((3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq));
    feq[dPPM]  = c1o216 * ((3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq));
    feq[dMMP]  = c1o216 * ((3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq));
    feq[dPMP]  = c1o216 * ((3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq));
    feq[dMPM]  = c1o216 * ((3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq));
    feq[dPMM]  = c1o216 * ((3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq));
    feq[dMPP]  = c1o216 * ((3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq));
}
//////////////////////////////////////////////////////////////////////////
static void calcMultiphaseHeq(real *const &heq /*[27]*/, const real &phi, const real &vx1, const real &vx2,
                              const real &vx3)
{
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    real cu_sq = 1.5 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

    heq[d000] = c8o27 * phi * (1.0 - cu_sq);
    heq[dP00]    = c2o27 * phi * (1.0 + 3.0 * (vx1) + c9o2 * (vx1) * (vx1)-cu_sq);
    heq[dM00]    = c2o27 * phi * (1.0 + 3.0 * (-vx1) + c9o2 * (-vx1) * (-vx1) - cu_sq);
    heq[d0P0]    = c2o27 * phi * (1.0 + 3.0 * (vx2) + c9o2 * (vx2) * (vx2)-cu_sq);
    heq[d0M0]    = c2o27 * phi * (1.0 + 3.0 * (-vx2) + c9o2 * (-vx2) * (-vx2) - cu_sq);
    heq[d00P]    = c2o27 * phi * (1.0 + 3.0 * (vx3) + c9o2 * (vx3) * (vx3)-cu_sq);
    heq[d00M]    = c2o27 * phi * (1.0 + 3.0 * (-vx3) + c9o2 * (-vx3) * (-vx3) - cu_sq);
    heq[dPP0]   = c1o54 * phi * (1.0 + 3.0 * (vx1 + vx2) + c9o2 * (vx1 + vx2) * (vx1 + vx2) - cu_sq);
    heq[dMM0]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 - vx2) + c9o2 * (-vx1 - vx2) * (-vx1 - vx2) - cu_sq);
    heq[dPM0]   = c1o54 * phi * (1.0 + 3.0 * (vx1 - vx2) + c9o2 * (vx1 - vx2) * (vx1 - vx2) - cu_sq);
    heq[dMP0]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 + vx2) + c9o2 * (-vx1 + vx2) * (-vx1 + vx2) - cu_sq);
    heq[dP0P]   = c1o54 * phi * (1.0 + 3.0 * (vx1 + vx3) + c9o2 * (vx1 + vx3) * (vx1 + vx3) - cu_sq);
    heq[dM0M]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 - vx3) + c9o2 * (-vx1 - vx3) * (-vx1 - vx3) - cu_sq);
    heq[dP0M]   = c1o54 * phi * (1.0 + 3.0 * (vx1 - vx3) + c9o2 * (vx1 - vx3) * (vx1 - vx3) - cu_sq);
    heq[dM0P]   = c1o54 * phi * (1.0 + 3.0 * (-vx1 + vx3) + c9o2 * (-vx1 + vx3) * (-vx1 + vx3) - cu_sq);
    heq[d0PP]   = c1o54 * phi * (1.0 + 3.0 * (vx2 + vx3) + c9o2 * (vx2 + vx3) * (vx2 + vx3) - cu_sq);
    heq[d0MM]   = c1o54 * phi * (1.0 + 3.0 * (-vx2 - vx3) + c9o2 * (-vx2 - vx3) * (-vx2 - vx3) - cu_sq);
    heq[d0PM]   = c1o54 * phi * (1.0 + 3.0 * (vx2 - vx3) + c9o2 * (vx2 - vx3) * (vx2 - vx3) - cu_sq);
    heq[d0MP]   = c1o54 * phi * (1.0 + 3.0 * (-vx2 + vx3) + c9o2 * (-vx2 + vx3) * (-vx2 + vx3) - cu_sq);
    heq[dPPP]  = c1o216 * phi * (1.0 + 3.0 * (vx1 + vx2 + vx3) + c9o2 * (vx1 + vx2 + vx3) * (vx1 + vx2 + vx3) - cu_sq);
    heq[dMMM] = c1o216 * phi * (1.0 + 3.0 * (-vx1 - vx2 - vx3) + c9o2 * (-vx1 - vx2 - vx3) * (-vx1 - vx2 - vx3) - cu_sq);
    heq[dPPM] = c1o216 * phi * (1.0 + 3.0 * (vx1 + vx2 - vx3) + c9o2 * (vx1 + vx2 - vx3) * (vx1 + vx2 - vx3) - cu_sq);
    heq[dMMP] = c1o216 * phi * (1.0 + 3.0 * (-vx1 - vx2 + vx3) + c9o2 * (-vx1 - vx2 + vx3) * (-vx1 - vx2 + vx3) - cu_sq);
    heq[dPMP] = c1o216 * phi * (1.0 + 3.0 * (vx1 - vx2 + vx3) + c9o2 * (vx1 - vx2 + vx3) * (vx1 - vx2 + vx3) - cu_sq);
    heq[dMPM] = c1o216 * phi * (1.0 + 3.0 * (-vx1 + vx2 - vx3) + c9o2 * (-vx1 + vx2 - vx3) * (-vx1 + vx2 - vx3) - cu_sq);
    heq[dPMM] = c1o216 * phi * (1.0 + 3.0 * (vx1 - vx2 - vx3) + c9o2 * (vx1 - vx2 - vx3) * (vx1 - vx2 - vx3) - cu_sq);
    heq[dMPP] = c1o216 * phi * (1.0 + 3.0 * (-vx1 + vx2 + vx3) + c9o2 * (-vx1 + vx2 + vx3) * (-vx1 + vx2 + vx3) - cu_sq);
}

} // namespace D3Q27System
#endif
