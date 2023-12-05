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
//! \author Martin Schoenherr
//=======================================================================================
#include "LBM/LB.h" 
#include "lbm/constants/D3Q27.h"

#include <basics/constants/NumericConstants.h>

using namespace vf::basics::constant;
using namespace vf::lbm::dir;

////////////////////////////////////////////////////////////////////////////////
inline __device__ real calcDistributionBC_AD_interpol(real q, real weight, real v, real v_sq, real f, real finf, real omegaDiffusivity, real jTangential, real concentration) {
    real feq = weight * concentration * (c1o1 + c3o1 * v + c9o2 * v * v * concentration - v_sq * concentration);
    return (c1o1 - q) / (c1o1 + q) * ((f - feq * omegaDiffusivity) / (c1o1 - omegaDiffusivity)) + (q * (f + finf) - c6o1 * weight * (jTangential)) / (c1o1 + q);
}
////////////////////////////////////////////////////////////////////////////////
inline __device__ real calcDistributionBC_AD(real q, real weight, real v, real v_sq, real f, real finf, real omegaDiffusivity, real jTangential, real concentration) {
    return f - c6o1 * weight * jTangential;
}


// has to be excecuted before Fluid BCs
//////////////////////////////////////////////////////////////////////////////
__global__ void AdvectionDiffusionSlipVelocityCompressible_Device(
    real *normalX,
    real *normalY,
    real *normalZ,
    real *distributions,
    real *distributionsAD,
    int *QindexArray,
    real *Qarrays,
    uint numberOfBCnodes,
    real omegaDiffusivity,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    Distributions27 D;
    if (isEvenTimestep)
    {
        D.f[dP00] = &distributions[dP00 * numberOfLBnodes];
        D.f[dM00] = &distributions[dM00 * numberOfLBnodes];
        D.f[d0P0] = &distributions[d0P0 * numberOfLBnodes];
        D.f[d0M0] = &distributions[d0M0 * numberOfLBnodes];
        D.f[d00P] = &distributions[d00P * numberOfLBnodes];
        D.f[d00M] = &distributions[d00M * numberOfLBnodes];
        D.f[dPP0] = &distributions[dPP0 * numberOfLBnodes];
        D.f[dMM0] = &distributions[dMM0 * numberOfLBnodes];
        D.f[dPM0] = &distributions[dPM0 * numberOfLBnodes];
        D.f[dMP0] = &distributions[dMP0 * numberOfLBnodes];
        D.f[dP0P] = &distributions[dP0P * numberOfLBnodes];
        D.f[dM0M] = &distributions[dM0M * numberOfLBnodes];
        D.f[dP0M] = &distributions[dP0M * numberOfLBnodes];
        D.f[dM0P] = &distributions[dM0P * numberOfLBnodes];
        D.f[d0PP] = &distributions[d0PP * numberOfLBnodes];
        D.f[d0MM] = &distributions[d0MM * numberOfLBnodes];
        D.f[d0PM] = &distributions[d0PM * numberOfLBnodes];
        D.f[d0MP] = &distributions[d0MP * numberOfLBnodes];
        D.f[d000] = &distributions[d000 * numberOfLBnodes];
        D.f[dPPP] = &distributions[dPPP * numberOfLBnodes];
        D.f[dMMP] = &distributions[dMMP * numberOfLBnodes];
        D.f[dPMP] = &distributions[dPMP * numberOfLBnodes];
        D.f[dMPP] = &distributions[dMPP * numberOfLBnodes];
        D.f[dPPM] = &distributions[dPPM * numberOfLBnodes];
        D.f[dMMM] = &distributions[dMMM * numberOfLBnodes];
        D.f[dPMM] = &distributions[dPMM * numberOfLBnodes];
        D.f[dMPM] = &distributions[dMPM * numberOfLBnodes];
    }
    else
    {
        D.f[dM00] = &distributions[dP00 * numberOfLBnodes];
        D.f[dP00] = &distributions[dM00 * numberOfLBnodes];
        D.f[d0M0] = &distributions[d0P0 * numberOfLBnodes];
        D.f[d0P0] = &distributions[d0M0 * numberOfLBnodes];
        D.f[d00M] = &distributions[d00P * numberOfLBnodes];
        D.f[d00P] = &distributions[d00M * numberOfLBnodes];
        D.f[dMM0] = &distributions[dPP0 * numberOfLBnodes];
        D.f[dPP0] = &distributions[dMM0 * numberOfLBnodes];
        D.f[dMP0] = &distributions[dPM0 * numberOfLBnodes];
        D.f[dPM0] = &distributions[dMP0 * numberOfLBnodes];
        D.f[dM0M] = &distributions[dP0P * numberOfLBnodes];
        D.f[dP0P] = &distributions[dM0M * numberOfLBnodes];
        D.f[dM0P] = &distributions[dP0M * numberOfLBnodes];
        D.f[dP0M] = &distributions[dM0P * numberOfLBnodes];
        D.f[d0MM] = &distributions[d0PP * numberOfLBnodes];
        D.f[d0PP] = &distributions[d0MM * numberOfLBnodes];
        D.f[d0MP] = &distributions[d0PM * numberOfLBnodes];
        D.f[d0PM] = &distributions[d0MP * numberOfLBnodes];
        D.f[d000] = &distributions[d000 * numberOfLBnodes];
        D.f[dPPP] = &distributions[dMMM * numberOfLBnodes];
        D.f[dMMP] = &distributions[dPPM * numberOfLBnodes];
        D.f[dPMP] = &distributions[dMPM * numberOfLBnodes];
        D.f[dMPP] = &distributions[dPMM * numberOfLBnodes];
        D.f[dPPM] = &distributions[dMMP * numberOfLBnodes];
        D.f[dMMM] = &distributions[dPPP * numberOfLBnodes];
        D.f[dPMM] = &distributions[dMPP * numberOfLBnodes];
        D.f[dMPM] = &distributions[dPMP * numberOfLBnodes];
    }
    ////////////////////////////////////////////////////////////////////////////////
    Distributions27 DAD;
    if (isEvenTimestep)
    {
        DAD.f[dP00] = &distributionsAD[dP00 * numberOfLBnodes];
        DAD.f[dM00] = &distributionsAD[dM00 * numberOfLBnodes];
        DAD.f[d0P0] = &distributionsAD[d0P0 * numberOfLBnodes];
        DAD.f[d0M0] = &distributionsAD[d0M0 * numberOfLBnodes];
        DAD.f[d00P] = &distributionsAD[d00P * numberOfLBnodes];
        DAD.f[d00M] = &distributionsAD[d00M * numberOfLBnodes];
        DAD.f[dPP0] = &distributionsAD[dPP0 * numberOfLBnodes];
        DAD.f[dMM0] = &distributionsAD[dMM0 * numberOfLBnodes];
        DAD.f[dPM0] = &distributionsAD[dPM0 * numberOfLBnodes];
        DAD.f[dMP0] = &distributionsAD[dMP0 * numberOfLBnodes];
        DAD.f[dP0P] = &distributionsAD[dP0P * numberOfLBnodes];
        DAD.f[dM0M] = &distributionsAD[dM0M * numberOfLBnodes];
        DAD.f[dP0M] = &distributionsAD[dP0M * numberOfLBnodes];
        DAD.f[dM0P] = &distributionsAD[dM0P * numberOfLBnodes];
        DAD.f[d0PP] = &distributionsAD[d0PP * numberOfLBnodes];
        DAD.f[d0MM] = &distributionsAD[d0MM * numberOfLBnodes];
        DAD.f[d0PM] = &distributionsAD[d0PM * numberOfLBnodes];
        DAD.f[d0MP] = &distributionsAD[d0MP * numberOfLBnodes];
        DAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
        DAD.f[dPPP] = &distributionsAD[dPPP * numberOfLBnodes];
        DAD.f[dMMP] = &distributionsAD[dMMP * numberOfLBnodes];
        DAD.f[dPMP] = &distributionsAD[dPMP * numberOfLBnodes];
        DAD.f[dMPP] = &distributionsAD[dMPP * numberOfLBnodes];
        DAD.f[dPPM] = &distributionsAD[dPPM * numberOfLBnodes];
        DAD.f[dMMM] = &distributionsAD[dMMM * numberOfLBnodes];
        DAD.f[dPMM] = &distributionsAD[dPMM * numberOfLBnodes];
        DAD.f[dMPM] = &distributionsAD[dMPM * numberOfLBnodes];
    }
    else
    {
        DAD.f[dM00] = &distributionsAD[dP00 * numberOfLBnodes];
        DAD.f[dP00] = &distributionsAD[dM00 * numberOfLBnodes];
        DAD.f[d0M0] = &distributionsAD[d0P0 * numberOfLBnodes];
        DAD.f[d0P0] = &distributionsAD[d0M0 * numberOfLBnodes];
        DAD.f[d00M] = &distributionsAD[d00P * numberOfLBnodes];
        DAD.f[d00P] = &distributionsAD[d00M * numberOfLBnodes];
        DAD.f[dMM0] = &distributionsAD[dPP0 * numberOfLBnodes];
        DAD.f[dPP0] = &distributionsAD[dMM0 * numberOfLBnodes];
        DAD.f[dMP0] = &distributionsAD[dPM0 * numberOfLBnodes];
        DAD.f[dPM0] = &distributionsAD[dMP0 * numberOfLBnodes];
        DAD.f[dM0M] = &distributionsAD[dP0P * numberOfLBnodes];
        DAD.f[dP0P] = &distributionsAD[dM0M * numberOfLBnodes];
        DAD.f[dM0P] = &distributionsAD[dP0M * numberOfLBnodes];
        DAD.f[dP0M] = &distributionsAD[dM0P * numberOfLBnodes];
        DAD.f[d0MM] = &distributionsAD[d0PP * numberOfLBnodes];
        DAD.f[d0PP] = &distributionsAD[d0MM * numberOfLBnodes];
        DAD.f[d0MP] = &distributionsAD[d0PM * numberOfLBnodes];
        DAD.f[d0PM] = &distributionsAD[d0MP * numberOfLBnodes];
        DAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
        DAD.f[dPPP] = &distributionsAD[dMMM * numberOfLBnodes];
        DAD.f[dMMP] = &distributionsAD[dPPM * numberOfLBnodes];
        DAD.f[dPMP] = &distributionsAD[dMPM * numberOfLBnodes];
        DAD.f[dMPP] = &distributionsAD[dPMM * numberOfLBnodes];
        DAD.f[dPPM] = &distributionsAD[dMMP * numberOfLBnodes];
        DAD.f[dMMM] = &distributionsAD[dPPP * numberOfLBnodes];
        DAD.f[dPMM] = &distributionsAD[dMPP * numberOfLBnodes];
        DAD.f[dMPM] = &distributionsAD[dPMP * numberOfLBnodes];
    }
    ////////////////////////////////////////////////////////////////////////////////
    const unsigned  x = threadIdx.x;
    const unsigned  y = blockIdx.x; 
    const unsigned  z = blockIdx.y; 

    const unsigned nx = blockDim.x;
    const unsigned ny = gridDim.x;

    const unsigned k = nx * (ny * z + y) + x;
    //////////////////////////////////////////////////////////////////////////

    if (k < numberOfBCnodes)
    {
        ////////////////////////////////////////////////////////////////////////////////
        real NormX = normalX[k];
        real NormY = normalY[k];
        real NormZ = normalZ[k];
        ////////////////////////////////////////////////////////////////////////////////
        real* q_dirE, * q_dirW, * q_dirN, * q_dirS, * q_dirT, * q_dirB,
            * q_dirNE, * q_dirSW, * q_dirSE, * q_dirNW, * q_dirTE, * q_dirBW,
            * q_dirBE, * q_dirTW, * q_dirTN, * q_dirBS, * q_dirBN, * q_dirTS,
            * q_dirTNE, * q_dirTSW, * q_dirTSE, * q_dirTNW, * q_dirBNE, * q_dirBSW,
            * q_dirBSE, * q_dirBNW;
        q_dirE   = &Qarrays[dP00 * numberOfBCnodes];
        q_dirW   = &Qarrays[dM00 * numberOfBCnodes];
        q_dirN   = &Qarrays[d0P0 * numberOfBCnodes];
        q_dirS   = &Qarrays[d0M0 * numberOfBCnodes];
        q_dirT   = &Qarrays[d00P * numberOfBCnodes];
        q_dirB   = &Qarrays[d00M * numberOfBCnodes];
        q_dirNE  = &Qarrays[dPP0 * numberOfBCnodes];
        q_dirSW  = &Qarrays[dMM0 * numberOfBCnodes];
        q_dirSE  = &Qarrays[dPM0 * numberOfBCnodes];
        q_dirNW  = &Qarrays[dMP0 * numberOfBCnodes];
        q_dirTE  = &Qarrays[dP0P * numberOfBCnodes];
        q_dirBW  = &Qarrays[dM0M * numberOfBCnodes];
        q_dirBE  = &Qarrays[dP0M * numberOfBCnodes];
        q_dirTW  = &Qarrays[dM0P * numberOfBCnodes];
        q_dirTN  = &Qarrays[d0PP * numberOfBCnodes];
        q_dirBS  = &Qarrays[d0MM * numberOfBCnodes];
        q_dirBN  = &Qarrays[d0PM * numberOfBCnodes];
        q_dirTS  = &Qarrays[d0MP * numberOfBCnodes];
        q_dirTNE = &Qarrays[dPPP * numberOfBCnodes];
        q_dirTSW = &Qarrays[dMMP * numberOfBCnodes];
        q_dirTSE = &Qarrays[dPMP * numberOfBCnodes];
        q_dirTNW = &Qarrays[dMPP * numberOfBCnodes];
        q_dirBNE = &Qarrays[dPPM * numberOfBCnodes];
        q_dirBSW = &Qarrays[dMMM * numberOfBCnodes];
        q_dirBSE = &Qarrays[dPMM * numberOfBCnodes];
        q_dirBNW = &Qarrays[dMPM * numberOfBCnodes];
        ////////////////////////////////////////////////////////////////////////////////
        //index
        unsigned int KQK   = QindexArray[k];
        unsigned int kzero = KQK;
        unsigned int ke    = KQK;
        unsigned int kw    = neighborX[KQK];
        unsigned int kn    = KQK;
        unsigned int ks    = neighborY[KQK];
        unsigned int kt    = KQK;
        unsigned int kb    = neighborZ[KQK];
        unsigned int ksw   = neighborY[kw];
        unsigned int kne   = KQK;
        unsigned int kse   = ks;
        unsigned int knw   = kw;
        unsigned int kbw   = neighborZ[kw];
        unsigned int kte   = KQK;
        unsigned int kbe   = kb;
        unsigned int ktw   = kw;
        unsigned int kbs   = neighborZ[ks];
        unsigned int ktn   = KQK;
        unsigned int kbn   = kb;
        unsigned int kts   = ks;
        unsigned int ktse  = ks;
        unsigned int kbnw  = kbw;
        unsigned int ktnw  = kw;
        unsigned int kbse  = kbs;
        unsigned int ktsw  = ksw;
        unsigned int kbne  = kb;
        unsigned int ktne  = KQK;
        unsigned int kbsw  = neighborZ[ksw];
        ////////////////////////////////////////////////////////////////////////////////
        real f_E, f_W, f_N, f_S, f_T, f_B, f_NE, f_SW, f_SE, f_NW, f_TE, f_BW, f_BE,
            f_TW, f_TN, f_BS, f_BN, f_TS, f_TNE, f_TSW, f_TSE, f_TNW, f_BNE, f_BSW, f_BSE, f_BNW;

        f_W   = (D.f[dP00])[ke];
        f_E   = (D.f[dM00])[kw];
        f_S   = (D.f[d0P0])[kn];
        f_N   = (D.f[d0M0])[ks];
        f_B   = (D.f[d00P])[kt];
        f_T   = (D.f[d00M])[kb];
        f_SW  = (D.f[dPP0])[kne];
        f_NE  = (D.f[dMM0])[ksw];
        f_NW  = (D.f[dPM0])[kse];
        f_SE  = (D.f[dMP0])[knw];
        f_BW  = (D.f[dP0P])[kte];
        f_TE  = (D.f[dM0M])[kbw];
        f_TW  = (D.f[dP0M])[kbe];
        f_BE  = (D.f[dM0P])[ktw];
        f_BS  = (D.f[d0PP])[ktn];
        f_TN  = (D.f[d0MM])[kbs];
        f_TS  = (D.f[d0PM])[kbn];
        f_BN  = (D.f[d0MP])[kts];
        f_BSW = (D.f[dPPP])[ktne];
        f_BNE = (D.f[dMMP])[ktsw];
        f_BNW = (D.f[dPMP])[ktse];
        f_BSE = (D.f[dMPP])[ktnw];
        f_TSW = (D.f[dPPM])[kbne];
        f_TNE = (D.f[dMMM])[kbsw];
        f_TNW = (D.f[dPMM])[kbse];
        f_TSE = (D.f[dMPM])[kbnw];
        ////////////////////////////////////////////////////////////////////////////////
        real vx1, vx2, vx3, drho, q;
        drho = f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]);

        vx1 = (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) +
            (f_E - f_W)) / (c1o1 + drho);


        vx2 = ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) +
            (f_N - f_S)) / (c1o1 + drho);

        vx3 = (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
            (-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) +
            (f_T - f_B)) / (c1o1 + drho);

        real cu_sq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) * (c1o1 + drho);

        ////////////////////////////////////////////////////////////////////////////////
        f_W   = (DAD.f[dP00])[ke];
        f_E   = (DAD.f[dM00])[kw];
        f_S   = (DAD.f[d0P0])[kn];
        f_N   = (DAD.f[d0M0])[ks];
        f_B   = (DAD.f[d00P])[kt];
        f_T   = (DAD.f[d00M])[kb];
        f_SW  = (DAD.f[dPP0])[kne];
        f_NE  = (DAD.f[dMM0])[ksw];
        f_NW  = (DAD.f[dPM0])[kse];
        f_SE  = (DAD.f[dMP0])[knw];
        f_BW  = (DAD.f[dP0P])[kte];
        f_TE  = (DAD.f[dM0M])[kbw];
        f_TW  = (DAD.f[dP0M])[kbe];
        f_BE  = (DAD.f[dM0P])[ktw];
        f_BS  = (DAD.f[d0PP])[ktn];
        f_TN  = (DAD.f[d0MM])[kbs];
        f_TS  = (DAD.f[d0PM])[kbn];
        f_BN  = (DAD.f[d0MP])[kts];
        f_BSW = (DAD.f[dPPP])[ktne];
        f_BNE = (DAD.f[dMMP])[ktsw];
        f_BNW = (DAD.f[dPMP])[ktse];
        f_BSE = (DAD.f[dMPP])[ktnw];
        f_TSW = (DAD.f[dPPM])[kbne];
        f_TNE = (DAD.f[dMMM])[kbsw];
        f_TNW = (DAD.f[dPMM])[kbse];
        f_TSE = (DAD.f[dMPM])[kbnw];
        //////////////////////////////////////////////////////////////////////////
        if (!isEvenTimestep)
        {
            DAD.f[dP00] = &distributionsAD[dP00 * numberOfLBnodes];
            DAD.f[dM00] = &distributionsAD[dM00 * numberOfLBnodes];
            DAD.f[d0P0] = &distributionsAD[d0P0 * numberOfLBnodes];
            DAD.f[d0M0] = &distributionsAD[d0M0 * numberOfLBnodes];
            DAD.f[d00P] = &distributionsAD[d00P * numberOfLBnodes];
            DAD.f[d00M] = &distributionsAD[d00M * numberOfLBnodes];
            DAD.f[dPP0] = &distributionsAD[dPP0 * numberOfLBnodes];
            DAD.f[dMM0] = &distributionsAD[dMM0 * numberOfLBnodes];
            DAD.f[dPM0] = &distributionsAD[dPM0 * numberOfLBnodes];
            DAD.f[dMP0] = &distributionsAD[dMP0 * numberOfLBnodes];
            DAD.f[dP0P] = &distributionsAD[dP0P * numberOfLBnodes];
            DAD.f[dM0M] = &distributionsAD[dM0M * numberOfLBnodes];
            DAD.f[dP0M] = &distributionsAD[dP0M * numberOfLBnodes];
            DAD.f[dM0P] = &distributionsAD[dM0P * numberOfLBnodes];
            DAD.f[d0PP] = &distributionsAD[d0PP * numberOfLBnodes];
            DAD.f[d0MM] = &distributionsAD[d0MM * numberOfLBnodes];
            DAD.f[d0PM] = &distributionsAD[d0PM * numberOfLBnodes];
            DAD.f[d0MP] = &distributionsAD[d0MP * numberOfLBnodes];
            DAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
            DAD.f[dPPP] = &distributionsAD[dPPP * numberOfLBnodes];
            DAD.f[dMMP] = &distributionsAD[dMMP * numberOfLBnodes];
            DAD.f[dPMP] = &distributionsAD[dPMP * numberOfLBnodes];
            DAD.f[dMPP] = &distributionsAD[dMPP * numberOfLBnodes];
            DAD.f[dPPM] = &distributionsAD[dPPM * numberOfLBnodes];
            DAD.f[dMMM] = &distributionsAD[dMMM * numberOfLBnodes];
            DAD.f[dPMM] = &distributionsAD[dPMM * numberOfLBnodes];
            DAD.f[dMPM] = &distributionsAD[dMPM * numberOfLBnodes];
        }
        else
        {
            DAD.f[dM00] = &distributionsAD[dP00 * numberOfLBnodes];
            DAD.f[dP00] = &distributionsAD[dM00 * numberOfLBnodes];
            DAD.f[d0M0] = &distributionsAD[d0P0 * numberOfLBnodes];
            DAD.f[d0P0] = &distributionsAD[d0M0 * numberOfLBnodes];
            DAD.f[d00M] = &distributionsAD[d00P * numberOfLBnodes];
            DAD.f[d00P] = &distributionsAD[d00M * numberOfLBnodes];
            DAD.f[dMM0] = &distributionsAD[dPP0 * numberOfLBnodes];
            DAD.f[dPP0] = &distributionsAD[dMM0 * numberOfLBnodes];
            DAD.f[dMP0] = &distributionsAD[dPM0 * numberOfLBnodes];
            DAD.f[dPM0] = &distributionsAD[dMP0 * numberOfLBnodes];
            DAD.f[dM0M] = &distributionsAD[dP0P * numberOfLBnodes];
            DAD.f[dP0P] = &distributionsAD[dM0M * numberOfLBnodes];
            DAD.f[dM0P] = &distributionsAD[dP0M * numberOfLBnodes];
            DAD.f[dP0M] = &distributionsAD[dM0P * numberOfLBnodes];
            DAD.f[d0MM] = &distributionsAD[d0PP * numberOfLBnodes];
            DAD.f[d0PP] = &distributionsAD[d0MM * numberOfLBnodes];
            DAD.f[d0MP] = &distributionsAD[d0PM * numberOfLBnodes];
            DAD.f[d0PM] = &distributionsAD[d0MP * numberOfLBnodes];
            DAD.f[d000] = &distributionsAD[d000 * numberOfLBnodes];
            DAD.f[dPPP] = &distributionsAD[dMMM * numberOfLBnodes];
            DAD.f[dMMP] = &distributionsAD[dPPM * numberOfLBnodes];
            DAD.f[dPMP] = &distributionsAD[dMPM * numberOfLBnodes];
            DAD.f[dMPP] = &distributionsAD[dPMM * numberOfLBnodes];
            DAD.f[dPPM] = &distributionsAD[dMMP * numberOfLBnodes];
            DAD.f[dMMM] = &distributionsAD[dPPP * numberOfLBnodes];
            DAD.f[dPMM] = &distributionsAD[dMPP * numberOfLBnodes];
            DAD.f[dMPM] = &distributionsAD[dPMP * numberOfLBnodes];
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        real concentration =
            f_TSE + f_TNW + f_TNE + f_TSW + f_BSE + f_BNW + f_BNE + f_BSW +
            f_BN + f_TS + f_TN + f_BS + f_BE + f_TW + f_TE + f_BW + f_SE + f_NW + f_NE + f_SW +
            f_T + f_B + f_N + f_S + f_E + f_W + ((D.f[d000])[kzero]);

        real jx1 =
            (((f_TSE - f_BNW) - (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BE - f_TW) + (f_TE - f_BW)) + ((f_SE - f_NW) + (f_NE - f_SW)) +
            (f_E - f_W)) - (vx1 * concentration);

        real jx2 =
            ((-(f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) - (f_TSW - f_BNE)) +
            ((f_BN - f_TS) + (f_TN - f_BS)) + (-(f_SE - f_NW) + (f_NE - f_SW)) +
            (f_N - f_S)) - (vx2 * concentration);

        real jx3 =
            (((f_TSE - f_BNW) + (f_TNW - f_BSE)) + ((f_TNE - f_BSW) + (f_TSW - f_BNE)) +
            (-(f_BN - f_TS) + (f_TN - f_BS)) + ((f_TE - f_BW) - (f_BE - f_TW)) +
            (f_T - f_B)) - (vx3 * concentration);

        real NormJ = jx1 * NormX + jx2 * NormY + jx3 * NormZ;

        real jTan1 = jx1 - NormJ * NormX;
        real jTan2 = jx2 - NormJ * NormY;
        real jTan3 = jx3 - NormJ * NormZ;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        q = q_dirE[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dM00])[kw  ] = calcDistributionBC_AD(q, c2o27,   vx1,         cu_sq, f_E,   f_W,   omegaDiffusivity,        jTan1,       concentration); }
        q = q_dirW[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[dP00])[ke  ] = calcDistributionBC_AD(q, c2o27,  -vx1,         cu_sq, f_W,   f_E,   omegaDiffusivity,       -jTan1,       concentration); }
        q = q_dirN[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[d0M0])[ks  ] = calcDistributionBC_AD(q, c2o27,   vx2,         cu_sq, f_N,   f_S,   omegaDiffusivity,        jTan2,       concentration); }
        q = q_dirS[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[d0P0])[kn  ] = calcDistributionBC_AD(q, c2o27,  -vx2,         cu_sq, f_S,   f_N,   omegaDiffusivity,       -jTan2,       concentration); }
        q = q_dirT[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[d00M])[kb  ] = calcDistributionBC_AD(q, c2o27,   vx3,         cu_sq, f_T,   f_B,   omegaDiffusivity,        jTan3,       concentration); }
        q = q_dirB[k];   if (q >= c0o1 && q <= c1o1) { (DAD.f[d00P])[kt  ] = calcDistributionBC_AD(q, c2o27,  -vx3,         cu_sq, f_B,   f_T,   omegaDiffusivity,       -jTan3,       concentration); }
        q = q_dirNE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dMM0])[ksw ] = calcDistributionBC_AD(q, c1o54,   vx1+vx2,     cu_sq, f_NE,  f_SW,  omegaDiffusivity,  jTan1+jTan2,       concentration); }
        q = q_dirSW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dPP0])[kne ] = calcDistributionBC_AD(q, c1o54,  -vx1-vx2,     cu_sq, f_SW,  f_NE,  omegaDiffusivity, -jTan1-jTan2,       concentration); }
        q = q_dirSE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dMP0])[knw ] = calcDistributionBC_AD(q, c1o54,   vx1-vx2,     cu_sq, f_SE,  f_NW,  omegaDiffusivity,  jTan1-jTan2,       concentration); }
        q = q_dirNW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dPM0])[kse ] = calcDistributionBC_AD(q, c1o54,  -vx1+vx2,     cu_sq, f_NW,  f_SE,  omegaDiffusivity, -jTan1+jTan2,       concentration); }
        q = q_dirTE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dM0M])[kbw ] = calcDistributionBC_AD(q, c1o54,   vx1    +vx3, cu_sq, f_TE,  f_BW,  omegaDiffusivity,  jTan1      +jTan3, concentration); }
        q = q_dirBW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dP0P])[kte ] = calcDistributionBC_AD(q, c1o54,  -vx1    -vx3, cu_sq, f_BW,  f_TE,  omegaDiffusivity, -jTan1      -jTan3, concentration); }
        q = q_dirBE[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dM0P])[ktw ] = calcDistributionBC_AD(q, c1o54,   vx1    -vx3, cu_sq, f_BE,  f_TW,  omegaDiffusivity,  jTan1      -jTan3, concentration); }
        q = q_dirTW[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[dP0M])[kbe ] = calcDistributionBC_AD(q, c1o54,  -vx1    +vx3, cu_sq, f_TW,  f_BE,  omegaDiffusivity, -jTan1      +jTan3, concentration); }
        q = q_dirTN[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[d0MM])[kbs ] = calcDistributionBC_AD(q, c1o54,       vx2+vx3, cu_sq, f_TN,  f_BS,  omegaDiffusivity,        jTan2+jTan3, concentration); }
        q = q_dirBS[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[d0PP])[ktn ] = calcDistributionBC_AD(q, c1o54,      -vx2-vx3, cu_sq, f_BS,  f_TN,  omegaDiffusivity,       -jTan2-jTan3, concentration); }
        q = q_dirBN[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[d0MP])[kts ] = calcDistributionBC_AD(q, c1o54,       vx2-vx3, cu_sq, f_BN,  f_TS,  omegaDiffusivity,        jTan2-jTan3, concentration); }
        q = q_dirTS[k];  if (q >= c0o1 && q <= c1o1) { (DAD.f[d0PM])[kbn ] = calcDistributionBC_AD(q, c1o54,      -vx2+vx3, cu_sq, f_TS,  f_BN,  omegaDiffusivity,       -jTan2+jTan3, concentration); }
        q = q_dirTNE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dMMM])[kbsw] = calcDistributionBC_AD(q, c1o216,  vx1+vx2+vx3, cu_sq, f_TNE, f_BSW, omegaDiffusivity,  jTan1+jTan2+jTan3, concentration); }
        q = q_dirBSW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dPPP])[ktne] = calcDistributionBC_AD(q, c1o216, -vx1-vx2-vx3, cu_sq, f_BSW, f_TNE, omegaDiffusivity, -jTan1-jTan2-jTan3, concentration); }
        q = q_dirBNE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dMMP])[ktsw] = calcDistributionBC_AD(q, c1o216,  vx1+vx2-vx3, cu_sq, f_BNE, f_TSW, omegaDiffusivity,  jTan1+jTan2-jTan3, concentration); }
        q = q_dirTSW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dPPM])[kbne] = calcDistributionBC_AD(q, c1o216, -vx1-vx2+vx3, cu_sq, f_TSW, f_BNE, omegaDiffusivity, -jTan1-jTan2+jTan3, concentration); }
        q = q_dirTSE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dMPM])[kbnw] = calcDistributionBC_AD(q, c1o216,  vx1-vx2+vx3, cu_sq, f_TSE, f_BNW, omegaDiffusivity,  jTan1-jTan2+jTan3, concentration); }
        q = q_dirBNW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dPMP])[ktse] = calcDistributionBC_AD(q, c1o216, -vx1+vx2-vx3, cu_sq, f_BNW, f_TSE, omegaDiffusivity, -jTan1+jTan2-jTan3, concentration); }
        q = q_dirBSE[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dMPP])[ktnw] = calcDistributionBC_AD(q, c1o216,  vx1-vx2-vx3, cu_sq, f_BSE, f_TNW, omegaDiffusivity,  jTan1-jTan2-jTan3, concentration); }
        q = q_dirTNW[k]; if (q >= c0o1 && q <= c1o1) { (DAD.f[dPMM])[kbse] = calcDistributionBC_AD(q, c1o216, -vx1+vx2+vx3, cu_sq, f_TNW, f_BSE, omegaDiffusivity, -jTan1+jTan2+jTan3, concentration); }
    }
}
