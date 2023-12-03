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
//! \author Henry Korb, Henrik Asmuth, Martin Schoenherr
//======================================================================================
#include "LBM/LB.h"
#include <basics/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>
#include <lbm/MacroscopicQuantities.h>

#include "LBM/GPUHelperFunctions/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void PrecursorEquilibrium_Device(
    int *subgridDistanceIndices,
    int numberOfBCnodes,
    int numberOfPrecursorNodes,
    real omega,
    real* distributions,
    uint* neighborX,
    uint* neighborY,
    uint* neighborZ,
    uint* neighbors0PP,
    uint* neighbors0PM,
    uint* neighbors0MP,
    uint* neighbors0MM,
    real* weights0PP,
    real* weights0PM,
    real* weights0MP,
    real* weights0MM,
    real* vLast,
    real* vCurrent,
    real velocityX,
    real velocityY,
    real velocityZ,
    real timeRatio,
    real velocityRatio,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    if(nodeIndex>=numberOfBCnodes) return;

    ////////////////////////////////////////////////////////////////////////////////
    // interpolation of velocity
    real vxLastInterpd, vyLastInterpd, vzLastInterpd;
    real vxNextInterpd, vyNextInterpd, vzNextInterpd;

    uint kNeighbor0PP = neighbors0PP[nodeIndex];
    real d0PP = weights0PP[nodeIndex];

    real* vxLast = vLast;
    real* vyLast = &vLast[numberOfPrecursorNodes];
    real* vzLast = &vLast[2*numberOfPrecursorNodes];

    real* vxCurrent = vCurrent;
    real* vyCurrent = &vCurrent[numberOfPrecursorNodes];
    real* vzCurrent = &vCurrent[2*numberOfPrecursorNodes];

    if(d0PP < 1e6)
    {
        uint kNeighbor0PM = neighbors0PM[nodeIndex];
        uint kNeighbor0MP = neighbors0MP[nodeIndex];
        uint kNeighbor0MM = neighbors0MM[nodeIndex];

        real d0PM = weights0PM[nodeIndex];
        real d0MP = weights0MP[nodeIndex];
        real d0MM = weights0MM[nodeIndex];

        real invWeightSum = 1.f/(d0PP+d0PM+d0MP+d0MM);

        vxLastInterpd = (vxLast[kNeighbor0PP]*d0PP + vxLast[kNeighbor0PM]*d0PM + vxLast[kNeighbor0MP]*d0MP + vxLast[kNeighbor0MM]*d0MM)*invWeightSum;
        vyLastInterpd = (vyLast[kNeighbor0PP]*d0PP + vyLast[kNeighbor0PM]*d0PM + vyLast[kNeighbor0MP]*d0MP + vyLast[kNeighbor0MM]*d0MM)*invWeightSum;
        vzLastInterpd = (vzLast[kNeighbor0PP]*d0PP + vzLast[kNeighbor0PM]*d0PM + vzLast[kNeighbor0MP]*d0MP + vzLast[kNeighbor0MM]*d0MM)*invWeightSum;

        vxNextInterpd = (vxCurrent[kNeighbor0PP]*d0PP + vxCurrent[kNeighbor0PM]*d0PM + vxCurrent[kNeighbor0MP]*d0MP + vxCurrent[kNeighbor0MM]*d0MM)*invWeightSum;
        vyNextInterpd = (vyCurrent[kNeighbor0PP]*d0PP + vyCurrent[kNeighbor0PM]*d0PM + vyCurrent[kNeighbor0MP]*d0MP + vyCurrent[kNeighbor0MM]*d0MM)*invWeightSum;
        vzNextInterpd = (vzCurrent[kNeighbor0PP]*d0PP + vzCurrent[kNeighbor0PM]*d0PM + vzCurrent[kNeighbor0MP]*d0MP + vzCurrent[kNeighbor0MM]*d0MM)*invWeightSum;
    }
    else
    {
        vxLastInterpd = vxLast[kNeighbor0PP];
        vyLastInterpd = vyLast[kNeighbor0PP];
        vzLastInterpd = vzLast[kNeighbor0PP];

        vxNextInterpd = vxCurrent[kNeighbor0PP];
        vyNextInterpd = vyCurrent[kNeighbor0PP];
        vzNextInterpd = vzCurrent[kNeighbor0PP];
    }

    // if(k==16300) printf("%f %f %f\n", vxLastInterpd, vyLastInterpd, vzLastInterpd);
    real VeloX = (velocityX + (1.f-timeRatio)*vxLastInterpd + timeRatio*vxNextInterpd)/velocityRatio;
    real VeloY = (velocityY + (1.f-timeRatio)*vyLastInterpd + timeRatio*vyNextInterpd)/velocityRatio;
    real VeloZ = (velocityZ + (1.f-timeRatio)*vzLastInterpd + timeRatio*vzNextInterpd)/velocityRatio;
    // From here on just a copy of QVelDeviceCompZeroPress
    ////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep
    //! is based on the esoteric twist algorithm \ref <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier
    //! et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
    //!
    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

    unsigned int KQK  = subgridDistanceIndices[nodeIndex]; //QK
    unsigned int k000 = KQK; //000
    unsigned int kP00 = KQK; //P00
    unsigned int kM00 = neighborX[KQK]; //M00
    unsigned int k0P0   = KQK; //n
    unsigned int k0M0   = neighborY[KQK]; //s
    unsigned int k00P   = KQK; //t
    unsigned int k00M   = neighborZ[KQK]; //b
    unsigned int kMM0  = neighborY[kM00]; //sw
    unsigned int kPP0  = KQK; //ne
    unsigned int kPM0  = k0M0; //se
    unsigned int kMP0  = kM00; //nw
    unsigned int kM0M  = neighborZ[kM00]; //bw
    unsigned int kP0P  = KQK; //te
    unsigned int kP0M  = k00M; //be
    unsigned int k0PP  = KQK; //tn
    unsigned int k0MM  = neighborZ[k0M0]; //bs
    unsigned int kM0P  = kM00; //tw
    unsigned int k0PM  = k00M; //bn
    unsigned int k0MP  = k0M0; //ts
    unsigned int kPMP = k0M0; //tse
    unsigned int kMPM = kM0M; //bnw
    unsigned int kMPP = kM00; //tnw
    unsigned int kPMM = k0MM; //bse
    unsigned int kMMP = kMM0; //tsw
    unsigned int kPPM = k00M; //bne
    unsigned int kPPP = KQK; //tne
    unsigned int kMMM = neighborZ[kMM0]; //bsw

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // based on BGK Plus Comp
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    real f_M00 = (dist.f[dP00])[kP00];
    real f_P00 = (dist.f[dM00])[kM00];
    real f_0M0 = (dist.f[d0P0])[k0P0];
    real f_0P0 = (dist.f[d0M0])[k0M0];
    real f_00M = (dist.f[d00P])[k00P];
    real f_00P = (dist.f[d00M])[k00M];
    real f_MM0 = (dist.f[dPP0])[kPP0];
    real f_PP0 = (dist.f[dMM0])[kMM0];
    real f_MP0 = (dist.f[dPM0])[kPM0];
    real f_PM0 = (dist.f[dMP0])[kMP0];
    real f_M0M = (dist.f[dP0P])[kP0P];
    real f_P0P = (dist.f[dM0M])[kM0M];
    real f_M0P = (dist.f[dP0M])[kP0M];
    real f_P0M = (dist.f[dM0P])[kM0P];
    real f_0MM = (dist.f[vf::lbm::dir::d0PP])[k0PP];
    real f_0PP = (dist.f[d0MM])[k0MM];
    real f_0PM = (dist.f[d0MP])[k0MP];
    real f_0MP = (dist.f[d0PM])[k0PM];
    real f_000 = (dist.f[d000])[k000];
    real f_MMM = (dist.f[dPPP])[kPPP];
    real f_PPM = (dist.f[dMMP])[kMMP];
    real f_MPM = (dist.f[dPMP])[kPMP];
    real f_PMM = (dist.f[dMPP])[kMPP];
    real f_MMP = (dist.f[dPPM])[kPPM];
    real f_PPP = (dist.f[dMMM])[kMMM];
    real f_MPP = (dist.f[dPMM])[kPMM];
    real f_PMP = (dist.f[dMPM])[kMPM];

      ////////////////////////////////////////////////////////////////////////////////
      //! - Set macroscopic quantities
      //!
      real drho = c0o1;

      real vx1  = VeloX;

      real vx2  = VeloY;

      real vx3  = VeloZ;

      real cusq = c3o2 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3);

      ////////////////////////////////////////////////////////////////////////////////
      f_000 = c8o27* (drho-(drho+c1o1)*cusq);
      f_P00 = c2o27* (drho+(drho+c1o1)*(c3o1*( vx1        )+c9o2*( vx1        )*( vx1        )-cusq));
      f_M00 = c2o27* (drho+(drho+c1o1)*(c3o1*(-vx1        )+c9o2*(-vx1        )*(-vx1        )-cusq));
      f_0P0 = c2o27* (drho+(drho+c1o1)*(c3o1*(    vx2     )+c9o2*(     vx2    )*(     vx2    )-cusq));
      f_0M0 = c2o27* (drho+(drho+c1o1)*(c3o1*(   -vx2     )+c9o2*(    -vx2    )*(    -vx2    )-cusq));
      f_00P = c2o27* (drho+(drho+c1o1)*(c3o1*(         vx3)+c9o2*(         vx3)*(         vx3)-cusq));
      f_00M = c2o27* (drho+(drho+c1o1)*(c3o1*(        -vx3)+c9o2*(        -vx3)*(        -vx3)-cusq));
      f_PP0 = c1o54* (drho+(drho+c1o1)*(c3o1*( vx1+vx2    )+c9o2*( vx1+vx2    )*( vx1+vx2    )-cusq));
      f_MM0 = c1o54* (drho+(drho+c1o1)*(c3o1*(-vx1-vx2    )+c9o2*(-vx1-vx2    )*(-vx1-vx2    )-cusq));
      f_PM0 = c1o54* (drho+(drho+c1o1)*(c3o1*( vx1-vx2    )+c9o2*( vx1-vx2    )*( vx1-vx2    )-cusq));
      f_MP0 = c1o54* (drho+(drho+c1o1)*(c3o1*(-vx1+vx2    )+c9o2*(-vx1+vx2    )*(-vx1+vx2    )-cusq));
      f_P0P = c1o54* (drho+(drho+c1o1)*(c3o1*( vx1    +vx3)+c9o2*( vx1    +vx3)*( vx1    +vx3)-cusq));
      f_M0M = c1o54* (drho+(drho+c1o1)*(c3o1*(-vx1    -vx3)+c9o2*(-vx1    -vx3)*(-vx1    -vx3)-cusq));
      f_P0M = c1o54* (drho+(drho+c1o1)*(c3o1*( vx1    -vx3)+c9o2*( vx1    -vx3)*( vx1    -vx3)-cusq));
      f_M0P = c1o54* (drho+(drho+c1o1)*(c3o1*(-vx1    +vx3)+c9o2*(-vx1    +vx3)*(-vx1    +vx3)-cusq));
      f_0PP = c1o54* (drho+(drho+c1o1)*(c3o1*(     vx2+vx3)+c9o2*(     vx2+vx3)*(     vx2+vx3)-cusq));
      f_0MM = c1o54* (drho+(drho+c1o1)*(c3o1*(    -vx2-vx3)+c9o2*(    -vx2-vx3)*(    -vx2-vx3)-cusq));
      f_0PM = c1o54* (drho+(drho+c1o1)*(c3o1*(     vx2-vx3)+c9o2*(     vx2-vx3)*(     vx2-vx3)-cusq));
      f_0MP = c1o54* (drho+(drho+c1o1)*(c3o1*(    -vx2+vx3)+c9o2*(    -vx2+vx3)*(    -vx2+vx3)-cusq));
      f_PPP = c1o216*(drho+(drho+c1o1)*(c3o1*( vx1+vx2+vx3)+c9o2*( vx1+vx2+vx3)*( vx1+vx2+vx3)-cusq));
      f_MMM = c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1-vx2-vx3)+c9o2*(-vx1-vx2-vx3)*(-vx1-vx2-vx3)-cusq));
      f_PPM = c1o216*(drho+(drho+c1o1)*(c3o1*( vx1+vx2-vx3)+c9o2*( vx1+vx2-vx3)*( vx1+vx2-vx3)-cusq));
      f_MMP = c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1-vx2+vx3)+c9o2*(-vx1-vx2+vx3)*(-vx1-vx2+vx3)-cusq));
      f_PMP = c1o216*(drho+(drho+c1o1)*(c3o1*( vx1-vx2+vx3)+c9o2*( vx1-vx2+vx3)*( vx1-vx2+vx3)-cusq));
      f_MPM = c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1+vx2-vx3)+c9o2*(-vx1+vx2-vx3)*(-vx1+vx2-vx3)-cusq));
      f_PMM = c1o216*(drho+(drho+c1o1)*(c3o1*( vx1-vx2-vx3)+c9o2*( vx1-vx2-vx3)*( vx1-vx2-vx3)-cusq));
      f_MPP = c1o216*(drho+(drho+c1o1)*(c3o1*(-vx1+vx2+vx3)+c9o2*(-vx1+vx2+vx3)*(-vx1+vx2+vx3)-cusq));

      ////////////////////////////////////////////////////////////////////////////////
      //! write the new distributions to the bc nodes
      //!
      (dist.f[dP00])[kP00] = f_M00;
      (dist.f[dPP0])[kPP0] = f_MM0;
      (dist.f[dP0M])[kP0M] = f_M0P;
      (dist.f[dPM0])[kPM0] = f_MP0;
      (dist.f[dPMP])[kPMP] = f_MPM;
      (dist.f[dP0P])[kP0P] = f_M0M;
      (dist.f[dPPM])[kPPM] = f_MMP;
      (dist.f[dPPP])[kPPP] = f_MMM;
      (dist.f[dPMM])[kPMM] = f_MPP;

      (dist.f[dM00])[kM00] = f_P00;
      (dist.f[dMM0])[kMM0] = f_PP0;
      (dist.f[dM0M])[kM0M] = f_P0P;
      (dist.f[dMP0])[kMP0] = f_PM0;
      (dist.f[dM0P])[kM0P] = f_P0M;
      (dist.f[dMMM])[kMMM] = f_PPP;
      (dist.f[dMMP])[kMMP] = f_PPM;
      (dist.f[dMPP])[kMPP] = f_PMM;
      (dist.f[dMPM])[kMPM] = f_PMP;

      (dist.f[d0P0])[k0P0] = f_0M0;
      (dist.f[d0M0])[k0M0] = f_0P0;
      (dist.f[d00P])[k00P] = f_00M;
      (dist.f[d00M])[k00M] = f_00P;
      (dist.f[vf::lbm::dir::d0PP])[k0PP] = f_0MM;
      (dist.f[d0MM])[k0MM] = f_0PP;
      (dist.f[d0PM])[k0PM] = f_0MP;
      (dist.f[d0MP])[k0MP] = f_0PM;
      (dist.f[d000])[k000] = f_000;
}

