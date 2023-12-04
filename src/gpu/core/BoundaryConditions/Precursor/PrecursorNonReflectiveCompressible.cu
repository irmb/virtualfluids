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

__global__ void PrecursorNonReflectiveCompressible_Device(
    int* subgridDistanceIndices,
    int numberOfBCnodes,
    int numberOfPrecursorNodes,
    int sizeQ,
    real omega,
    real* distributions,
    real* subgridDistances,
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

    // if(k==16300)s printf("%f %f %f\n", vxLastInterpd, vyLastInterpd, vzLastInterpd);
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
    getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);

    unsigned int KQK  = subgridDistanceIndices[nodeIndex];
    unsigned int k000= KQK;
    unsigned int kP00   = KQK;
    unsigned int kM00   = neighborX[KQK];
    unsigned int k0P0   = KQK;
    unsigned int k0M0   = neighborY[KQK];
    unsigned int k00P   = KQK;
    unsigned int k00M   = neighborZ[KQK];
    unsigned int kMM0  = neighborY[kM00];
    unsigned int kPP0  = KQK;
    unsigned int kPM0  = k0M0;
    unsigned int kMP0  = kM00;
    unsigned int kM0M  = neighborZ[kM00];
    unsigned int kP0P  = KQK;
    unsigned int kP0M  = k00M;
    unsigned int kM0P  = kM00;
    unsigned int k0PP  = KQK;
    unsigned int k0MM  = neighborZ[k0M0];
    unsigned int k0PM  = k00M;
    unsigned int k0MP  = k0M0;
    unsigned int kPMP = k0M0;
    unsigned int kMPM = kM0M;
    unsigned int kMPP = kM00;
    unsigned int kPMM = k0MM;
    unsigned int kMMP = kMM0;
    unsigned int kPPM = k00M;
    unsigned int kPPP = KQK;
    unsigned int kMMM = neighborZ[kMM0];

    ////////////////////////////////////////////////////////////////////////////////
    //! - Set local distributions
    //!
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
    real f_0MP = (dist.f[d0PM])[k0PM];
    real f_0PM = (dist.f[d0MP])[k0MP];
    real f_MMM = (dist.f[dPPP])[kPPP];
    real f_PPM = (dist.f[dMMP])[kMMP];
    real f_MPM = (dist.f[dPMP])[kPMP];
    real f_PMM = (dist.f[dMPP])[kMPP];
    real f_MMP = (dist.f[dPPM])[kPPM];
    real f_PPP = (dist.f[dMMM])[kMMM];
    real f_MPP = (dist.f[dPMM])[kPMM];
    real f_PMP = (dist.f[dMPM])[kMPM];

    SubgridDistances27 subgridD;
    getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);

    ////////////////////////////////////////////////////////////////////////////////
      real drho   =  f_PMP + f_MPP + f_PPP + f_MMP + f_PMM + f_MPM + f_PPM + f_MMM +
                     f_0PM + f_0PP + f_0MP + f_0MM + f_P0M + f_M0P + f_P0P + f_M0M + f_PM0 + f_MP0 + f_PP0 + f_MM0 +
                     f_00P + f_00M + f_0P0 + f_0M0 + f_P00 + f_M00 + ((dist.f[d000])[k000]);

      real vx1 =  (((f_PMP - f_MPM) - (f_MPP - f_PMM)) + ((f_PPP - f_MMM) - (f_MMP - f_PPM)) +
                      ((f_P0M - f_M0P)   + (f_P0P - f_M0M))   + ((f_PM0 - f_MP0)   + (f_PP0 - f_MM0)) +
                      (f_P00 - f_M00)) / (c1o1 + drho);


      real vx2 =   ((-(f_PMP - f_MPM) + (f_MPP - f_PMM)) + ((f_PPP - f_MMM) - (f_MMP - f_PPM)) +
                       ((f_0PM - f_0MP)   + (f_0PP - f_0MM))    + (-(f_PM0 - f_MP0)  + (f_PP0 - f_MM0)) +
                       (f_0P0 - f_0M0)) / (c1o1 + drho);

      real vx3 =   (((f_PMP - f_MPM) + (f_MPP - f_PMM)) + ((f_PPP - f_MMM) + (f_MMP - f_PPM)) +
                       (-(f_0PM - f_0MP)  + (f_0PP - f_0MM))   + ((f_P0P - f_M0M)   - (f_P0M - f_M0P)) +
                       (f_00P - f_00M)) / (c1o1 + drho);


    // if(k==16383 || k==0) printf("k %d kQ %d drho = %f u %f v %f w %f\n",k, KQK, drho, vx1, vx2, vx3);
      real cu_sq=c3o2*(vx1*vx1+vx2*vx2+vx3*vx3) * (c1o1 + drho);
    //////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    //! - Update distributions with subgrid distance (q) between zero and one
    real feq, q, velocityLB, velocityBC;
    q = (subgridD.q[dP00])[nodeIndex];
    if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
    {
        velocityLB = vx1;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloX;
        (dist.f[dM00])[kM00] = getInterpolatedDistributionForVeloWithPressureBC(q, f_P00, f_M00, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[dM00])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloX;
        (dist.f[dP00])[kP00] = getInterpolatedDistributionForVeloWithPressureBC(q, f_M00, f_P00, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[d0P0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloY;
        (dist.f[d0M0])[d0M0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0P0, f_0M0, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[d0M0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloY;
        (dist.f[d0P0])[k0P0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0M0, f_0P0, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[d00P])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloZ;
        (dist.f[d00M])[k00M] = getInterpolatedDistributionForVeloWithPressureBC(q, f_00P, f_00M, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[d00M])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloZ;
        (dist.f[d00P])[k00P] = getInterpolatedDistributionForVeloWithPressureBC(q, f_00M, f_00P, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[dPP0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX + VeloY;
        (dist.f[dMM0])[kMM0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PP0, f_MM0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dMM0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloY;
        (dist.f[dPP0])[kPP0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MM0, f_PP0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dPM0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX - VeloY;
        (dist.f[dMP0])[kMP0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PM0, f_MP0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dMP0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX + VeloY;
        (dist.f[dPM0])[kPM0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MP0, f_PM0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dP0P])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX + VeloZ;
        (dist.f[dM0M])[kM0M] = getInterpolatedDistributionForVeloWithPressureBC(q, f_P0P, f_M0M, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dM0M])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloZ;
        (dist.f[dP0P])[kP0P] = getInterpolatedDistributionForVeloWithPressureBC(q, f_M0M, f_P0P, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dP0M])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX - VeloZ;
        (dist.f[dM0P])[kM0P] = getInterpolatedDistributionForVeloWithPressureBC(q, f_P0M, f_M0P, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dM0P])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX + VeloZ;
        (dist.f[dP0M])[kP0M] = getInterpolatedDistributionForVeloWithPressureBC(q, f_M0P, f_P0M, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[vf::lbm::dir::d0PP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloY + VeloZ;
        (dist.f[d0MM])[k0MM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0PP, f_0MM, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[d0MM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloY - VeloZ;
        (dist.f[vf::lbm::dir::d0PP])[k0PP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0MM, f_0PP, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[d0PM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloY - VeloZ;
        (dist.f[d0MP])[k0MP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0PM, f_0PP, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[d0MP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloY + VeloZ;
        (dist.f[d0PM])[k0PM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0PP, f_0PM, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[dPPP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX + VeloY + VeloZ;
        (dist.f[dMMM])[kMMM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PPP, f_MMM, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[dMMM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX - VeloY - VeloZ;
        (dist.f[dPPP])[kPPP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MMM, f_PPP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[dPPM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX + VeloY - VeloZ;
        (dist.f[dMMP])[kMMP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PPM, f_MMP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[dMMP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX - VeloY + VeloZ;
        (dist.f[dPPM])[kPPM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MMP, f_PPM, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[dPMP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX - VeloY + VeloZ;
        (dist.f[dMPM])[kMPM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PMP, f_MPM, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[dMPM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX + VeloY - VeloZ;
        (dist.f[dPMP])[kPMP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MPM, f_PMP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[dPMM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX - VeloY - VeloZ;
        (dist.f[dMPP])[kMPP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PMM, f_MPP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[dMPP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX + VeloY + VeloZ;
        (dist.f[dPMM])[kPMM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MPP, f_PMM, feq, omega, drho, velocityBC, c1o216);
    }
}
