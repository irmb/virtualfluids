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
//! \file PrecursorBCs27.cu
//! \ingroup GPU
//! \author Henry Korb, Henrik Asmuth
//======================================================================================
#include "LBM/LB.h"
#include <lbm/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>
#include <lbm/MacroscopicQuantities.h>

#include "LBM/GPUHelperFunctions/KernelUtilities.h"

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void QPrecursorDeviceCompZeroPress(
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
    real f_M00 = (dist.f[DIR_P00])[kP00];
    real f_P00 = (dist.f[DIR_M00])[kM00];
    real f_0M0 = (dist.f[DIR_0P0])[k0P0];
    real f_0P0 = (dist.f[DIR_0M0])[k0M0];
    real f_00M = (dist.f[DIR_00P])[k00P];
    real f_00P = (dist.f[DIR_00M])[k00M];
    real f_MM0 = (dist.f[DIR_PP0])[kPP0];
    real f_PP0 = (dist.f[DIR_MM0])[kMM0];
    real f_MP0 = (dist.f[DIR_PM0])[kPM0];
    real f_PM0 = (dist.f[DIR_MP0])[kMP0];
    real f_M0M = (dist.f[DIR_P0P])[kP0P];
    real f_P0P = (dist.f[DIR_M0M])[kM0M];
    real f_M0P = (dist.f[DIR_P0M])[kP0M];
    real f_P0M = (dist.f[DIR_M0P])[kM0P];
    real f_0MM = (dist.f[DIR_0PP])[k0PP];
    real f_0PP = (dist.f[DIR_0MM])[k0MM];
    real f_0MP = (dist.f[DIR_0PM])[k0PM];
    real f_0PM = (dist.f[DIR_0MP])[k0MP];
    real f_MMM = (dist.f[DIR_PPP])[kPPP];
    real f_PPM = (dist.f[DIR_MMP])[kMMP];
    real f_MPM = (dist.f[DIR_PMP])[kPMP];
    real f_PMM = (dist.f[DIR_MPP])[kMPP];
    real f_MMP = (dist.f[DIR_PPM])[kPPM];
    real f_PPP = (dist.f[DIR_MMM])[kMMM];
    real f_MPP = (dist.f[DIR_PMM])[kPMM];
    real f_PMP = (dist.f[DIR_MPM])[kMPM];

    SubgridDistances27 subgridD;
    getPointersToSubgridDistances(subgridD, subgridDistances, numberOfBCnodes);

    ////////////////////////////////////////////////////////////////////////////////
      real drho   =  f_PMP + f_MPP + f_PPP + f_MMP + f_PMM + f_MPM + f_PPM + f_MMM +
                     f_0PM + f_0PP + f_0MP + f_0MM + f_P0M + f_M0P + f_P0P + f_M0M + f_PM0 + f_MP0 + f_PP0 + f_MM0 +
                     f_00P + f_00M + f_0P0 + f_0M0 + f_P00 + f_M00 + ((dist.f[DIR_000])[k000]);

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
    q = (subgridD.q[DIR_P00])[nodeIndex];
    if (q>=c0o1 && q<=c1o1) // only update distribution for q between zero and one
    {
        velocityLB = vx1;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloX;
        (dist.f[DIR_M00])[kM00] = getInterpolatedDistributionForVeloWithPressureBC(q, f_P00, f_M00, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_M00])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloX;
        (dist.f[DIR_P00])[kP00] = getInterpolatedDistributionForVeloWithPressureBC(q, f_M00, f_P00, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_0P0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloY;
        (dist.f[DIR_0M0])[DIR_0M0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0P0, f_0M0, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_0M0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloY;
        (dist.f[DIR_0P0])[k0P0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0M0, f_0P0, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_00P])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = VeloZ;
        (dist.f[DIR_00M])[k00M] = getInterpolatedDistributionForVeloWithPressureBC(q, f_00P, f_00M, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_00M])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c2o27);
        velocityBC = -VeloZ;
        (dist.f[DIR_00P])[k00P] = getInterpolatedDistributionForVeloWithPressureBC(q, f_00M, f_00P, feq, omega, drho, velocityBC, c2o27);
    }

    q = (subgridD.q[DIR_PP0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX + VeloY;
        (dist.f[DIR_MM0])[kMM0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PP0, f_MM0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_MM0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloY;
        (dist.f[DIR_PP0])[kPP0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MM0, f_PP0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_PM0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX - VeloY;
        (dist.f[DIR_MP0])[kMP0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PM0, f_MP0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_MP0])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX + VeloY;
        (dist.f[DIR_PM0])[kPM0] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MP0, f_PM0, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_P0P])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX + VeloZ;
        (dist.f[DIR_M0M])[kM0M] = getInterpolatedDistributionForVeloWithPressureBC(q, f_P0P, f_M0M, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_M0M])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX - VeloZ;
        (dist.f[DIR_P0P])[kP0P] = getInterpolatedDistributionForVeloWithPressureBC(q, f_M0M, f_P0P, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_P0M])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloX - VeloZ;
        (dist.f[DIR_M0P])[kM0P] = getInterpolatedDistributionForVeloWithPressureBC(q, f_P0M, f_M0P, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_M0P])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloX + VeloZ;
        (dist.f[DIR_P0M])[kP0M] = getInterpolatedDistributionForVeloWithPressureBC(q, f_M0P, f_P0M, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0PP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloY + VeloZ;
        (dist.f[DIR_0MM])[k0MM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0PP, f_0MM, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0MM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloY - VeloZ;
        (dist.f[DIR_0PP])[k0PP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0MM, f_0PP, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0PM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = VeloY - VeloZ;
        (dist.f[DIR_0MP])[k0MP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0PM, f_0PP, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_0MP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o54);
        velocityBC = -VeloY + VeloZ;
        (dist.f[DIR_0PM])[k0PM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_0PP, f_0PM, feq, omega, drho, velocityBC, c1o54);
    }

    q = (subgridD.q[DIR_PPP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX + VeloY + VeloZ;
        (dist.f[DIR_MMM])[kMMM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PPP, f_MMM, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MMM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX - VeloY - VeloZ;
        (dist.f[DIR_PPP])[kPPP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MMM, f_PPP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_PPM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 + vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX + VeloY - VeloZ;
        (dist.f[DIR_MMP])[kMMP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PPM, f_MMP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MMP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 - vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX - VeloY + VeloZ;
        (dist.f[DIR_PPM])[kPPM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MMP, f_PPM, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_PMP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX - VeloY + VeloZ;
        (dist.f[DIR_MPM])[kMPM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PMP, f_MPM, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MPM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX + VeloY - VeloZ;
        (dist.f[DIR_PMP])[kPMP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MPM, f_PMP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_PMM])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = vx1 - vx2 - vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = VeloX - VeloY - VeloZ;
        (dist.f[DIR_MPP])[kMPP] = getInterpolatedDistributionForVeloWithPressureBC(q, f_PMM, f_MPP, feq, omega, drho, velocityBC, c1o216);
    }

    q = (subgridD.q[DIR_MPP])[nodeIndex];
    if (q>=c0o1 && q<=c1o1)
    {
        velocityLB = -vx1 + vx2 + vx3;
        feq = getEquilibriumForBC(drho, velocityLB, cu_sq, c1o216);
        velocityBC = -VeloX + VeloY + VeloZ;
        (dist.f[DIR_PMM])[kPMM] = getInterpolatedDistributionForVeloWithPressureBC(q, f_MPP, f_PMM, feq, omega, drho, velocityBC, c1o216);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////











































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void PrecursorDeviceEQ27(
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
    real f_M00 = (dist.f[DIR_P00])[kP00];
    real f_P00 = (dist.f[DIR_M00])[kM00];
    real f_0M0 = (dist.f[DIR_0P0])[k0P0];
    real f_0P0 = (dist.f[DIR_0M0])[k0M0];
    real f_00M = (dist.f[DIR_00P])[k00P];
    real f_00P = (dist.f[DIR_00M])[k00M];
    real f_MM0 = (dist.f[DIR_PP0])[kPP0];
    real f_PP0 = (dist.f[DIR_MM0])[kMM0];
    real f_MP0 = (dist.f[DIR_PM0])[kPM0];
    real f_PM0 = (dist.f[DIR_MP0])[kMP0];
    real f_M0M = (dist.f[DIR_P0P])[kP0P];
    real f_P0P = (dist.f[DIR_M0M])[kM0M];
    real f_M0P = (dist.f[DIR_P0M])[kP0M];
    real f_P0M = (dist.f[DIR_M0P])[kM0P];
    real f_0MM = (dist.f[DIR_0PP])[k0PP];
    real f_0PP = (dist.f[DIR_0MM])[k0MM];
    real f_0PM = (dist.f[DIR_0MP])[k0MP];
    real f_0MP = (dist.f[DIR_0PM])[k0PM];
    real f_000 = (dist.f[DIR_000])[k000];
    real f_MMM = (dist.f[DIR_PPP])[kPPP];
    real f_PPM = (dist.f[DIR_MMP])[kMMP];
    real f_MPM = (dist.f[DIR_PMP])[kPMP];
    real f_PMM = (dist.f[DIR_MPP])[kMPP];
    real f_MMP = (dist.f[DIR_PPM])[kPPM];
    real f_PPP = (dist.f[DIR_MMM])[kMMM];
    real f_MPP = (dist.f[DIR_PMM])[kPMM];
    real f_PMP = (dist.f[DIR_MPM])[kMPM];

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
      (dist.f[DIR_P00])[kP00] = f_M00;
      (dist.f[DIR_PP0])[kPP0] = f_MM0;
      (dist.f[DIR_P0M])[kP0M] = f_M0P;
      (dist.f[DIR_PM0])[kPM0] = f_MP0;
      (dist.f[DIR_PMP])[kPMP] = f_MPM;
      (dist.f[DIR_P0P])[kP0P] = f_M0M;
      (dist.f[DIR_PPM])[kPPM] = f_MMP;
      (dist.f[DIR_PPP])[kPPP] = f_MMM;
      (dist.f[DIR_PMM])[kPMM] = f_MPP;

      (dist.f[DIR_M00])[kM00] = f_P00;
      (dist.f[DIR_MM0])[kMM0] = f_PP0;
      (dist.f[DIR_M0M])[kM0M] = f_P0P;
      (dist.f[DIR_MP0])[kMP0] = f_PM0;
      (dist.f[DIR_M0P])[kM0P] = f_P0M;
      (dist.f[DIR_MMM])[kMMM] = f_PPP;
      (dist.f[DIR_MMP])[kMMP] = f_PPM;
      (dist.f[DIR_MPP])[kMPP] = f_PMM;
      (dist.f[DIR_MPM])[kMPM] = f_PMP;

      (dist.f[DIR_0P0])[k0P0] = f_0M0;
      (dist.f[DIR_0M0])[k0M0] = f_0P0;
      (dist.f[DIR_00P])[k00P] = f_00M;
      (dist.f[DIR_00M])[k00M] = f_00P;
      (dist.f[DIR_0PP])[k0PP] = f_0MM;
      (dist.f[DIR_0MM])[k0MM] = f_0PP;
      (dist.f[DIR_0PM])[k0PM] = f_0MP;
      (dist.f[DIR_0MP])[k0MP] = f_0PM;
      (dist.f[DIR_000])[k000] = f_000;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void PrecursorDeviceDistributions(
    int *subgridDistanceIndices,
    int numberOfBCnodes,
    int numberOfPrecursorNodes,
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
    real* fsLast,
    real* fsNext,
    real timeRatio,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    if(nodeIndex>=numberOfBCnodes) return;

    uint kNeighbor0PP = neighbors0PP[nodeIndex];
    real d0PP = weights0PP[nodeIndex];

    real f0LastInterp, f1LastInterp, f2LastInterp, f3LastInterp, f4LastInterp, f5LastInterp, f6LastInterp, f7LastInterp, f8LastInterp;
    real f0NextInterp, f1NextInterp, f2NextInterp, f3NextInterp, f4NextInterp, f5NextInterp, f6NextInterp, f7NextInterp, f8NextInterp;

    real* f0Last = fsLast;
    real* f1Last = &fsLast[  numberOfPrecursorNodes];
    real* f2Last = &fsLast[2*numberOfPrecursorNodes];
    real* f3Last = &fsLast[3*numberOfPrecursorNodes];
    real* f4Last = &fsLast[4*numberOfPrecursorNodes];
    real* f5Last = &fsLast[5*numberOfPrecursorNodes];
    real* f6Last = &fsLast[6*numberOfPrecursorNodes];
    real* f7Last = &fsLast[7*numberOfPrecursorNodes];
    real* f8Last = &fsLast[8*numberOfPrecursorNodes];

    real* f0Next = fsNext;
    real* f1Next = &fsNext[  numberOfPrecursorNodes];
    real* f2Next = &fsNext[2*numberOfPrecursorNodes];
    real* f3Next = &fsNext[3*numberOfPrecursorNodes];
    real* f4Next = &fsNext[4*numberOfPrecursorNodes];
    real* f5Next = &fsNext[5*numberOfPrecursorNodes];
    real* f6Next = &fsNext[6*numberOfPrecursorNodes];
    real* f7Next = &fsNext[7*numberOfPrecursorNodes];
    real* f8Next = &fsNext[8*numberOfPrecursorNodes];


    if(d0PP<1e6)
    {
        uint kNeighbor0PM = neighbors0PM[nodeIndex];
        uint kNeighbor0MP = neighbors0MP[nodeIndex];
        uint kNeighbor0MM = neighbors0MM[nodeIndex];

        real d0PM = weights0PM[nodeIndex];
        real d0MP = weights0MP[nodeIndex];
        real d0MM = weights0MM[nodeIndex];

        real invWeightSum = 1.f/(d0PP+d0PM+d0MP+d0MM);

        f0LastInterp = (f0Last[kNeighbor0PP]*d0PP + f0Last[kNeighbor0PM]*d0PM + f0Last[kNeighbor0MP]*d0MP + f0Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f0NextInterp = (f0Next[kNeighbor0PP]*d0PP + f0Next[kNeighbor0PM]*d0PM + f0Next[kNeighbor0MP]*d0MP + f0Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f1LastInterp = (f1Last[kNeighbor0PP]*d0PP + f1Last[kNeighbor0PM]*d0PM + f1Last[kNeighbor0MP]*d0MP + f1Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f1NextInterp = (f1Next[kNeighbor0PP]*d0PP + f1Next[kNeighbor0PM]*d0PM + f1Next[kNeighbor0MP]*d0MP + f1Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f2LastInterp = (f2Last[kNeighbor0PP]*d0PP + f2Last[kNeighbor0PM]*d0PM + f2Last[kNeighbor0MP]*d0MP + f2Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f2NextInterp = (f2Next[kNeighbor0PP]*d0PP + f2Next[kNeighbor0PM]*d0PM + f2Next[kNeighbor0MP]*d0MP + f2Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f3LastInterp = (f3Last[kNeighbor0PP]*d0PP + f3Last[kNeighbor0PM]*d0PM + f3Last[kNeighbor0MP]*d0MP + f3Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f3NextInterp = (f3Next[kNeighbor0PP]*d0PP + f3Next[kNeighbor0PM]*d0PM + f3Next[kNeighbor0MP]*d0MP + f3Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f4LastInterp = (f4Last[kNeighbor0PP]*d0PP + f4Last[kNeighbor0PM]*d0PM + f4Last[kNeighbor0MP]*d0MP + f4Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f4NextInterp = (f4Next[kNeighbor0PP]*d0PP + f4Next[kNeighbor0PM]*d0PM + f4Next[kNeighbor0MP]*d0MP + f4Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f5LastInterp = (f5Last[kNeighbor0PP]*d0PP + f5Last[kNeighbor0PM]*d0PM + f5Last[kNeighbor0MP]*d0MP + f5Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f5NextInterp = (f5Next[kNeighbor0PP]*d0PP + f5Next[kNeighbor0PM]*d0PM + f5Next[kNeighbor0MP]*d0MP + f5Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f6LastInterp = (f6Last[kNeighbor0PP]*d0PP + f6Last[kNeighbor0PM]*d0PM + f6Last[kNeighbor0MP]*d0MP + f6Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f6NextInterp = (f6Next[kNeighbor0PP]*d0PP + f6Next[kNeighbor0PM]*d0PM + f6Next[kNeighbor0MP]*d0MP + f6Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f7LastInterp = (f7Last[kNeighbor0PP]*d0PP + f7Last[kNeighbor0PM]*d0PM + f7Last[kNeighbor0MP]*d0MP + f7Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f7NextInterp = (f7Next[kNeighbor0PP]*d0PP + f7Next[kNeighbor0PM]*d0PM + f7Next[kNeighbor0MP]*d0MP + f7Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f8LastInterp = (f8Last[kNeighbor0PP]*d0PP + f8Last[kNeighbor0PM]*d0PM + f8Last[kNeighbor0MP]*d0MP + f8Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f8NextInterp = (f8Next[kNeighbor0PP]*d0PP + f8Next[kNeighbor0PM]*d0PM + f8Next[kNeighbor0MP]*d0MP + f8Next[kNeighbor0MM]*d0MM)*invWeightSum;

    } else {
        f0LastInterp = f0Last[kNeighbor0PP];
        f1LastInterp = f1Last[kNeighbor0PP];
        f2LastInterp = f2Last[kNeighbor0PP];
        f3LastInterp = f3Last[kNeighbor0PP];
        f4LastInterp = f4Last[kNeighbor0PP];
        f5LastInterp = f5Last[kNeighbor0PP];
        f6LastInterp = f6Last[kNeighbor0PP];
        f7LastInterp = f7Last[kNeighbor0PP];
        f8LastInterp = f8Last[kNeighbor0PP];

        f0NextInterp = f0Next[kNeighbor0PP];
        f1NextInterp = f1Next[kNeighbor0PP];
        f2NextInterp = f2Next[kNeighbor0PP];
        f3NextInterp = f3Next[kNeighbor0PP];
        f4NextInterp = f4Next[kNeighbor0PP];
        f5NextInterp = f5Next[kNeighbor0PP];
        f6NextInterp = f6Next[kNeighbor0PP];
        f7NextInterp = f7Next[kNeighbor0PP];
        f8NextInterp = f8Next[kNeighbor0PP];
    }
    //////////////////////////////////////////////////////////////////////////
    //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep
    //! is based on the esoteric twist algorithm \ref <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier
    //! et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
    //!
    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

    unsigned int KQK  = subgridDistanceIndices[nodeIndex];
    // unsigned int k000= KQK;
    unsigned int kP00   = KQK;
    // unsigned int kM00   = neighborX[KQK];
    // unsigned int k0P0   = KQK;
    unsigned int k0M0   = neighborY[KQK];
    // unsigned int k00P   = KQK;
    unsigned int k00M   = neighborZ[KQK];
    // unsigned int kMM0  = neighborY[kM00];
    unsigned int kPP0  = KQK;
    unsigned int kPM0  = k0M0;
    // unsigned int kMP0  = kM00;
    // unsigned int kM0M  = neighborZ[kM00];
    unsigned int kP0P  = KQK;
    unsigned int kP0M  = k00M;
    // unsigned int kM0P  = kM00;
    unsigned int k0MM  = neighborZ[k0M0];
    // unsigned int k0PM  = k00M;
    // unsigned int k0MP  = k0M0;
    unsigned int kPMP = k0M0;
    // unsigned int kMPM = kM0M;
    // unsigned int kMPP = kM00;
    unsigned int kPMM = k0MM;
    // unsigned int kMMP = kMM0;
    unsigned int kPPM = k00M;
    unsigned int kPPP = KQK;
    // unsigned int kMMM = neighborZ[kMM0];

    dist.f[DIR_P00][kP00] = f0LastInterp*(1.f-timeRatio) + f0NextInterp*timeRatio;
    dist.f[DIR_PP0][kPP0] = f1LastInterp*(1.f-timeRatio) + f1NextInterp*timeRatio;
    dist.f[DIR_PM0][kPM0] = f2LastInterp*(1.f-timeRatio) + f2NextInterp*timeRatio;
    dist.f[DIR_P0P][kP0P] = f3LastInterp*(1.f-timeRatio) + f3NextInterp*timeRatio;
    dist.f[DIR_P0M][kP0M] = f4LastInterp*(1.f-timeRatio) + f4NextInterp*timeRatio;
    dist.f[DIR_PPP][kPPP] = f5LastInterp*(1.f-timeRatio) + f5NextInterp*timeRatio;
    dist.f[DIR_PMP][kPMP] = f6LastInterp*(1.f-timeRatio) + f6NextInterp*timeRatio;
    dist.f[DIR_PPM][kPPM] = f7LastInterp*(1.f-timeRatio) + f7NextInterp*timeRatio;
    dist.f[DIR_PMM][kPMM] = f8LastInterp*(1.f-timeRatio) + f8NextInterp*timeRatio;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////












































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// NOTE: Has not been tested after bug fix!
__global__ void QPrecursorDeviceDistributions(
    int* subgridDistanceIndices,
    real* subgridDistances,
    int sizeQ,
    int numberOfBCnodes,
    int numberOfPrecursorNodes,
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
    real* fsLast,
    real* fsNext,
    real timeRatio,
    unsigned long long numberOfLBnodes,
    bool isEvenTimestep)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = getNodeIndex();

    if(nodeIndex>=numberOfBCnodes) return;

    uint kNeighbor0PP = neighbors0PP[nodeIndex];
    real d0PP = weights0PP[nodeIndex];

    real f0LastInterp, f1LastInterp, f2LastInterp, f3LastInterp, f4LastInterp, f5LastInterp, f6LastInterp, f7LastInterp, f8LastInterp;
    real f0NextInterp, f1NextInterp, f2NextInterp, f3NextInterp, f4NextInterp, f5NextInterp, f6NextInterp, f7NextInterp, f8NextInterp;

    real* f0Last = fsLast;
    real* f1Last = &fsLast[  numberOfPrecursorNodes];
    real* f2Last = &fsLast[2*numberOfPrecursorNodes];
    real* f3Last = &fsLast[3*numberOfPrecursorNodes];
    real* f4Last = &fsLast[4*numberOfPrecursorNodes];
    real* f5Last = &fsLast[5*numberOfPrecursorNodes];
    real* f6Last = &fsLast[6*numberOfPrecursorNodes];
    real* f7Last = &fsLast[7*numberOfPrecursorNodes];
    real* f8Last = &fsLast[8*numberOfPrecursorNodes];

    real* f0Next = fsNext;
    real* f1Next = &fsNext[  numberOfPrecursorNodes];
    real* f2Next = &fsNext[2*numberOfPrecursorNodes];
    real* f3Next = &fsNext[3*numberOfPrecursorNodes];
    real* f4Next = &fsNext[4*numberOfPrecursorNodes];
    real* f5Next = &fsNext[5*numberOfPrecursorNodes];
    real* f6Next = &fsNext[6*numberOfPrecursorNodes];
    real* f7Next = &fsNext[7*numberOfPrecursorNodes];
    real* f8Next = &fsNext[8*numberOfPrecursorNodes];


    if(d0PP<1e6)
    {
        uint kNeighbor0PM = neighbors0PM[nodeIndex];
        uint kNeighbor0MP = neighbors0MP[nodeIndex];
        uint kNeighbor0MM = neighbors0MM[nodeIndex];

        real d0PM = weights0PM[nodeIndex];
        real d0MP = weights0MP[nodeIndex];
        real d0MM = weights0MM[nodeIndex];

        real invWeightSum = 1.f/(d0PP+d0PM+d0MP+d0MM);

        f0LastInterp = (f0Last[kNeighbor0PP]*d0PP + f0Last[kNeighbor0PM]*d0PM + f0Last[kNeighbor0MP]*d0MP + f0Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f0NextInterp = (f0Next[kNeighbor0PP]*d0PP + f0Next[kNeighbor0PM]*d0PM + f0Next[kNeighbor0MP]*d0MP + f0Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f1LastInterp = (f1Last[kNeighbor0PP]*d0PP + f1Last[kNeighbor0PM]*d0PM + f1Last[kNeighbor0MP]*d0MP + f1Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f1NextInterp = (f1Next[kNeighbor0PP]*d0PP + f1Next[kNeighbor0PM]*d0PM + f1Next[kNeighbor0MP]*d0MP + f1Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f2LastInterp = (f2Last[kNeighbor0PP]*d0PP + f2Last[kNeighbor0PM]*d0PM + f2Last[kNeighbor0MP]*d0MP + f2Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f2NextInterp = (f2Next[kNeighbor0PP]*d0PP + f2Next[kNeighbor0PM]*d0PM + f2Next[kNeighbor0MP]*d0MP + f2Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f3LastInterp = (f3Last[kNeighbor0PP]*d0PP + f3Last[kNeighbor0PM]*d0PM + f3Last[kNeighbor0MP]*d0MP + f3Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f3NextInterp = (f3Next[kNeighbor0PP]*d0PP + f3Next[kNeighbor0PM]*d0PM + f3Next[kNeighbor0MP]*d0MP + f3Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f4LastInterp = (f4Last[kNeighbor0PP]*d0PP + f4Last[kNeighbor0PM]*d0PM + f4Last[kNeighbor0MP]*d0MP + f4Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f4NextInterp = (f4Next[kNeighbor0PP]*d0PP + f4Next[kNeighbor0PM]*d0PM + f4Next[kNeighbor0MP]*d0MP + f4Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f5LastInterp = (f5Last[kNeighbor0PP]*d0PP + f5Last[kNeighbor0PM]*d0PM + f5Last[kNeighbor0MP]*d0MP + f5Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f5NextInterp = (f5Next[kNeighbor0PP]*d0PP + f5Next[kNeighbor0PM]*d0PM + f5Next[kNeighbor0MP]*d0MP + f5Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f6LastInterp = (f6Last[kNeighbor0PP]*d0PP + f6Last[kNeighbor0PM]*d0PM + f6Last[kNeighbor0MP]*d0MP + f6Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f6NextInterp = (f6Next[kNeighbor0PP]*d0PP + f6Next[kNeighbor0PM]*d0PM + f6Next[kNeighbor0MP]*d0MP + f6Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f7LastInterp = (f7Last[kNeighbor0PP]*d0PP + f7Last[kNeighbor0PM]*d0PM + f7Last[kNeighbor0MP]*d0MP + f7Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f7NextInterp = (f7Next[kNeighbor0PP]*d0PP + f7Next[kNeighbor0PM]*d0PM + f7Next[kNeighbor0MP]*d0MP + f7Next[kNeighbor0MM]*d0MM)*invWeightSum;

        f8LastInterp = (f8Last[kNeighbor0PP]*d0PP + f8Last[kNeighbor0PM]*d0PM + f8Last[kNeighbor0MP]*d0MP + f8Last[kNeighbor0MM]*d0MM)*invWeightSum;
        f8NextInterp = (f8Next[kNeighbor0PP]*d0PP + f8Next[kNeighbor0PM]*d0PM + f8Next[kNeighbor0MP]*d0MP + f8Next[kNeighbor0MM]*d0MM)*invWeightSum;

    } else {
        f0LastInterp = f0Last[kNeighbor0PP];
        f1LastInterp = f1Last[kNeighbor0PP];
        f2LastInterp = f2Last[kNeighbor0PP];
        f3LastInterp = f3Last[kNeighbor0PP];
        f4LastInterp = f4Last[kNeighbor0PP];
        f5LastInterp = f5Last[kNeighbor0PP];
        f6LastInterp = f6Last[kNeighbor0PP];
        f7LastInterp = f7Last[kNeighbor0PP];
        f8LastInterp = f8Last[kNeighbor0PP];

        f0NextInterp = f0Next[kNeighbor0PP];
        f1NextInterp = f1Next[kNeighbor0PP];
        f2NextInterp = f2Next[kNeighbor0PP];
        f3NextInterp = f3Next[kNeighbor0PP];
        f4NextInterp = f4Next[kNeighbor0PP];
        f5NextInterp = f5Next[kNeighbor0PP];
        f6NextInterp = f6Next[kNeighbor0PP];
        f7NextInterp = f7Next[kNeighbor0PP];
        f8NextInterp = f8Next[kNeighbor0PP];
    }
    //////////////////////////////////////////////////////////////////////////
    //! - Read distributions: style of reading and writing the distributions from/to stored arrays dependent on timestep
    //! is based on the esoteric twist algorithm \ref <a href="https://doi.org/10.3390/computation5020019"><b>[ M. Geier
    //! et al. (2017), DOI:10.3390/computation5020019 ]</b></a>
    //!
    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, !isEvenTimestep);

    unsigned int KQK  = subgridDistanceIndices[nodeIndex];
    // unsigned int k000= KQK;
    unsigned int kP00   = KQK;
    // unsigned int kM00   = neighborX[KQK];
    // unsigned int k0P0   = KQK;
    unsigned int k0M0   = neighborY[KQK];
    // unsigned int k00P   = KQK;
    unsigned int k00M   = neighborZ[KQK];
    // unsigned int kMM0  = neighborY[kM00];
    unsigned int kPP0  = KQK;
    unsigned int kPM0  = k0M0;
    // unsigned int kMP0  = kM00;
    // unsigned int kM0M  = neighborZ[kM00];
    unsigned int kP0P  = KQK;
    unsigned int kP0M  = k00M;
    // unsigned int kM0P  = kM00;
    unsigned int k0MM  = neighborZ[k0M0];
    // unsigned int k0PM  = k00M;
    // unsigned int k0MP  = k0M0;
    unsigned int kPMP = k0M0;
    // unsigned int kMPM = kM0M;
    // unsigned int kMPP = kM00;
    unsigned int kPMM = k0MM;
    // unsigned int kMMP = kMM0;
    unsigned int kPPM = k00M;
    unsigned int kPPP = KQK;
    // unsigned int kMMM = neighborZ[kMM0];
    SubgridDistances27 qs;
    getPointersToSubgridDistances(qs, subgridDistances, sizeQ);

    real q;
    q = qs.q[DIR_P00][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_P00][kP00] = f0LastInterp*(1.f-timeRatio) + f0NextInterp*timeRatio;
    q = qs.q[DIR_PP0][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_PP0][kPP0] = f1LastInterp*(1.f-timeRatio) + f1NextInterp*timeRatio;
    q = qs.q[DIR_PM0][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_PM0][kPM0] = f2LastInterp*(1.f-timeRatio) + f2NextInterp*timeRatio;
    q = qs.q[DIR_P0P][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_P0P][kP0P] = f3LastInterp*(1.f-timeRatio) + f3NextInterp*timeRatio;
    q = qs.q[DIR_P0M][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_P0M][kP0M] = f4LastInterp*(1.f-timeRatio) + f4NextInterp*timeRatio;
    q = qs.q[DIR_PPP][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_PPP][kPPP] = f5LastInterp*(1.f-timeRatio) + f5NextInterp*timeRatio;
    q = qs.q[DIR_PMP][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_PMP][kPMP] = f6LastInterp*(1.f-timeRatio) + f6NextInterp*timeRatio;
    q = qs.q[DIR_PPM][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_PPM][kPPM] = f7LastInterp*(1.f-timeRatio) + f7NextInterp*timeRatio;
    q = qs.q[DIR_PMM][nodeIndex]; if(q>= c0o1 && q <= c1o1) dist.f[DIR_PMM][kPMM] = f8LastInterp*(1.f-timeRatio) + f8NextInterp*timeRatio;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
