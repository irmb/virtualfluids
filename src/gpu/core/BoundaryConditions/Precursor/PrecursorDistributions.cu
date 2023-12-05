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
#include "Calculation/Calculation.h"
#include <basics/constants/NumericConstants.h>
#include <lbm/constants/D3Q27.h>
#include <lbm/MacroscopicQuantities.h>

#include "Utilities/KernelUtilities.h"

using namespace vf::basics::constant;
using namespace vf::lbm::dir;
using namespace vf::gpu;

__global__ void PrecursorDistributions_Device(
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

    dist.f[dP00][kP00] = f0LastInterp*(1.f-timeRatio) + f0NextInterp*timeRatio;
    dist.f[dPP0][kPP0] = f1LastInterp*(1.f-timeRatio) + f1NextInterp*timeRatio;
    dist.f[dPM0][kPM0] = f2LastInterp*(1.f-timeRatio) + f2NextInterp*timeRatio;
    dist.f[dP0P][kP0P] = f3LastInterp*(1.f-timeRatio) + f3NextInterp*timeRatio;
    dist.f[dP0M][kP0M] = f4LastInterp*(1.f-timeRatio) + f4NextInterp*timeRatio;
    dist.f[dPPP][kPPP] = f5LastInterp*(1.f-timeRatio) + f5NextInterp*timeRatio;
    dist.f[dPMP][kPMP] = f6LastInterp*(1.f-timeRatio) + f6NextInterp*timeRatio;
    dist.f[dPPM][kPPM] = f7LastInterp*(1.f-timeRatio) + f7NextInterp*timeRatio;
    dist.f[dPMM][kPMM] = f8LastInterp*(1.f-timeRatio) + f8NextInterp*timeRatio;
}
