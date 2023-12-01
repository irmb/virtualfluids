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
//! \file CompressibleOffsetMomentsInterpolator.cpp
//! \ingroup Interpolation
//! \author Konstantin Kutscher
//=======================================================================================

#include "CompressibleOffsetMomentsInterpolator.h"

#include <algorithm>

#include <basics/constants/NumericConstants.h>

#include <lbm/refinement/InterpolationCF.h>
#include <lbm/refinement/InterpolationFC.h>
#include <lbm/interpolation/InterpolationCoefficients.h>


void calculateCoefficients(vf::lbm::InterpolationCoefficients& coefficients, const D3Q27ICell& icell, real omega, real xoff, real yoff, real zoff)
{
    vf::lbm::MomentsOnSourceNodeSet momentsSet;

    momentsSet.calculateMMM(icell.BSW, omega);
    momentsSet.calculateMMP(icell.TSW, omega);
    momentsSet.calculateMPP(icell.TNW, omega);
    momentsSet.calculateMPM(icell.BNW, omega);
    momentsSet.calculatePMM(icell.BSE, omega);
    momentsSet.calculatePPP(icell.TNE, omega);
    momentsSet.calculatePMP(icell.TSE, omega);
    momentsSet.calculatePPM(icell.BNE, omega);

    momentsSet.calculateCoefficients(coefficients, xoff, yoff, zoff);
}

CompressibleOffsetMomentsInterpolator::CompressibleOffsetMomentsInterpolator(real omegaC, real omegaF)
   : omegaC(omegaC), omegaF(omegaF)
{
}

InterpolationProcessorPtr CompressibleOffsetMomentsInterpolator::clone()
{
   return InterpolationProcessorPtr (new CompressibleOffsetMomentsInterpolator(this->omegaC, this->omegaF));
}

void CompressibleOffsetMomentsInterpolator::setOmegas(real omegaC, real omegaF)
{
   this->omegaC = omegaC;
   this->omegaF = omegaF;
}

void CompressibleOffsetMomentsInterpolator::interpolateCoarseToFine(D3Q27ICell& icellC, D3Q27ICell& icellF, real xoff, real yoff, real zoff)
{
    vf::lbm::InterpolationCoefficients coefficients;
    calculateCoefficients(coefficients, icellC, omegaC, xoff, yoff, zoff);

     vf::lbm::interpolateCF(icellF.BSW, omegaF, vf::basics::constant::c1o2, coefficients, -c1o4, -c1o4, -c1o4);
     vf::lbm::interpolateCF(icellF.BNE, omegaF, vf::basics::constant::c1o2, coefficients,  c1o4,  c1o4, -c1o4);
     vf::lbm::interpolateCF(icellF.TNW, omegaF, vf::basics::constant::c1o2, coefficients, -c1o4,  c1o4,  c1o4);
     vf::lbm::interpolateCF(icellF.TSE, omegaF, vf::basics::constant::c1o2, coefficients,  c1o4, -c1o4,  c1o4);
     vf::lbm::interpolateCF(icellF.BNW, omegaF, vf::basics::constant::c1o2, coefficients, -c1o4,  c1o4, -c1o4);
     vf::lbm::interpolateCF(icellF.BSE, omegaF, vf::basics::constant::c1o2, coefficients,  c1o4, -c1o4, -c1o4);
     vf::lbm::interpolateCF(icellF.TSW, omegaF, vf::basics::constant::c1o2, coefficients, -c1o4, -c1o4,  c1o4);
     vf::lbm::interpolateCF(icellF.TNE, omegaF, vf::basics::constant::c1o2, coefficients,  c1o4,  c1o4,  c1o4);
}

void CompressibleOffsetMomentsInterpolator::interpolateFineToCoarse(D3Q27ICell& icellF, real* icellC, real xoff, real yoff, real zoff)
{
    vf::lbm::InterpolationCoefficients coefficients;
    calculateCoefficients(coefficients, icellF, omegaF, xoff, yoff, zoff);

    vf::lbm::interpolateFC(icellC, vf::basics::constant::c2o1, omegaC, coefficients);
}
