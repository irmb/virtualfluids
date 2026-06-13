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
//! \addtogroup interpolation
//! \ingroup lbm
//! \{
//=======================================================================================
#ifndef ADVECTION_DIFFUSION_SCALING_HELPER_FUNCTIONS_H
#define ADVECTION_DIFFUSION_SCALING_HELPER_FUNCTIONS_H

#include <basics/constants/NumericConstants.h>

#include "lbm/constants/D3Q27.h"

#include "lbm/MacroscopicQuantities.h"

//using namespace vf::basics::constant;

namespace vf::lbm::ad
{

// The Coefficients struct needs be created like this:
// MomentsOnSourceNodeSet momentsSet;
// momentsSet.calculatePPP(f, omega);
// ... and so on
// InterpolationCoefficients coeffs;
// momentsSet.calculateCoefficients(coeffs);

// Coefficients of the interpolation polynomial
struct InterpolationCoefficients
{
    real d000, d100, d010, d001, d200, d020, d002, d110, d101, d011, d111;
};

// Private struct - is only used within the MomentsOnSourceNodeSet
struct MomentsOnSourceNode
{
    real concentration;
    real dxConcentration;
    real dyConcentration;
    real dzConcentration;

    constexpr void calculate(const real* const populations, const real* const populationsAD, const real omega)
    {
        using namespace vf::basics::constant;
        real drho = c0o1;
        real velocity100 = c0o1;
        real velocity010 = c0o1;
        real velocity001 = c0o1;
        real oneOverRho = c0o1;
        vf::lbm::getCompressibleMacroscopicValues(populations, drho, oneOverRho, velocity100, velocity010, velocity001);
        this->concentration = vf::lbm::getDensity(populationsAD);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Calculate first order moments for interpolation
        //!
        const real momentAdvectionDiffusion100 = vf::lbm::getIncompressibleVelocityX1(populationsAD);
        const real momentAdvectionDiffusion010 = vf::lbm::getIncompressibleVelocityX2(populationsAD);
        const real momentAdvectionDiffusion001 = vf::lbm::getIncompressibleVelocityX3(populationsAD);

        ////////////////////////////////////////////////////////////////////////////////////
        //! - Calculate gradients for interpolation
        //!
        this->dxConcentration = (this->concentration * velocity100 - momentAdvectionDiffusion100) * (c3o1 * omega);
        this->dyConcentration = (this->concentration * velocity010 - momentAdvectionDiffusion010) * (c3o1 * omega);
        this->dzConcentration = (this->concentration * velocity001 - momentAdvectionDiffusion001) * (c3o1 * omega);
    }

};

class MomentsOnSourceNodeSet
{
private:
    vf::lbm::ad::MomentsOnSourceNode momentsPPP;
    vf::lbm::ad::MomentsOnSourceNode momentsMPP;
    vf::lbm::ad::MomentsOnSourceNode momentsPMP;
    vf::lbm::ad::MomentsOnSourceNode momentsMMP;
    vf::lbm::ad::MomentsOnSourceNode momentsPPM;
    vf::lbm::ad::MomentsOnSourceNode momentsMPM;
    vf::lbm::ad::MomentsOnSourceNode momentsPMM;
    vf::lbm::ad::MomentsOnSourceNode momentsMMM;

public:
    constexpr void calculatePPP(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsPPP.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculateMPP(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsMPP.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculatePMP(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsPMP.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculateMMP(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsMMP.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculatePPM(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsPPM.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculateMPM(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsMPM.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculatePMM(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsPMM.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculateMMM(const real* const populations, const real* const populationsAD, const real omegaDiffusivity)
    {
        momentsMMM.calculate(populations, populationsAD, omegaDiffusivity);
    }

    constexpr void calculateCoefficients(InterpolationCoefficients &coefficients, real xoff, real yoff, real zoff) const
    {
        using namespace vf::basics::constant;
        real &d000 = coefficients.d000;
        real &d100 = coefficients.d100, &d010 = coefficients.d010, &d001 = coefficients.d001;
        real &d200 = coefficients.d200, &d020 = coefficients.d020, &d002 = coefficients.d002;
        real &d110 = coefficients.d110, &d101 = coefficients.d101, &d011 = coefficients.d011;
        real &d111 = coefficients.d111;

        const real xoffsq = xoff * xoff;
        const real yoffsq = yoff * yoff;
        const real zoffsq = zoff * zoff;

        const real dxxConcentration = (((momentsPPM.dxConcentration + momentsPMP.dxConcentration) + (momentsPMM.dxConcentration + momentsPPP.dxConcentration)) -
                                       ((momentsMPM.dxConcentration + momentsMMP.dxConcentration) + (momentsMMM.dxConcentration + momentsMPP.dxConcentration))) * c1o4;
        const real dyyConcentration = (((momentsPPM.dyConcentration + momentsMPP.dyConcentration) + (momentsMPM.dyConcentration + momentsPPP.dyConcentration)) -
                                       ((momentsPMM.dyConcentration + momentsMMP.dyConcentration) + (momentsMMM.dyConcentration + momentsPMP.dyConcentration))) * c1o4;
        const real dzzConcentration = (((momentsPMP.dzConcentration + momentsMPP.dzConcentration) + (momentsMMP.dzConcentration + momentsPPP.dzConcentration)) -
                                       ((momentsPMM.dzConcentration + momentsMPM.dzConcentration) + (momentsMMM.dzConcentration + momentsPPM.dzConcentration))) * c1o4;
            


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //!- Calculate coefficients for the polynomial interpolation of the pressure
        //! 

        d100 = c1o4 * (((momentsPPM.concentration - momentsMMP.concentration) + (momentsPMP.concentration - momentsMPM.concentration)) + 
                       ((momentsPPP.concentration - momentsMMM.concentration) + (momentsPMM.concentration - momentsMPP.concentration)));
        d010 = c1o4 * (((momentsPPM.concentration - momentsMMP.concentration) + (momentsMPP.concentration - momentsPMM.concentration)) + 
                       ((momentsPPP.concentration - momentsMMM.concentration) + (momentsMPM.concentration - momentsPMP.concentration)));
        d001 = c1o4 * (((momentsPMP.concentration - momentsMPM.concentration) + (momentsMPP.concentration - momentsPMM.concentration)) + 
                       ((momentsPPP.concentration - momentsMMM.concentration) + (momentsMMP.concentration - momentsPPM.concentration)));
        d200 = c1o2 * dxxConcentration;
        d020 = c1o2 * dyyConcentration;
        d002 = c1o2 * dzzConcentration;
        d110 = c1o2 * (((momentsPPP.concentration - momentsPMP.concentration) + (momentsPPM.concentration - momentsPMM.concentration)) + 
                       ((momentsMMP.concentration - momentsMPP.concentration) + (momentsMMM.concentration - momentsMPM.concentration)));
        d101 = c1o2 * (((momentsPPP.concentration - momentsPPM.concentration) + (momentsPMP.concentration - momentsPMM.concentration)) + 
                       ((momentsMPM.concentration - momentsMPP.concentration) + (momentsMMM.concentration - momentsMMP.concentration)));
        d011 = c1o2 * (((momentsPPP.concentration - momentsPPM.concentration) + (momentsMPP.concentration - momentsMPM.concentration)) + 
                       ((momentsPMM.concentration - momentsPMP.concentration) + (momentsMMM.concentration - momentsMMP.concentration)));

        d111 = ((momentsPPP.concentration - momentsMMM.concentration) + (momentsPMM.concentration - momentsMPP.concentration)) +
               ((momentsMPM.concentration - momentsPMP.concentration) + (momentsMMP.concentration - momentsPPM.concentration));
            
        d000 = c1o8 * (( - dxxConcentration - dyyConcentration - dzzConcentration) + 
                (((momentsPPP.concentration + momentsMMM.concentration) + (momentsPMM.concentration + momentsMPP.concentration)) +
                ((momentsMPM.concentration + momentsPMP.concentration) + (momentsMMP.concentration + momentsPPM.concentration))));

        //////////////////////////////////////////////////////////////////////////
        //! - Extrapolation for refinement in to the wall (polynomial coefficients)
        //!
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Coarse to Fine
        // X------X
        // |      | x---x
        // |   ---+-+-> |    ----> off-vector
        // |      | x---x
        // X------X
        ///////////////////////////////////////////
        // Fine to Coarse
        // x------x
        // |      |
        // |   ---+--->X
        // |      |  |
        // x------x  |
        //          offset-vector
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        d000 = d000 + 
               xoff * d100 + yoff * d010 + zoff * d001 + 
               xoffsq * d200 + yoffsq * d020 + zoffsq * d002 +
               xoff * yoff * d110 + xoff * zoff * d101 + yoff * zoff * d011;
        d100 = d100 + c2o1 * xoff * d200 + yoff * d110 + zoff * d101;
        d010 = d010 + c2o1 * yoff * d020 + xoff * d110 + zoff * d011;
        d001 = d001 + c2o1 * zoff * d002 + xoff * d101 + yoff * d011;

    }

};

} // namespace vf::lbm::ad

#endif

//! \}
