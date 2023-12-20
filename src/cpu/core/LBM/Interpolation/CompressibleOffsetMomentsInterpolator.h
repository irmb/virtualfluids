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
//! \addtogroup cpu_LBM LBM
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef CompressibleOffsetMomentsInterpolationProcessor_H_
#define CompressibleOffsetMomentsInterpolationProcessor_H_

#include "Interpolator.h"
#include "D3Q27System.h"

//////////////////////////////////////////////////////////////////////////
//it works only for cascaded LBM
//super compact interpolation method by Martin Geier
//////////////////////////////////////////////////////////////////////////

class CompressibleOffsetMomentsInterpolator : public Interpolator
{
public:
    CompressibleOffsetMomentsInterpolator() = default;
    CompressibleOffsetMomentsInterpolator(real omegaC, real omegaF);
    InterpolationProcessorPtr clone() override;

    void setOmegas(real omegaC, real omegaF) override;

    void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF) override;
    void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF, real xoff, real yoff, real zoff) override;
    void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC) override;
    void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC, real xoff, real yoff, real zoff) override;

private:
    real omegaC{ 0.0 }, omegaF{ 0.0 };
};

//////////////////////////////////////////////////////////////////////////
inline void CompressibleOffsetMomentsInterpolator::interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF)
{
    this->interpolateCoarseToFine(icellC, icellF, 0.0, 0.0, 0.0);
}
//////////////////////////////////////////////////////////////////////////
inline void CompressibleOffsetMomentsInterpolator::interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC)
{
    this->interpolateFineToCoarse(icellF, icellC, 0.0, 0.0, 0.0);
}

#endif

//! \}
