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

#ifndef CPU_INTERPOLATER_H
#define CPU_INTERPOLATER_H

#include <memory>

#include "BCArray3D.h"
#include "BoundaryConditions.h"
#include "DistributionArray3D.h"
#include "LBMSystem.h"

struct D3Q27ICell {
    real TSW[27];
    real TNW[27];
    real TNE[27];
    real TSE[27];
    real BSW[27];
    real BNW[27];
    real BNE[27];
    real BSE[27];
};

class Interpolator;
using InterpolationProcessorPtr = std::shared_ptr<Interpolator>;


class Interpolator
{
public:
    virtual ~Interpolator() = default;
    virtual InterpolationProcessorPtr clone() = 0;

    virtual void setOmegas(real omegaC, real omegaF)                       = 0;
    virtual void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF) = 0;
    virtual void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF, real xoff, real yoff,
                                         real zoff)                           = 0;
    virtual void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC)    = 0;
    virtual void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC, real xoff, real yoff,
                                         real zoff)                           = 0;

    static void readICell(SPtr<DistributionArray3D> f, D3Q27ICell &icell, int x1, int x2, int x3);
    static void writeICell(SPtr<DistributionArray3D> f, const D3Q27ICell &icell, int x1, int x2, int x3);
    static void writeICellInv(SPtr<DistributionArray3D> f, const D3Q27ICell &icell, int x1, int x2, int x3);
    static void writeINode(SPtr<DistributionArray3D> f, const real *const inode, int x1, int x2, int x3);
    static void writeINodeInv(SPtr<DistributionArray3D> f, const real *const inode, int x1, int x2, int x3);
    static bool iCellHasSolid(const SPtr<BCArray3D> bcArray, int x1, int x2, int x3);
    static int iCellHowManySolids(const SPtr<BCArray3D> bcArray, int x1, int x2, int x3);

    bool findNeighborICell(const SPtr<BCArray3D> bcArray, SPtr<DistributionArray3D> f, D3Q27ICell &icell, int maxX1,
                           int maxX2, int maxX3, int x1, int x2, int x3, real &xoff, real &yoff, real &zoff);

protected:
    virtual void calcInterpolatedCoefficiets(const D3Q27ICell &icell, real omega, real eps_new) {}
    virtual void calcInterpolatedNodeFC(real *f, real omega) {}
    virtual void calcInterpolatedVelocity(real x, real y, real z, real &vx1, real &vx2, real &vx3) {}
    virtual void calcInterpolatedShearStress(real x, real y, real z, real &tauxx, real &tauyy,
                                             real &tauzz, real &tauxy, real &tauxz, real &tauyz) {}
    virtual void setOffsets(real xoff, real yoff, real zoff) {}

};

#endif

//! \}
