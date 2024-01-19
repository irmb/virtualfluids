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
//! \addtogroup cpu_BoundaryConditions BoundaryConditions
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef BCStrategy_H
#define BCStrategy_H

#include <PointerDefinitions.h>

#include "D3Q27System.h"

class DistributionArray3D;
class BCArray3D;
class BoundaryConditions;
class Block3D;

//! \brief Abstract class of baundary conditions strategy
//! \details  BCStrategy provides interface for implementation of diferent boundary conditions
class BCStrategy
{
public:
    BCStrategy() = default;
    virtual ~BCStrategy() = default;

    virtual void addDistributions(SPtr<DistributionArray3D> distributions)   = 0;
    void setBlock(SPtr<Block3D> block);
    void setNodeIndex(int x1, int x2, int x3);
    void setBcPointer(SPtr<BoundaryConditions> bcPtr);
    void setCompressible(bool c);
    void setCollFactor(real cf);

    bool isPreCollision();
    virtual SPtr<BCStrategy> clone() = 0;
    SPtr<BCArray3D> getBcArray();
    void setBcArray(SPtr<BCArray3D> bcarray);
    virtual void applyBC() = 0;

protected:
    bool compressible { false };
    bool preCollision;

    SPtr<BoundaryConditions> bcPtr;
    SPtr<DistributionArray3D> distributions;
    SPtr<BCArray3D> bcArray;
    SPtr<Block3D> block;

    real collFactor;
    int x1, x2, x3;

    real compressibleFactor;

    using CalcMacrosFct    = void (*)(const real *const &, real &, real &, real &, real &);
    using CalcFeqForDirFct = real (*)(const int &, const real &, const real &, const real &,
                                         const real &);
    using CalcFeqFct = void (*)(real *const &, const real &, const real &, const real &, const real &);

    CalcFeqForDirFct calcFeqsForDirFct;
    CalcMacrosFct calcMacrosFct;
    CalcFeqFct calcFeqFct;
};

#endif

//! \}
