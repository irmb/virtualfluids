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
//! \addtogroup cpu_Visitors Visitors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef InitDistributionsWithCoarseGridBlockVisitor_h__
#define InitDistributionsWithCoarseGridBlockVisitor_h__

#include <PointerDefinitions.h>

#include "Grid3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;
class Block3D;
class Interpolator;

class InitDistributionsWithInterpolationGridVisitor : public Grid3DVisitor
{
public:
    InitDistributionsWithInterpolationGridVisitor(SPtr<Grid3D> oldGrid, SPtr<Interpolator> iProcessor,
                                                  real nu);
    ~InitDistributionsWithInterpolationGridVisitor() override;
    void visit(SPtr<Grid3D> grid) override;

private:
    void copyLocalBlock(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
    void copyRemoteBlock(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
    void interpolateLocalBlockCoarseToFine(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
    void interpolateRemoteBlockCoarseToFine(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
    void interpolateLocalBlockFineToCoarse(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);
    void interpolateRemoteBlockFineToCoarse(SPtr<Block3D> oldBlock, SPtr<Block3D> newBlock);

    SPtr<Grid3D> newGrid;
    SPtr<Grid3D> oldGrid;
    real nu;

    SPtr<Interpolator> iProcessor;
};

#endif // InitDistributionsWithVelocityProfileBlockVisitor_h__

//! \}
