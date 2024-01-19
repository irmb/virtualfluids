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
#include "SetUndefinedNodesBlockVisitor.h"
#include "BCArray3D.h"
#include "BCSet.h"
#include "Block3D.h"
#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "Grid3D.h"
#include "D3Q27System.h"
#include "ILBMKernel.h"

SetUndefinedNodesBlockVisitor::SetUndefinedNodesBlockVisitor(bool twoTypeOfConnectorsCheck)
    : Block3DVisitor(0, D3Q27System::MAXLEVEL), twoTypeOfConnectorsCheck(twoTypeOfConnectorsCheck)
{
}
//////////////////////////////////////////////////////////////////////////
void SetUndefinedNodesBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
    using namespace vf::lbm::dir;

    if (!block->hasInterpolationFlag())
        return;

    SPtr<ILBMKernel> kernel = block->getKernel();

    if (!kernel && (block->getRank() != grid->getRank()))
        return;

    // width of ghost layer
    // int gl = kernel->getGhostLayerWidth();
    int gl = 0;

    SPtr<BCArray3D> bcMatrix = kernel->getBCSet()->getBCArray();

    int minX1 = gl;
    int minX2 = gl;
    int minX3 = gl;

    int maxX1 = static_cast<int>(bcMatrix->getNX1()) - 1 - gl;
    int maxX2 = static_cast<int>(bcMatrix->getNX2()) - 1 - gl;
    int maxX3 = static_cast<int>(bcMatrix->getNX3()) - 1 - gl;

    // int offset = 2;
    int offset = 3;

    if (block->hasInterpolationFlag(dP00)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dM00)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d0P0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d0M0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d00P)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d00M)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dPP0)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dMM0)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dPM0)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dMP0)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dP0P)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dM0M)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dP0M)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dM0P)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d0PP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d0MM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d0PM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(d0MP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dPPP)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dMPP)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dPMP)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dMMP)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = maxX3;
        int endix3   = maxX3;
        if (block->hasInterpolationFlagCF())
            startix3 = startix3 - offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dPPM)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dMPM)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = maxX2;
        int endix2   = maxX2;
        if (block->hasInterpolationFlagCF())
            startix2 = startix2 - offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dPMM)) {
        int startix1 = maxX1;
        int endix1   = maxX1;
        if (block->hasInterpolationFlagCF())
            startix1 = startix1 - offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlag(dMMM)) {
        int startix1 = minX1;
        int endix1   = minX1;
        if (block->hasInterpolationFlagCF())
            endix1 = endix1 + offset;
        int startix2 = minX2;
        int endix2   = minX2;
        if (block->hasInterpolationFlagCF())
            endix2 = endix2 + offset;
        int startix3 = minX3;
        int endix3   = minX3;
        if (block->hasInterpolationFlagCF())
            endix3 = endix3 + offset;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }

    //////////////////////////////////////////////////////////////////////////
    int offset2 = 1;
    int ll      = 0;

    minX1 = ll;
    minX2 = ll;
    minX3 = ll;

    maxX1 = static_cast<int>(bcMatrix->getNX1()) - 1 - ll;
    maxX2 = static_cast<int>(bcMatrix->getNX2()) - 1 - ll;
    maxX3 = static_cast<int>(bcMatrix->getNX3()) - 1 - ll;

    if (block->hasInterpolationFlagFC(dP00)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dM00)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d0P0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d0M0)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d00P)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d00M)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dPP0)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dMM0)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dPM0)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dMP0)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dP0P)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dM0M)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dP0M)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dM0P)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d0PP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d0MM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d0PM)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(d0MP)) {
        int startix1 = minX1;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dPPP)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dMPP)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dPMP)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dMMP)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = maxX3 - offset2;
        int endix3   = maxX3;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dPPM)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dMPM)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = maxX2 - offset2;
        int endix2   = maxX2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dPMM)) {
        int startix1 = maxX1 - offset2;
        int endix1   = maxX1;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }
    if (block->hasInterpolationFlagFC(dMMM)) {
        int startix1 = minX1;
        int endix1   = minX1 + offset2;
        int startix2 = minX2;
        int endix2   = minX2 + offset2;
        int startix3 = minX3;
        int endix3   = minX3 + offset2;
        this->setNodesUndefined(startix1, endix1, startix2, endix2, startix3, endix3, bcMatrix);
    }

    // invert scaleCF blocks
    if (block->hasInterpolationFlagCF()) {
        if (block->hasInterpolationFlagFC() && twoTypeOfConnectorsCheck) {
            for (int i = (int)dP00; i <= (int)dMMM; i++) {
                UBLOG(logINFO, "FC in dir=" << i << " " << block->hasInterpolationFlagFC(i));
            }
            for (int i = (int)dP00; i <= (int)dMMM; i++) {
                UBLOG(logINFO, "CF in dir=" << i << " " << block->hasInterpolationFlagCF(i));
            }
            throw UbException(UB_EXARGS, "block " + block->toString() + " has CF and FC");
        }

        minX1 = gl;
        minX2 = gl;
        minX3 = gl;

        maxX1 = static_cast<int>(bcMatrix->getNX1()) - 1 - gl;
        maxX2 = static_cast<int>(bcMatrix->getNX2()) - 1 - gl;
        maxX3 = static_cast<int>(bcMatrix->getNX3()) - 1 - gl;

        for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
            for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
                for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                    if (bcMatrix->isUndefined(ix1, ix2, ix3))
                        bcMatrix->setFluid(ix1, ix2, ix3);
                    else
                        bcMatrix->setUndefined(ix1, ix2, ix3);
                }
            }
        }

        return;
    }
}
//////////////////////////////////////////////////////////////////////////
void SetUndefinedNodesBlockVisitor::setNodesUndefined(int startix1, int endix1, int startix2, int endix2, int startix3,
                                                      int endix3, SPtr<BCArray3D> bcMatrix)
{
    for (int ix3 = startix3; ix3 <= endix3; ix3++)
        for (int ix2 = startix2; ix2 <= endix2; ix2++)
            for (int ix1 = startix1; ix1 <= endix1; ix1++) {
                bcMatrix->setUndefined(ix1, ix2, ix3);
            }
}

//! \}
