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

#include "GenBlocksGridVisitor.h"
#include "Block3D.h"
#include "CoordinateTransformation3D.h"
#include "Grid3D.h"
#include "D3Q27System.h"

#include <geometry3d/GbObject3D.h>

GenBlocksGridVisitor::GenBlocksGridVisitor(SPtr<GbObject3D> boundingBox) : boundingBox(boundingBox) {}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::visit(const SPtr<Grid3D> grid)
{
    real orgX1 = boundingBox->getX1Minimum();
    real orgX2 = boundingBox->getX2Minimum();
    real orgX3 = boundingBox->getX3Minimum();

    real dx = grid->getDeltaX(0);

    UbTupleInt3 blockNX = grid->getBlockNX();

    real blockLentghX1 = (real)val<1>(blockNX) * dx;
    real blockLentghX2 = (real)val<2>(blockNX) * dx;
    real blockLentghX3 = (real)val<3>(blockNX) * dx;

    SPtr<CoordinateTransformation3D> trafo(
        new CoordinateTransformation3D(orgX1, orgX2, orgX3, blockLentghX1, blockLentghX2, blockLentghX3));
    grid->setCoordinateTransformator(trafo);

    genBlocks(grid);
}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::fillExtentWithBlocks(SPtr<Grid3D> grid)
{
    for (int x3 = val<3>(minInd); x3 < val<3>(maxInd); x3++) {
        for (int x2 = val<2>(minInd); x2 < val<2>(maxInd); x2++) {
            for (int x1 = val<1>(minInd); x1 < val<1>(maxInd); x1++) {
                SPtr<Block3D> block(new Block3D(x1, x2, x3, 0));
                grid->addBlock(block);
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void GenBlocksGridVisitor::genBlocks(SPtr<Grid3D> grid)
{
    minInd =
        grid->getBlockIndexes(boundingBox->getX1Minimum(), boundingBox->getX2Minimum(), boundingBox->getX3Minimum());
    real geoMaxX1           = boundingBox->getX1Maximum();
    real geoMaxX2           = boundingBox->getX2Maximum();
    real geoMaxX3           = boundingBox->getX3Maximum();
    maxInd                    = grid->getBlockIndexes(geoMaxX1, geoMaxX2, geoMaxX3);
    UbTupleDouble3 blockCoord = grid->getBlockWorldCoordinates(
        static_cast<int>(val<1>(maxInd)), static_cast<int>(val<2>(maxInd)), static_cast<int>(val<3>(maxInd)), 0);
    // if (geoMaxX1 > val<1>(blockCoord))
    //    val<1>(maxInd) += 1;
    // if (geoMaxX2 > val<2>(blockCoord))
    //    val<2>(maxInd) += 1;
    // if (geoMaxX3 > val<3>(blockCoord))
    //    val<3>(maxInd) += 1;

    real dx = grid->getDeltaX(0);
    if (fabs(geoMaxX1 - val<1>(blockCoord)) > dx)
        val<1>(maxInd) += 1;
    if (fabs(geoMaxX2 - val<2>(blockCoord)) > dx)
        val<2>(maxInd) += 1;
    if (fabs(geoMaxX3 - val<3>(blockCoord)) > dx)
        val<3>(maxInd) += 1;

    this->fillExtentWithBlocks(grid);

    grid->setNX1(val<1>(maxInd));
    grid->setNX2(val<2>(maxInd));
    grid->setNX3(val<3>(maxInd));
}

//! \}
