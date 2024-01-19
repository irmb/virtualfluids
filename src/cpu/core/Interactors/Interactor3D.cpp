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
//! \addtogroup cpu_Interactors Interactors
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================

#include "Interactor3D.h"

#include "UbException.h"
//#include <basics/utilities/UbMath.h>
#include <fstream>
#include <geometry3d/GbCuboid3D.h>

#include "Block3D.h"
#include "GbObject3D.h"
#include "Grid3D.h"

using namespace std;

const int Interactor3D::SOLID              = (1 << 0); // 1
const int Interactor3D::INVERSESOLID       = (1 << 1); // 2
const int Interactor3D::TIMEDEPENDENT      = (1 << 2); // 4   //zeitlich
const int Interactor3D::FLUID              = (1 << 3); // 8
const int Interactor3D::MOVEABLE           = (1 << 4); // 16  // geometrisch
const int Interactor3D::CHANGENOTNECESSARY = (1 << 5); // 32

//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D() : type(SOLID) {}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(SPtr<Grid3D> grid, int type) : grid(grid), type(type) {}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type)
    : geoObject3D(geoObject3D), grid(grid), type(type), accuracy(SIMPLE)
{
}
//////////////////////////////////////////////////////////////////////////
Interactor3D::Interactor3D(SPtr<GbObject3D> geoObject3D, SPtr<Grid3D> grid, int type, Interactor3D::Accuracy a)
    : geoObject3D(geoObject3D), grid(grid), type(type), accuracy(a)
{
}
//////////////////////////////////////////////////////////////////////////
Interactor3D::~Interactor3D() = default;
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::arePointsInsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2,
                                            real maxX3, real delta)
{
    bool result = true;
    for (real ix3 = minX3; ix3 <= maxX3; ix3 += delta)
        for (real ix2 = minX2; ix2 <= maxX2; ix2 += delta)
            for (real ix1 = minX1; ix1 <= maxX1; ix1 += delta)
                result = result && this->geoObject3D->isPointInGbObject3D(ix1, ix2, ix3);

    return result;
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::arePointsOutsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2,
                                             real maxX3, real delta)
{
    bool result = true;
    for (real ix3 = minX3; ix3 <= maxX3; ix3 += delta)
        for (real ix2 = minX2; ix2 <= maxX2; ix2 += delta)
            for (real ix1 = minX1; ix1 <= maxX1; ix1 += delta)
                result = result && (!this->geoObject3D->isPointInGbObject3D(ix1, ix2, ix3));

    return result;
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::arePointsCuttingGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2,
                                             real maxX3, real delta)
{
    bool result = true;
    for (real ix3 = minX3; ix3 <= maxX3; ix3 += delta)
        for (real ix2 = minX2; ix2 <= maxX2; ix2 += delta)
            for (real ix1 = minX1; ix1 <= maxX1; ix1 += delta)
                result = result || this->geoObject3D->isPointInGbObject3D(ix1, ix2, ix3);

    return result;
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isBlockOutsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2,
                                           real maxX3, real delta)
{
    switch (accuracy) {
            // simple duff
        case SIMPLE:
            return !this->geoObject3D->isCellInsideOrCuttingGbObject3D(minX1, minX2, minX3, maxX1, maxX2, maxX3);
            // test only edges
        case EDGES:
            return arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, delta) &&
                   arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, delta) &&
                   arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, delta) &&

                   arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, delta) &&
                   arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
                   arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, delta) &&

                   arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, delta);
            // test only faces
        case FACES:
            return arePointsOutsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, delta) &&
                   arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
                   arePointsOutsideGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, delta);
            // test all points
        case POINTS:
            return arePointsOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, delta);
        default:
            UB_THROW(UbException(UB_EXARGS, "Accuracy isn't correct"));
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isBlockInsideGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2,
                                          real maxX3, real delta)
{
    switch (accuracy) {
            // simple duff
        case SIMPLE:
            return this->geoObject3D->isCellInsideGbObject3D(minX1, minX2, minX3, maxX1, maxX2, maxX3);
            // test only edges
        case EDGES:
            return arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, delta) &&
                   arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, delta) &&
                   arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
                   arePointsInsideGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, delta) &&

                   arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, delta) &&
                   arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
                   arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) &&
                   arePointsInsideGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, delta) &&

                   arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
                   arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
                   arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, delta) &&
                   arePointsInsideGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, delta);
            // test only faces
        case FACES:
            return arePointsInsideGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) &&
                   arePointsInsideGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) &&
                   arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, delta) &&
                   arePointsInsideGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, delta) &&
                   arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, delta) &&
                   arePointsInsideGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, delta);
            // test all points
        case POINTS:
            return arePointsInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, delta);
        default:
            UB_THROW(UbException(UB_EXARGS, "Accuracy isn't correct"));
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isBlockCuttingGeoObject(real minX1, real minX2, real minX3, real maxX1, real maxX2,
                                           real maxX3, real delta)
{
    switch (accuracy) {
            // simple duff
        case SIMPLE:
            return this->geoObject3D->isCellCuttingGbObject3D(minX1, minX2, minX3, maxX1, maxX2, maxX3);
            // test only edges
        case EDGES:
            return arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, minX2, minX3, delta) ||
                   arePointsCuttingGeoObject(minX1, maxX2, minX3, maxX1, maxX2, minX3, delta) ||
                   arePointsCuttingGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(minX1, maxX2, maxX3, maxX1, maxX2, maxX3, delta) ||

                   arePointsCuttingGeoObject(minX1, minX2, minX3, minX1, maxX2, minX3, delta) ||
                   arePointsCuttingGeoObject(maxX1, minX2, minX3, maxX1, maxX2, minX3, delta) ||
                   arePointsCuttingGeoObject(minX1, minX2, maxX3, maxX1, minX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(maxX1, minX2, maxX3, maxX1, maxX2, maxX3, delta) ||

                   arePointsCuttingGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(minX1, maxX2, minX3, maxX1, minX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(maxX1, maxX2, minX3, maxX1, maxX2, maxX3, delta);
            // test only faceCutting
        case FACES:
            return arePointsCuttingGeoObject(minX1, minX2, minX3, minX1, maxX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(maxX1, minX2, minX3, maxX1, maxX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, minX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(minX1, maxX2, minX3, maxX1, maxX2, maxX3, delta) ||
                   arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, maxX2, minX3, delta) ||
                   arePointsCuttingGeoObject(minX1, minX2, maxX3, maxX1, maxX2, maxX3, delta);
            // test all pointCutting
        case POINTS:
            return arePointsCuttingGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, delta);
        default:
            UB_THROW(UbException(UB_EXARGS, "Accuracy isn't correct"));
            break;
    }
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setSolidBlock(SPtr<Block3D> block)
{
    real minX1, minX2, minX3, maxX1, maxX2, maxX3;

    real deltaX               = grid.lock()->getDeltaX(block);
    UbTupleDouble3 blockLengths = grid.lock()->getBlockLengths(block);
    UbTupleDouble3 org          = grid.lock()->getBlockWorldCoordinates(block);
    UbTupleDouble3 nodeOffset   = grid.lock()->getNodeOffset(block);

    // coordinates of block without ghost layer
    minX1 = val<1>(org) + val<1>(nodeOffset);
    minX2 = val<2>(org) + val<2>(nodeOffset);
    minX3 = val<3>(org) + val<3>(nodeOffset);
    maxX1 = val<1>(org) + val<1>(blockLengths) - val<1>(nodeOffset);
    maxX2 = val<2>(org) + val<2>(blockLengths) - val<2>(nodeOffset);
    maxX3 = val<3>(org) + val<3>(blockLengths) - val<3>(nodeOffset);

    if (this->isInverseSolid()) {
        if (isBlockOutsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX)) {
            block->setActive(false);
            this->solidBlocks.push_back(block);
        }
    } else // solid
    {
        if (isBlockInsideGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX)) {
            block->setActive(false);
            this->solidBlocks.push_back(block);
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setBCBlock(SPtr<Block3D> block)
{
    real minX1, minX2, minX3, maxX1, maxX2, maxX3;

    real deltaX               = grid.lock()->getDeltaX(block);
    UbTupleDouble3 blockLengths = grid.lock()->getBlockLengths(block);
    UbTupleDouble3 org          = grid.lock()->getBlockWorldCoordinates(block);
    UbTupleDouble3 nodeOffset   = grid.lock()->getNodeOffset(block);

    // coordinates of block with ghost layer
    minX1 = val<1>(org) - val<1>(nodeOffset);
    minX2 = val<2>(org) - val<2>(nodeOffset);
    minX3 = val<3>(org) - val<3>(nodeOffset);
    maxX1 = val<1>(org) + val<1>(blockLengths) + val<1>(nodeOffset);
    maxX2 = val<2>(org) + val<2>(blockLengths) + val<2>(nodeOffset);
    maxX3 = val<3>(org) + val<3>(blockLengths) + val<3>(nodeOffset);

    if (isBlockCuttingGeoObject(minX1, minX2, minX3, maxX1, maxX2, maxX3, deltaX))
        this->bcBlocks.push_back(block);
}

UbTupleDouble3 Interactor3D::getForces()
{
    UB_THROW(UbException("UbTupleDouble3 getForces() - gehoert in die abgeleitete klasse"));
}
void Interactor3D::setID(int id) { this->id = id; }
//////////////////////////////////////////////////////////////////////////
int Interactor3D::getID() { return id; }
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setActive() { active = true; }
//////////////////////////////////////////////////////////////////////////
void Interactor3D::setInactive() { active = false; }
//////////////////////////////////////////////////////////////////////////
bool Interactor3D::isActive() { return active; }
//////////////////////////////////////////////////////////////////////////
void Interactor3D::updateBlocks()
{
    for (SPtr<Block3D> block : bcBlocks) 
    {
        this->setDifferencesToGbObject3D(block);
    }
}
//////////////////////////////////////////////////////////////////////////
void Interactor3D::updateInteractor(const real & /*timeStep*/)
{
    UB_THROW(UbException("Interactor3D::updateInteractor - toDo"));
}
//////////////////////////////////////////////////////////////////////////

//! \}
