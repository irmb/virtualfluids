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
//! \addtogroup gpu_grid grid
//! \ingroup gpu_GridGenerator GridGenerator
//! \{
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "Side.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/Grid.h"
#include "grid/NodeValues.h"

#include "utilities/math/Math.h"
#include <array>
#include <cstddef>
#include <vector>

using namespace gg;

std::array<real, 3> Side::getNormal() const
{
    std::array<real, 3> normal;
    if(this->getCoordinate()==X_INDEX)
        normal = {(real)this->getDirection(), 0.0, 0.0};
    if(this->getCoordinate()==Y_INDEX)
        normal = {0.0, (real)this->getDirection(), 0.0};
    if(this->getCoordinate()==Z_INDEX)
        normal = {0.0, 0.0, (real)this->getDirection()};
    return normal;
}

void Side::addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::string coord, real constant,
                      real startInner, real endInner, real startOuter, real endOuter)
{
    for (real v2 = startOuter; v2 <= endOuter; v2 += grid->getDelta())
    {
        for (real v1 = startInner; v1 <= endInner; v1 += grid->getDelta())
        {
            const uint index = getIndex(grid, coord, constant, v1, v2);

            if(index == INVALID_INDEX)
                continue;

            if (   grid->getFieldEntry(index) == vf::gpu::FLUID
                                            ||  grid->getFieldEntry(index) == vf::gpu::FLUID_CFC
                                            ||  grid->getFieldEntry(index) == vf::gpu::FLUID_CFF
                                            ||  grid->getFieldEntry(index) == vf::gpu::FLUID_FCC
                                            ||  grid->getFieldEntry(index) == vf::gpu::FLUID_FCF
                                            ||  grid->getFieldEntry(index) == vf::gpu::FLUID_FCF
                                            // Overlap of BCs on edge nodes
                                            || grid->nodeHasBC(index) )
            {   
                grid->setFieldEntry(index, boundaryCondition->getType());
                boundaryCondition->indices.push_back(index);
                setPressureNeighborIndices(boundaryCondition, grid, index);
                setStressSamplingIndices(boundaryCondition, grid, index);
                // if(grid->getFieldEntry(index)==26) printf("index = %u, v1 = %f, v2 = %f, field entry=%u \n", index, v1, v2, grid->getFieldEntry(index) );
                setQs(grid, boundaryCondition, index);
                boundaryCondition->patches.push_back(0);
            }
        }
    }

    const auto currentBCSide = this->whoAmI();
    if(currentBCSide != SideType::GEOMETRY)
        grid->addBCalreadySet(currentBCSide);
}

void Side::setPressureNeighborIndices(SPtr<BoundaryCondition> boundaryCondition, SPtr<Grid> grid, const uint index)
{
    auto pressureBoundaryCondition = std::dynamic_pointer_cast<PressureBoundaryCondition>(boundaryCondition);
    if (pressureBoundaryCondition)
    {
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);

        real nx = x;
        real ny = y;
        real nz = z;

        if (boundaryCondition->side->getCoordinate() == X_INDEX)
            nx = -boundaryCondition->side->getDirection() * grid->getDelta() + x;
        if (boundaryCondition->side->getCoordinate() == Y_INDEX)
            ny = -boundaryCondition->side->getDirection() * grid->getDelta() + y;
        if (boundaryCondition->side->getCoordinate() == Z_INDEX)
            nz = -boundaryCondition->side->getDirection() * grid->getDelta() + z;

        int neighborIndex = grid->transCoordToIndex(nx, ny, nz);
        pressureBoundaryCondition->neighborIndices.push_back(neighborIndex);
    }
}

void Side::setStressSamplingIndices(SPtr<BoundaryCondition> boundaryCondition, SPtr<Grid> grid, const uint index)
{
    auto stressBoundaryCondition = std::dynamic_pointer_cast<StressBoundaryCondition>(boundaryCondition);
    if (stressBoundaryCondition)
    {
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);

        real nx = x;
        real ny = y;
        real nz = z;

        if (boundaryCondition->side->getCoordinate() == X_INDEX)
            nx = -boundaryCondition->side->getDirection() * stressBoundaryCondition->samplingOffset * grid->getDelta() + x;
        if (boundaryCondition->side->getCoordinate() == Y_INDEX)
            ny = -boundaryCondition->side->getDirection() * stressBoundaryCondition->samplingOffset * grid->getDelta() + y;
        if (boundaryCondition->side->getCoordinate() == Z_INDEX)
            nz = -boundaryCondition->side->getDirection() * stressBoundaryCondition->samplingOffset * grid->getDelta() + z;

        uint samplingIndex = grid->transCoordToIndex(nx, ny, nz);
        stressBoundaryCondition->velocitySamplingIndices.push_back(samplingIndex);
    }
}

void Side::setQs(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, uint index)
{
    std::vector<real> qNode(grid->getEndDirection() + 1);

    for (int dir = 0; dir <= grid->getEndDirection(); dir++) {
        real x, y, z;
        grid->transIndexToCoords(index, x, y, z);

        std::array<real, 3> coords = { x, y, z };
        std::array<real, 3> neighborCoords = getNeighborCoordinates(grid.get(), coords, (size_t)dir);

        correctNeighborForPeriodicBoundaries(grid.get(), coords, neighborCoords);

        const uint neighborIndex = grid->transCoordToIndex(neighborCoords[0], neighborCoords[1], neighborCoords[2]);

        //! Only setting q's that partially point in the Side-normal direction
        const bool alignedWithNormal = this->isAlignedWithMyNormal(grid.get(), dir);
        if (grid->isStopperForBC(neighborIndex) && alignedWithNormal) {
            qNode[dir] = 0.5;
        } else {
            qNode[dir] = -1.0;
        }

        // reset diagonals in case they were set by another bc
        resetDiagonalsInCaseOfOtherBC(grid.get(), qNode, dir, coords);
    }

    boundaryCondition->qs.push_back(qNode);
}

std::array<real, 3> Side::getNeighborCoordinates(Grid *grid, const std::array<real, 3> &coordinates, size_t direction) const
{
    return { coordinates[0] + grid->getDirection()[direction * DIMENSION + 0] * grid->getDelta(),
             coordinates[1] + grid->getDirection()[direction * DIMENSION + 1] * grid->getDelta(),
             coordinates[2] + grid->getDirection()[direction * DIMENSION + 2] * grid->getDelta() };
}

bool Side::neighborNormalToSideIsAStopper(Grid *grid, const std::array<real, 3> &coordinates, SideType side) const
{
    const auto neighborCoords = getNeighborCoordinates(grid, coordinates, sideToD3Q27.at(side));
    const auto neighborIndex = grid->transCoordToIndex(neighborCoords[0], neighborCoords[1], neighborCoords[2]);
    return grid->isStopperForBC(neighborIndex);
}

void Side::resetDiagonalsInCaseOfOtherBC(Grid *grid, std::vector<real> &qNode, int dir,
                                         const std::array<real, 3> &coordinates) const
{
    // When to reset a diagonal q to -1:
    // - it is normal to another boundary condition which was already set
    // - and it actually is influenced by the other bc:
    //   We check if its neighbor in the regular direction to the other bc is a stopper. If it is a stopper, it is influenced by the other bc.

    if (qNode[dir] == 0.5 && grid->getBCAlreadySet().size() > 0) {
        for (int i = 0; i < (int)grid->getBCAlreadySet().size(); i++) {
            SideType otherDir = grid->getBCAlreadySet()[i];

            // only reset normals for nodes on edges and corners, not on faces
            if (!neighborNormalToSideIsAStopper(grid, coordinates, otherDir))
                continue;

            const auto otherNormal = normals.at(otherDir);
            if (isAlignedWithNormal(grid, dir, otherNormal)) {
                qNode[dir] = -1.0;
            }
        }
    }
}

bool Side::isAlignedWithMyNormal(const Grid *grid, int dir) const
{
    std::array<real, 3> normal = this->getNormal();
    return isAlignedWithNormal(grid, dir, normal);
}

bool Side::isAlignedWithNormal(const Grid *grid, int dir, const std::array<real, 3> &normal) const
{
    return (normal[0] * grid->getDirection()[dir * DIMENSION + 0] +
            normal[1] * grid->getDirection()[dir * DIMENSION + 1] +
            normal[2] * grid->getDirection()[dir * DIMENSION + 2]) > 0;
}

void Side::correctNeighborForPeriodicBoundaries(const Grid *grid, std::array<real, 3>& coords, std::array<real, 3>& neighborCoords) const
{
    // correct neighbor coordinates in case of periodic boundaries
    if (grid->getPeriodicityX() &&
        grid->getFieldEntry(grid->transCoordToIndex(neighborCoords[0], coords[1], coords[2])) == vf::gpu::STOPPER_OUT_OF_GRID_BOUNDARY) {
        if (neighborCoords[0] > coords[0])
            neighborCoords[0] = grid->getFirstFluidNode(coords.data(), 0, grid->getStartX());
        else
            neighborCoords[0] = grid->getLastFluidNode(coords.data(), 0, grid->getEndX());
    }

    if (grid->getPeriodicityY() &&
        grid->getFieldEntry(grid->transCoordToIndex(coords[0], neighborCoords[1], coords[2])) == vf::gpu::STOPPER_OUT_OF_GRID_BOUNDARY) {
        if (neighborCoords[1] > coords[1])
            neighborCoords[1] = grid->getFirstFluidNode(coords.data(), 1, grid->getStartY());
        else
            neighborCoords[1] = grid->getLastFluidNode(coords.data(), 1, grid->getEndY());
    }

    if (grid->getPeriodicityZ() &&
        grid->getFieldEntry(grid->transCoordToIndex(coords[0], coords[1], neighborCoords[2])) == vf::gpu::STOPPER_OUT_OF_GRID_BOUNDARY) {
        if (neighborCoords[2] > coords[2])
            neighborCoords[2] = grid->getFirstFluidNode(coords.data(), 2, grid->getStartZ());
        else
            neighborCoords[2] = grid->getLastFluidNode(coords.data(), 2, grid->getEndZ());
    }
}

uint Side::getIndex(SPtr<Grid> grid, std::string coord, real constant, real v1, real v2)
{
    if (coord == "x")
        return grid->transCoordToIndex(constant, v1, v2);
    if (coord == "y")
        return grid->transCoordToIndex(v1, constant, v2);
    if (coord == "z")
        return grid->transCoordToIndex(v1, v2, constant);
    return INVALID_INDEX;
}


void Geometry::addIndices(const std::vector<SPtr<Grid>> &grids, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    auto geometryBoundaryCondition = std::dynamic_pointer_cast<GeometryBoundaryCondition>(boundaryCondition);

    std::vector<real> qNode(grids[level]->getEndDirection() + 1);

    for (uint index = 0; index < grids[level]->getSize(); index++)
    {
        if (grids[level]->getFieldEntry(index) != vf::gpu::BC_SOLID)
            continue;

        for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
        {
            const real q = grids[level]->getQValue(index, dir);

            qNode[dir] = q;

            // also the neighbor if any Qs are required
            real x,y,z;
            grids[level]->transIndexToCoords( index, x, y, z );

            x += grids[level]->getDirection()[dir * DIMENSION + 0] * grids[level]->getDelta();
            y += grids[level]->getDirection()[dir * DIMENSION + 1] * grids[level]->getDelta();
            z += grids[level]->getDirection()[dir * DIMENSION + 2] * grids[level]->getDelta();

            uint neighborIndex = grids[level]->transCoordToIndex( x, y, z );

            if( qNode[dir] < -0.5 && ( grids[level]->getFieldEntry(neighborIndex) == vf::gpu::STOPPER_OUT_OF_GRID_BOUNDARY ||
                                       grids[level]->getFieldEntry(neighborIndex) == vf::gpu::STOPPER_OUT_OF_GRID ||
                                       grids[level]->getFieldEntry(neighborIndex) == vf::gpu::STOPPER_SOLID ) )
                qNode[dir] = 0.5;
        }

        geometryBoundaryCondition->indices.push_back(index);
        geometryBoundaryCondition->qs.push_back(qNode);
        geometryBoundaryCondition->patches.push_back(grids[level]->getQPatch(index));
    }
}



void MX::addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartY();
    real endInner = grid[level]->getEndY();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getStartX() + grid[level]->getDelta();

    if( coordinateNormal > grid[0]->getStartX() + grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "x", coordinateNormal, startInner, endInner, startOuter, endOuter);

}

void PX::addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartY();
    real endInner = grid[level]->getEndY();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getEndX() - grid[level]->getDelta();

    if( coordinateNormal < grid[0]->getEndX() - grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "x", coordinateNormal, startInner, endInner, startOuter, endOuter);
}

void MY::addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getStartY() + grid[level]->getDelta();

    if( coordinateNormal > grid[0]->getStartY() + grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "y", coordinateNormal, startInner, endInner, startOuter, endOuter);
}


void PY::addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartZ();
    real endOuter = grid[level]->getEndZ();

    real coordinateNormal = grid[level]->getEndY() - grid[level]->getDelta();

    if( coordinateNormal < grid[0]->getEndY() - grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "y", coordinateNormal, startInner, endInner, startOuter, endOuter);
}


void MZ::addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartY();
    real endOuter = grid[level]->getEndY();

    real coordinateNormal = grid[level]->getStartZ() + grid[level]->getDelta();

    if( coordinateNormal > grid[0]->getStartZ() + grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "z", coordinateNormal, startInner, endInner, startOuter, endOuter);
}

void PZ::addIndices(const std::vector<SPtr<Grid>> &grid, uint level, SPtr<BoundaryCondition> boundaryCondition)
{
    real startInner = grid[level]->getStartX();
    real endInner = grid[level]->getEndX();

    real startOuter = grid[level]->getStartY();
    real endOuter = grid[level]->getEndY();

    real coordinateNormal = grid[level]->getEndZ() - grid[level]->getDelta();

    if( coordinateNormal < grid[0]->getEndZ() - grid[0]->getDelta() ) return;

    Side::addIndices(grid[level], boundaryCondition, "z", coordinateNormal, startInner, endInner, startOuter, endOuter);
}

//! \}
