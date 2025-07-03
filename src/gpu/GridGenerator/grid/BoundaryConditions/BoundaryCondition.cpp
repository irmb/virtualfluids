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
//=======================================================================================
#include "BoundaryCondition.h"

#include <cmath>
#include <stdexcept>

#include "DataTypes.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/Grid.h"
#include "GridGenerator/TransientBCSetter/TransientBCSetter.h"

bool grid_generator::BoundaryCondition::isSide( SideType side ) const
{
    return this->side->whoAmI() == side;
}
//////////////////////////////////////////////////////////////////////////

void VelocityBoundaryCondition::setVelocityProfile(
    SPtr<Grid> grid, std::function<void(real, real, real, real &, real &, real &)> velocityProfile)
{
    for (uint index = 0; index < this->indices.size(); index++) {

        real x, y, z;

        grid->transIndexToCoords(this->indices[index], x, y, z);

        velocityProfile(x, y, z, this->vxList[index], this->vyList[index], this->vzList[index]);
    }
}

//////////////////////////////////////////////////////////////////////////

void GeometryBoundaryCondition::setTangentialVelocityForPatch(SPtr<Grid> grid, uint patch, 
                                                              real p1x, real p1y, real p1z, 
                                                              real p2x, real p2y, real p2z, 
                                                              real v, real r)
{
    for( uint index = 0; index < this->indices.size(); index++ ){
        if( this->patches[index] == patch ){

            real x, y, z;

            grid->transIndexToCoords( this->indices[index], x, y, z );

            real pVecX = p2x - p1x;
            real pVecY = p2y - p1y;
            real pVecZ = p2z - p1z;

            real xVecX =   x - p1x;
            real xVecY =   y - p1y;
            real xVecZ =   z - p1z;

            // compute and normalize tangent

            real tangentX = pVecY * xVecZ - pVecZ * xVecY;
            real tangentY = pVecZ * xVecX - pVecX * xVecZ; 
            real tangentZ = pVecX * xVecY - pVecY * xVecX;

            real tangentNorm = sqrt( tangentX*tangentX + tangentY*tangentY + tangentZ*tangentZ );

            tangentX /= tangentNorm;
            tangentY /= tangentNorm;
            tangentZ /= tangentNorm;

            // compute distance from rotation axis

            real projection = ( pVecX*xVecX + pVecY*xVecY + pVecZ*xVecZ ) / ( pVecX*pVecX + pVecY*pVecY + pVecZ*pVecZ );

            real d = sqrt( ( xVecX - projection * pVecX ) * ( xVecX - projection * pVecX )
                         + ( xVecY - projection * pVecY ) * ( xVecY - projection * pVecY )
                         + ( xVecZ - projection * pVecZ ) * ( xVecZ - projection * pVecZ ) );

            this->vxList[index] = tangentX * d / r * v;
            this->vyList[index] = tangentY * d / r * v;
            this->vzList[index] = tangentZ * d / r * v;

        }
    }
}

//////////////////////////////////////////////////////////////////////////

void StressBoundaryCondition::setSamplingIndices(const SPtr<Grid>& grid, const uint index)
{
    using namespace vf::lbm::dir;
    real nodeX, nodeY, nodeZ;
    grid->transIndexToCoords(index, nodeX, nodeY, nodeZ);
    const auto normal = side->getNormal();
    const real samplingOffsetFromNode = static_cast<real>(getSamplingOffset() - 1);
    const real samplingX = nodeX - normal[0] * grid->getDelta() * samplingOffsetFromNode;
    const real samplingY = nodeY - normal[1] * grid->getDelta() * samplingOffsetFromNode;
    const real samplingZ = nodeZ - normal[2] * grid->getDelta() * samplingOffsetFromNode;
    const uint samplingIndex = grid->transCoordToIndex(samplingX, samplingY, samplingZ);
    if(samplingIndex == INVALID_INDEX)
        throw std::runtime_error("StressBoundaryCondition::setSamplingIndices could not find sampling index");
    samplingIndices.push_back(samplingIndex);

    const auto qNode = qs.back();
    real distanceNodeToWall;
    switch (side->whoAmI()) {
        case SideType::PX:
            distanceNodeToWall = qNode[dP00];
            break;
        case SideType::MX:
            distanceNodeToWall = qNode[dM00];
            break;
        case SideType::PY:
            distanceNodeToWall = qNode[d0P0];
            break;
        case SideType::MY:
            distanceNodeToWall = qNode[d0M0];
            break;
        case SideType::PZ:
            distanceNodeToWall = qNode[d00P];
            break;
        case SideType::MZ:
            distanceNodeToWall = qNode[d00M];
            break;
        default:
            throw std::runtime_error("Side type not implemented!");
    }

    const real samplingDistance =
        std::hypot(samplingX - nodeX, samplingY - nodeY, samplingZ - nodeZ) / grid->getDelta() + distanceNodeToWall;
    samplingDistanceList.push_back(samplingDistance);
}

void PressureBoundaryCondition::setNeighborIndices(const SPtr<Grid> &grid, uint index)
{

    real x, y, z;
    grid->transIndexToCoords(index, x, y, z);

    real nx = x;
    real ny = y;
    real nz = z;

    if (side->getCoordinate() == X_INDEX)
        nx = -side->getDirection() * grid->getDelta() + x;
    if (side->getCoordinate() == Y_INDEX)
        ny = -side->getDirection() * grid->getDelta() + y;
    if (side->getCoordinate() == Z_INDEX)
        nz = -side->getDirection() * grid->getDelta() + z;

    uint neighborIndex = grid->transCoordToIndex(nx, ny, nz);
    neighborIndices.push_back(neighborIndex);
}

void ADNeumannBoundaryCondition::fillBoundaryValueLists()
{
    const real gradient = static_cast<real>(side->getDirection()) * this->gradient;
    std::fill_n(std::back_inserter(this->vxList), this->indices.size(), vx);
    std::fill_n(std::back_inserter(this->vyList), this->indices.size(), vy);
    std::fill_n(std::back_inserter(this->vzList), this->indices.size(), vz);
    std::fill_n(std::back_inserter(this->gradientList), this->indices.size(), gradient);
}

void ADFluxBoundaryCondition::fillBoundaryValueLists()
{
    const real gradient = static_cast<real>(side->getDirection()) * this->gradient;
    std::fill_n(std::back_inserter(this->normalXList), this->indices.size(), normalX);
    std::fill_n(std::back_inserter(this->normalYList), this->indices.size(), normalY);
    std::fill_n(std::back_inserter(this->normalZList), this->indices.size(), normalZ);
    std::fill_n(std::back_inserter(this->gradientList), this->indices.size(), gradient);
}
//! \}
