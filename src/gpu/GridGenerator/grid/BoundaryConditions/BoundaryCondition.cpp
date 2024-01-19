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

#include "grid/BoundaryConditions/Side.h"
#include "grid/Grid.h"
#include "GridGenerator/TransientBCSetter/TransientBCSetter.h"

bool gg::BoundaryCondition::isSide( SideType side ) const
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

void StressBoundaryCondition::fillSamplingIndices(std::vector<SPtr<Grid> > grid, uint level, uint samplingOffset)
{

    for( uint i = 0; i < this->indices.size(); i++ )
    {
        real x, y, z;
        grid[level]->transIndexToCoords(this->indices[i], x, y, z);

        real x_sampling = x + this->getNormalx(i)*samplingOffset*grid[level]->getDelta();
        real y_sampling = y + this->getNormaly(i)*samplingOffset*grid[level]->getDelta();
        real z_sampling = z + this->getNormalz(i)*samplingOffset*grid[level]->getDelta();

        this->velocitySamplingIndices.push_back( grid[level]->transCoordToIndex(x_sampling, y_sampling, z_sampling) );
    }
    
}
//! \}
