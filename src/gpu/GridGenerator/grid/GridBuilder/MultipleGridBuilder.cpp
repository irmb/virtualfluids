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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file MultipleGridBuilder.cpp
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "MultipleGridBuilder.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>

#include "utilities/math/Math.h"

#include "geometries/BoundingBox/BoundingBox.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/Grid.h"
#include "grid/GridFactory.h"

MultipleGridBuilder::MultipleGridBuilder(SPtr<GridFactory> gridFactory, Device device, const std::string &d3qxx) :
    LevelGridBuilder(device, d3qxx), gridFactory(gridFactory), solidObject(nullptr), numberOfLayersFine(12), numberOfLayersBetweenLevels(8), subDomainBox(nullptr)
{

}

SPtr<MultipleGridBuilder> MultipleGridBuilder::makeShared(SPtr<GridFactory> gridFactory)
{
    return SPtr<MultipleGridBuilder>(new MultipleGridBuilder(gridFactory));
}

void MultipleGridBuilder::addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta)
{
    boundaryConditions.push_back(SPtr<BoundaryConditions>(new BoundaryConditions));

    startX -= 0.5 * delta;
    startY -= 0.5 * delta;
    startZ -= 0.5 * delta;
    endX   += 0.5 * delta;
    endY   += 0.5 * delta;
    endZ   += 0.5 * delta;

    const auto grid = this->makeGrid(new Cuboid(startX, startY, startZ, endX, endY, endZ), startX, startY, startZ, endX, endY, endZ, delta, 0);
    addGridToList(grid);
}

SPtr<Grid> MultipleGridBuilder::makeGrid(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const
{
    return gridFactory->makeGrid(gridShape, startX, startY, startZ, endX, endY, endZ, delta, level);
}

bool MultipleGridBuilder::coarseGridExists() const
{
    return !grids.empty();
}

uint MultipleGridBuilder::getNumberOfLevels() const
{
    return uint(grids.size());
}

real MultipleGridBuilder::getDelta(uint level) const
{
    if (grids.size() <= level)
        throw std::runtime_error("delta from invalid level was required.");
    return grids[level]->getDelta();
}

real MultipleGridBuilder::getStartX(uint level) const
{
    return grids[level]->getStartX();
}

real MultipleGridBuilder::getStartY(uint level) const
{
    return grids[level]->getStartY();
}

real MultipleGridBuilder::getStartZ(uint level) const
{
    return grids[level]->getStartZ();
}

real MultipleGridBuilder::getEndX(uint level) const
{
    return grids[level]->getEndX();
}

real MultipleGridBuilder::getEndY(uint level) const
{
    return grids[level]->getEndY();
}

real MultipleGridBuilder::getEndZ(uint level) const
{
    return grids[level]->getEndZ();
}

void MultipleGridBuilder::addGridToList(SPtr<Grid> grid)
{
    grids.push_back(grid);
}

std::vector<SPtr<Grid> > MultipleGridBuilder::getGrids() const
{
    return this->grids;
}

void MultipleGridBuilder::buildGrids( LbmOrGks lbmOrGks, bool enableThinWalls )
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start initializing level " << 0 << "\n";

    grids[0]->inital( nullptr, 0 );

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Done initializing level " << 0 << "\n";

    //////////////////////////////////////////////////////////////////////////

	if (lbmOrGks == LBM)
		grids[0]->findSparseIndices(nullptr);
}

GRIDGENERATOR_EXPORT void MultipleGridBuilder::setNumberOfLayers(uint numberOfLayersFine, uint numberOfLayersBetweenLevels)
{
    this->numberOfLayersFine = numberOfLayersFine;
    this->numberOfLayersBetweenLevels = numberOfLayersBetweenLevels;
}

