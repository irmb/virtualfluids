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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file MultipleGridBuilder.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef MULTIPLE_GRID_BUILDER_H
#define MULTIPLE_GRID_BUILDER_H

#include <vector>
#include <array>

#include "Core/LbmOrGks.h"

#include "global.h"

#include "grid/GridBuilder/LevelGridBuilder.h"
#include "grid/GridFactory.h"

class Object;
class BoundingBox;

class MultipleGridBuilder : public LevelGridBuilder
{
private:
    GRIDGENERATOR_EXPORT MultipleGridBuilder(SPtr<GridFactory> gridFactory, Device device = Device::CPU, const std::string &d3qxx = "D3Q27");

public:
    GRIDGENERATOR_EXPORT static SPtr<MultipleGridBuilder> makeShared(SPtr<GridFactory> gridFactory);

    GRIDGENERATOR_EXPORT void addCoarseGrid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta);

    GRIDGENERATOR_EXPORT uint getNumberOfLevels() const;
    GRIDGENERATOR_EXPORT real getDelta(uint level) const;

    GRIDGENERATOR_EXPORT real getStartX(uint level) const;
    GRIDGENERATOR_EXPORT real getStartY(uint level) const;
    GRIDGENERATOR_EXPORT real getStartZ(uint level) const;

    GRIDGENERATOR_EXPORT real getEndX(uint level) const;
    GRIDGENERATOR_EXPORT real getEndY(uint level) const;
    GRIDGENERATOR_EXPORT real getEndZ(uint level) const;

    GRIDGENERATOR_EXPORT std::vector<SPtr<Grid> > getGrids() const;
    GRIDGENERATOR_EXPORT void buildGrids(LbmOrGks lbmOrGks, bool enableThinWalls = false);

    GRIDGENERATOR_EXPORT void setNumberOfLayers( uint numberOfLayersFine, uint numberOfLayersBetweenLevels );


private:
    void addGridToList(SPtr<Grid> grid);
    bool coarseGridExists() const;

    SPtr<Grid> makeGrid(Object* gridShape, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, uint level) const;


    SPtr<GridFactory> gridFactory;
    Object* solidObject;

    uint numberOfLayersFine;
    uint numberOfLayersBetweenLevels;

    SPtr<BoundingBox> subDomainBox;

};

#endif

