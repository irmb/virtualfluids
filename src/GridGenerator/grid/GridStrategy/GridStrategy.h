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
//! \file GridStrategy.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef GRID_STRATEGY_H
#define GRID_STRATEGY_H

#include "Core/LbmOrGks.h"

#include "global.h"

#include "grid/Field.h"

struct Vertex;
class GridImp;

class VF_PUBLIC GridStrategy
{
public:
    virtual ~GridStrategy() {}

    virtual void allocateFieldMemory(Field* field) = 0;
    virtual void freeFieldMemory(Field* field) = 0;

    virtual void allocateGridMemory(SPtr<GridImp> grid) = 0;

    virtual void initalNodesToOutOfGrid(SPtr<GridImp> grid) = 0;
    virtual void findInnerNodes(SPtr<GridImp> grid) = 0;

	virtual void findEndOfGridStopperNodes(SPtr<GridImp> grid) = 0;

    virtual void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) = 0;


    virtual void freeMemory(SPtr<GridImp> grid) = 0;

};

#endif
