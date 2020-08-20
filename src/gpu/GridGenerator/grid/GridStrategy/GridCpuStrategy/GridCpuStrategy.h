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
//! \file GridCpuStrategy.h
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#ifndef GRID_CPU_STRATEGY_H
#define GRID_CPU_STRATEGY_H

#include "global.h"

#include "grid/GridStrategy/GridStrategy.h"

class GridImp;
class TriangularMesh;

class GRIDGENERATOR_EXPORT GridCpuStrategy : public GridStrategy
{
public:
    virtual ~GridCpuStrategy() {};

    void allocateGridMemory(SPtr<GridImp> grid) override;

    void initalNodesToOutOfGrid(SPtr<GridImp> grid) override;
    void findInnerNodes(SPtr<GridImp> grid) override;
	void findEndOfGridStopperNodes(SPtr<GridImp> grid) override;

    void freeMemory(SPtr<GridImp> grid) override;



    virtual void copyDataFromGPU() {};

protected:
    static void findForNeighborsNewIndices(SPtr<GridImp> grid);
public:
    void allocateFieldMemory(Field* field) override;
    void freeFieldMemory(Field* field) override;
    void findSparseIndices(SPtr<GridImp> coarseGrid, SPtr<GridImp> fineGrid) override;

};

#endif