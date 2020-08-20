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
//! \file BoundaryCondition.cpp
//! \ingroup BoundaryCondition
//! \author Stephan Lenz
//=======================================================================================
#include "BoundaryCondition.h"

#include <memory>
#include <vector>

#include "GridGenerator/grid/NodeValues.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseAllocator.h"
#include "DataBase/DataBaseStruct.h"

BoundaryCondition::BoundaryCondition( SPtr<DataBase> dataBase )
    : myAllocator ( dataBase->myAllocator )
{
      numberOfCells = INVALID_INDEX;
      ghostCells    = nullptr;
      domainCells   = nullptr;
      secondCells   = nullptr;
}

BoundaryCondition::~BoundaryCondition()
{
    this->myAllocator->freeMemory( *this );
}

void BoundaryCondition::findBoundaryCells(GksMeshAdapter & adapter, bool allowGhostCells, std::function<bool(Vec3)> boundaryFinder)
{
    this->myAllocator->freeMemory( *this );

    std::vector<uint> ghostCells;
    std::vector<uint> domainCells;
    std::vector<uint> secondCells;

    numberOfCellsPerLevel.resize( adapter.numberOfLevels );
    startOfCellsPerLevel.resize ( adapter.numberOfLevels );

    for( auto& n : numberOfCellsPerLevel ) n = 0;

    for( uint level = 0; level < adapter.numberOfLevels; level++ )
    {
        uint startIdx = adapter.startOfCellsPerLevel[level] 
                      + adapter.numberOfBulkCellsPerLevel[level];

        uint endIdx   = adapter.startOfCellsPerLevel[level] 
                      + adapter.numberOfCellsPerLevel[level];

        for( uint cellIdx = startIdx ; cellIdx < endIdx; cellIdx++ )
        {
            MeshCell& cell = adapter.cells[ cellIdx ];

            if( !boundaryFinder( cell.cellCenter ) ) continue;

            if( cell.type != STOPPER_OUT_OF_GRID && cell.type != STOPPER_OUT_OF_GRID_BOUNDARY && cell.type != STOPPER_SOLID ) continue;

            // look in all directions
            uint maximalSearchDirection = 27;

            // in case of Flux BC look only at face neighbors
            if( this->isFluxBC() ) maximalSearchDirection = 6;

            for( uint idx = 0; idx < maximalSearchDirection; idx++ )
            {
                uint neighborCellIdx = cell.cellToCell[ idx ];

                if( neighborCellIdx == INVALID_INDEX ) continue;

                MeshCell& neighborCell = adapter.cells[ neighborCellIdx ];

                bool neighborCellIsFluid = neighborCell.type != STOPPER_OUT_OF_GRID && 
                                           neighborCell.type != STOPPER_OUT_OF_GRID_BOUNDARY && 
                                           neighborCell.type != STOPPER_SOLID;

                bool neighborCellIsValidGhostCell = !this->isFluxBC() && allowGhostCells && !boundaryFinder( neighborCell.cellCenter );

                if( neighborCellIsFluid || neighborCellIsValidGhostCell )
                {
                    ghostCells.push_back ( cellIdx );
                    domainCells.push_back( neighborCellIdx );

                    this->numberOfCellsPerLevel[ level ]++;

                    if( this->secondCellsNeeded() )
                    {
                        secondCells.push_back( neighborCell.cellToCell[ idx ] );
                    }

                    cell.isWall      = this->isWall();
                    cell.isFluxBC    = this->isFluxBC();
                    cell.isInsulated = this->isInsulated();

                    break;
                }
            }
        }
    }

    startOfCellsPerLevel[ 0 ] = 0;

    for( uint level = 1; level < adapter.numberOfLevels; level++ )
    {
        startOfCellsPerLevel[ level ] = startOfCellsPerLevel [ level - 1 ]
                                      + numberOfCellsPerLevel[ level - 1 ];
    }

    this->numberOfCells = ghostCells.size();

    this->myAllocator->allocateMemory( shared_from_this(), ghostCells, domainCells, secondCells );
}

bool BoundaryCondition::isFluxBC()
{
    return false;
}

bool BoundaryCondition::isInsulated()
{
    return false;
}

bool BoundaryCondition::secondCellsNeeded()
{
    return false;
}
