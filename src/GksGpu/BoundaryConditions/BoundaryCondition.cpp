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

                    cell.isWall   = this->isWall();
                    cell.isFluxBC = this->isFluxBC();

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

bool BoundaryCondition::secondCellsNeeded()
{
    return false;
}
