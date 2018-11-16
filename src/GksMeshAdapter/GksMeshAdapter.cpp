#include "GksMeshAdapter.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <fstream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <iomanip>

#include "core/Logger/Logger.h"

#include "GridGenerator/grid/distributions/D3Q27.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/NodeValues.h"
#include "GridGenerator/utilities/math/Math.h"

#include "MeshCell.h"
#include "MeshFace.h"

void GksMeshAdapter::inputGrid(SPtr<MultipleGridBuilder> gridBuilder)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "inputGrid()" << "\n";

    this->numberOfLevels = gridBuilder->getNumberOfGridLevels();

    std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    this->dxCoarse = grids[0]->getDelta();

    //////////////////////////////////////////////////////////////////////////

    this->gridToMesh.resize( gridBuilder->getNumberOfGridLevels() );

    for( uint level = 0; level < gridBuilder->getNumberOfGridLevels(); level++ ){
        this->gridToMesh[level].resize( grids[level]->getSize() );

        for( auto& cellIdx : this->gridToMesh[level] ) cellIdx = INVALID_INDEX;
    }

    //////////////////////////////////////////////////////////////////////////
    //
    //    I d e n t i f y    C e l l s    i n    L B - G r i d
    //
    //////////////////////////////////////////////////////////////////////////

    uint numberOfCells = 0;

    for( uint level = 0; level < gridBuilder->getNumberOfGridLevels(); level++ ){
        for( uint gridIdx = 0; gridIdx < grids[level]->getSize(); gridIdx++ ){
            if (grids[level]->getFieldEntry(gridIdx)  != STOPPER_COARSE_UNDER_FINE &&
                grids[level]->getFieldEntry(gridIdx)  != STOPPER_SOLID &&
                grids[level]->getFieldEntry(gridIdx)  != INVALID_COARSE_UNDER_FINE &&
                grids[level]->getFieldEntry(gridIdx)  != INVALID_OUT_OF_GRID &&
                grids[level]->getFieldEntry(gridIdx)  != INVALID_SOLID )
            {
                this->gridToMesh[level][gridIdx] = numberOfCells++;
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////
    //
    //    S e t    M e s h    t o    G r i d    i n f o r m a t i o n
    //
    //////////////////////////////////////////////////////////////////////////

    this->cells.resize( numberOfCells );

    for( uint level = 0; level < gridBuilder->getNumberOfGridLevels(); level++ ){
        for( uint gridIdx = 0; gridIdx < grids[level]->getSize(); gridIdx++ ){
            if ( this->gridToMesh[level][gridIdx] != INVALID_INDEX ){

                uint cellIdx = gridToMesh[level][gridIdx];

                MeshCell& cell = this->cells[ cellIdx ];

                cell.level   = level;
                cell.gridIdx = gridIdx;

                cell.type = grids[level]->getFieldEntry(gridIdx);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////
}

void GksMeshAdapter::findQuadtreeConnectivity(SPtr<MultipleGridBuilder> gridBuilder)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "findQuadtreeConnectivity()" << "\n";

    std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        if( cell.type == FLUID_FCC || cell.type == FLUID_CFC ){

            real x, y, z;
            grids[cell.level]->transIndexToCoords(cell.gridIdx, x, y, z);

            real d = 0.25 * grids[cell.level]->getDelta();

            for( uint idx = 0; idx < 8; idx++ )
            {
                Distribution dirs = DistributionHelper::getDistribution27();

                real xSign = dirs.directions[idx + 19][0];
                real ySign = dirs.directions[idx + 19][1];
                real zSign = dirs.directions[idx + 19][2];

                cell.children[ idx ] = this->gridToMesh[cell.level+1][ grids[cell.level+1]->transCoordToIndex( x + xSign * d, 
                                                                                                               y + ySign * d, 
                                                                                                               z + zSign * d ) ];
            }

            // register parent
            for( uint child = 0; child < 8; child++ )
                this->cells[ cell.children[child] ].parent = cellIdx;

            // set correct type for CFF cells
            for( uint child = 0; child < 8; child++ )
                if( this->cells[ cell.children[child] ].type != FLUID_FCF ) 
                    this->cells[ cell.children[child] ].type = FLUID_CFF;

        }
    }
}

void GksMeshAdapter::findCellToCellConnectivity(SPtr<MultipleGridBuilder> gridBuilder)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "findCellToCellConnectivity()" << "\n";

    std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        real x, y, z;
        grids[cell.level]->transIndexToCoords(cell.gridIdx, x, y, z);

        real d = grids[cell.level]->getDelta();

        for( uint idx = 0; idx < 27; idx++ )
        {
            if( idx == DIR_27_ZERO ) continue;

            Distribution dirs = DistributionHelper::getDistribution27();

            int xSign = dirs.directions[idx][0];
            int ySign = dirs.directions[idx][1];
            int zSign = dirs.directions[idx][2];

            uint neighborGridIdx = grids[cell.level]->transCoordToIndex( x + xSign * d, 
                                                                         y + ySign * d, 
                                                                         z + zSign * d );

            if( neighborGridIdx == INVALID_INDEX || this->gridToMesh[cell.level][neighborGridIdx] == INVALID_INDEX ){
                if( !cell.isCoarseGhostCell() && cell.type != BC_SOLID )
                    cell.isGhostCell = true;

                continue;
            }

            cell.cellToCell[ idx ] = this->gridToMesh[cell.level][neighborGridIdx];
        }
    }
}

void GksMeshAdapter::countCells()
{
    this->numberOfCellsPerLevel    .resize( this->numberOfLevels );
    this->numberOfBulkCellsPerLevel.resize( this->numberOfLevels );
    this->startOfCellsPerLevel     .resize( this->numberOfLevels );

    for( auto& i : this->numberOfCellsPerLevel     ) i = 0;
    for( auto& i : this->numberOfBulkCellsPerLevel ) i = 0;
    for( auto& i : this->startOfCellsPerLevel      ) i = 0;

    uint level = 0;
    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
        MeshCell& cell = this->cells[ cellIdx ];

        if( cell.level != level ) level++;

        this->numberOfCellsPerLevel[ level ]++; 

        if( ! ( cell.isGhostCell || cell.isCoarseGhostCell() ) )
            this->numberOfBulkCellsPerLevel[ level ]++;
    }

    for( uint level = 1; level < this->numberOfLevels; level++ )
        this->startOfCellsPerLevel[ level ] = this->startOfCellsPerLevel[ level-1 ] + this->numberOfCellsPerLevel[ level-1 ];
}

void GksMeshAdapter::partitionCells()
{
    for( uint level = 0; level < this->numberOfLevels; level++ ){

        std::vector<uint> idxMap( this->cells.size() );
        std::iota( idxMap.begin(), idxMap.end(), 0 );

        // partition idxMap
        std::stable_partition(  idxMap.begin() + this->startOfCellsPerLevel[level], 
                                idxMap.begin() + this->startOfCellsPerLevel[level] 
                                               + this->numberOfCellsPerLevel[level], 
                                [this](int lhs){ 
                                    return ! ( this->cells[ lhs ].isGhostCell || this->cells[ lhs ].isCoarseGhostCell() );
                                }
                             );

        // invert idxMap
        {
            std::vector<uint> buffer = idxMap;
            for( uint idx = 0; idx < idxMap.size(); idx ++ )
                idxMap[ buffer[idx] ] = idx;
        }

        // partition cell list
        std::stable_partition(  this->cells.begin() + this->startOfCellsPerLevel[level], 
                                this->cells.begin() + this->startOfCellsPerLevel[level] 
                                                    + this->numberOfCellsPerLevel[level], 
                                [this](MeshCell lhs){ 
                                    return ! ( lhs.isGhostCell || lhs.isCoarseGhostCell() );
                                }
                             );

        this->refreshCellConnectivity( idxMap );
    }
}

void GksMeshAdapter::refreshCellConnectivity(const std::vector<uint>& idxMap)
{
    for( auto& cell : this->cells ){
        for( uint idx = 0; idx < 27; idx++ )
            if( cell.cellToCell[ idx ] != INVALID_INDEX )
                cell.cellToCell[ idx ] = idxMap[ cell.cellToCell[ idx ] ];

        if( cell.parent != INVALID_INDEX )
            cell.parent = idxMap[ cell.parent ];

        for( uint idx = 0; idx < 8; idx++ )
            if( cell.children[ idx ] != INVALID_INDEX )
                cell.children[ idx ] = idxMap[ cell.children[ idx ] ];
    }

    for( auto& grid : this->gridToMesh ){
        for( auto& cellIdx : grid ){
            if( cellIdx != INVALID_INDEX )
                cellIdx = idxMap[ cellIdx ];
        }
    }
}

void GksMeshAdapter::findCornerCells(SPtr<MultipleGridBuilder> gridBuilder, real z0)
{
    //SPtr<Grid> grid = gridBuilder->getGrids()[0];
    //
    //this->cornerCells[0] = this->gridToMesh[ 0 ][ grid->transCoordToIndex( grid->getStartX(), grid->getStartY(), z0 ) ];
    //this->cornerCells[1] = this->gridToMesh[ 0 ][ grid->transCoordToIndex( grid->getEndX()  , grid->getStartY(), z0 ) ];
    //this->cornerCells[2] = this->gridToMesh[ 0 ][ grid->transCoordToIndex( grid->getEndX()  , grid->getEndY()  , z0 ) ];
    //this->cornerCells[3] = this->gridToMesh[ 0 ][ grid->transCoordToIndex( grid->getStartX(), grid->getEndY()  , z0 ) ];
}

void GksMeshAdapter::generateNodes(SPtr<MultipleGridBuilder> gridBuilder)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "generateNodes()" << "\n";

    std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    nodes.reserve( 2 * this->cells.size() );

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        if( cell.type == STOPPER_SOLID ) continue;

        real x, y, z;
        grids[cell.level]->transIndexToCoords(cell.gridIdx, x, y, z);

        real d = 0.5 * grids[cell.level]->getDelta();

        std::array<Vec3,8> dir;

        for( uint idx = 0; idx < 8; idx++ )
        {
            if( cell.cellToNode[idx] == INVALID_INDEX )
            {
                Distribution dirs = DistributionHelper::getDistribution27();

                real dx = dirs.directions[idx + 19][0] * d;
                real dy = dirs.directions[idx + 19][1] * d;
                real dz = dirs.directions[idx + 19][2] * d;

                nodes.push_back( Vec3( x + dx, y + dy, z + dz ) );

                cell.cellToNode[idx] = nodes.size()-1;

                //// register new node at neighbor cells on same level
                for (uint idx = 0; idx < 8; idx++)
                {
                    Distribution dirs = DistributionHelper::getDistribution27();

                    real dxNeighbor = -dirs.directions[idx + 19][0] * d;
                    real dyNeighbor = -dirs.directions[idx + 19][1] * d;
                    real dzNeighbor = -dirs.directions[idx + 19][2] * d;

                    real xNeighbor = nodes.back().x + dxNeighbor;
                    real yNeighbor = nodes.back().y + dyNeighbor;
                    real zNeighbor = nodes.back().z + dzNeighbor;

                    uint neighborGridIdx = grids[cell.level]->transCoordToIndex(xNeighbor, yNeighbor, zNeighbor);

                    if ( neighborGridIdx == INVALID_INDEX ) continue;

                    uint neighborIdx = gridToMesh[cell.level][neighborGridIdx];

                    if ( neighborIdx != INVALID_INDEX )
                    {
                        this->cells[ neighborIdx ].cellToNode[idx] = nodes.size() - 1;
                    }
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////

        //if( cell.type == FLUID_FCC || cell.type == FLUID_CFC ){

        //    // register edge nodes at child cells
        //    this->cells[ cell.children[0] ].cellToNode[0] = cell.cellToNode[0];
        //    this->cells[ cell.children[1] ].cellToNode[1] = cell.cellToNode[1];
        //    this->cells[ cell.children[2] ].cellToNode[2] = cell.cellToNode[2];
        //    this->cells[ cell.children[3] ].cellToNode[3] = cell.cellToNode[3];                    

        //    std::array<arraytypes::Vec2,4> dir;

        //    dir[0] = arraytypes::Vec2( - d, 0.0 );
        //    dir[1] = arraytypes::Vec2( 0.0, - d );
        //    dir[2] = arraytypes::Vec2(   d, 0.0 );
        //    dir[3] = arraytypes::Vec2( 0.0,   d );

        //    // introduce new face nodes
        //    if( cell.cellToFaceNode[0] == -1 ){
        //        nodes.push_back( arraytypes::Vec2( x + dir[0].x, y + dir[0].y ) );

        //        cell.cellToFaceNode[0] = nodes.size()-1;
        //            
        //        // register face node at adjacent cell
        //        if( cell.cellToCell[0] != -1 ) this->cells[ cell.cellToCell[0] ].cellToFaceNode[2] = nodes.size()-1;

        //    }
        //            
        //    if( cell.cellToFaceNode[1] == -1 ){
        //        nodes.push_back( arraytypes::Vec2( x + dir[1].x, y + dir[1].y ) );

        //        cell.cellToFaceNode[1] = nodes.size()-1;
        //            
        //        // register face node at adjacent cell
        //        if( cell.cellToCell[1] != -1 ) this->cells[ cell.cellToCell[1] ].cellToFaceNode[3] = nodes.size()-1;

        //    }
        //            
        //    if( cell.cellToFaceNode[2] == -1 ){
        //        nodes.push_back( arraytypes::Vec2( x + dir[2].x, y + dir[2].y ) );

        //        cell.cellToFaceNode[2] = nodes.size()-1;
        //            
        //        // register face node at adjacent cell
        //        if( cell.cellToCell[2] != -1 ) this->cells[ cell.cellToCell[2] ].cellToFaceNode[0] = nodes.size()-1;

        //    }
        //            
        //    if( cell.cellToFaceNode[3] == -1 ){
        //        nodes.push_back( arraytypes::Vec2( x + dir[3].x, y + dir[3].y ) );

        //        cell.cellToFaceNode[3] = nodes.size()-1;
        //            
        //        // register face node at adjacent cell
        //        if( cell.cellToCell[3] != -1 ) this->cells[ cell.cellToCell[3] ].cellToFaceNode[1] = nodes.size()-1;

        //    }

        //    // register face node at child cells
        //    this->cells[ cell.children[3] ].cellToNode[0] = cell.cellToFaceNode[0];
        //    this->cells[ cell.children[0] ].cellToNode[3] = cell.cellToFaceNode[0];

        //    this->cells[ cell.children[0] ].cellToNode[1] = cell.cellToFaceNode[1];
        //    this->cells[ cell.children[1] ].cellToNode[0] = cell.cellToFaceNode[1];

        //    this->cells[ cell.children[1] ].cellToNode[2] = cell.cellToFaceNode[2];
        //    this->cells[ cell.children[2] ].cellToNode[1] = cell.cellToFaceNode[2];

        //    this->cells[ cell.children[2] ].cellToNode[3] = cell.cellToFaceNode[3];
        //    this->cells[ cell.children[3] ].cellToNode[2] = cell.cellToFaceNode[3];
        //}
    }
}

void GksMeshAdapter::inputQs(SPtr<MultipleGridBuilder> gridBuilder)
{

    //std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    //for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    MeshCell& cell = this->cells[ cellIdx ];

    //    if( cell.type != BC_SOLID ) continue;

    //    cell.Qs[0] = grids[cell.level]->getQValue(cell.gridIdx, 7);
    //    cell.Qs[1] = grids[cell.level]->getQValue(cell.gridIdx, 8);
    //    cell.Qs[2] = grids[cell.level]->getQValue(cell.gridIdx, 6);
    //    cell.Qs[3] = grids[cell.level]->getQValue(cell.gridIdx, 9);
    //}
}

void GksMeshAdapter::morphNodes()
{
    //this->displacement.resize(   this->nodes.size() );
    //this->nDisplacements.resize( this->nodes.size() );

    //for( uint n : this->nDisplacements ) n = 0;

    //for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    MeshCell& cell = this->cells[ cellIdx ];

    //    if( cell.type != BC_SOLID ) continue;

    //    std::array<arraytypes::Vec2,4> dir;
    //    dir[0] = arraytypes::Vec2( - getDx(cell.level), - getDx(cell.level) );
    //    dir[1] = arraytypes::Vec2(   getDx(cell.level), - getDx(cell.level) );
    //    dir[2] = arraytypes::Vec2(   getDx(cell.level),   getDx(cell.level) );
    //    dir[3] = arraytypes::Vec2( - getDx(cell.level),   getDx(cell.level) );

    //    for( uint idx = 0; idx < 4; idx ++ ){
    //    
    //        //if( vf::Math::equal( cell.Qs[idx], 0.0 ) ) continue;
    //        if( cell.Qs[idx] < 0.0 || cell.Qs[idx] > 1.0 ) continue;

    //        this->displacement[ cell.cellToNode[idx] ].x += ( cell.Qs[idx] - 0.5 ) * dir[idx].x;
    //        this->displacement[ cell.cellToNode[idx] ].y += ( cell.Qs[idx] - 0.5 ) * dir[idx].y;
    //    
    //        this->nDisplacements[ cell.cellToNode[idx] ]++;
    //    }
    //}

    //for( uint nodeIdx = 0; nodeIdx < this->nodes.size(); nodeIdx++ ){
    //
    //    if( this->nDisplacements[nodeIdx] == 0 ) continue;
    //
    //    nodes[ nodeIdx ].x += this->displacement[ nodeIdx ].x / double( this->nDisplacements[ nodeIdx ] );
    //    nodes[ nodeIdx ].y += this->displacement[ nodeIdx ].y / double( this->nDisplacements[ nodeIdx ] );
    //}
}

void GksMeshAdapter::repairTriangularCells()
{
    //for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    MeshCell& cell = this->cells[ cellIdx ];

    //    if( cell.type != BC_SOLID ) continue;

    //    arraytypes::uint4 cellToNodeTMP = cell.cellToNode;

    //    for( uint idx = 0; idx < 4; idx++ ){
    //    
    //        arraytypes::Vec2 v1, v2;

    //        v1.x =   this->nodes[ cell.cellToNode[ idx ] ].x - this->nodes[ cell.cellToNode[ ( idx - 1 ) % 4 ] ].x;
    //        v1.y =   this->nodes[ cell.cellToNode[ idx ] ].y - this->nodes[ cell.cellToNode[ ( idx - 1 ) % 4 ] ].y;
    //        v2.x = - this->nodes[ cell.cellToNode[ idx ] ].x + this->nodes[ cell.cellToNode[ ( idx + 1 ) % 4 ] ].x;
    //        v2.y = - this->nodes[ cell.cellToNode[ idx ] ].y + this->nodes[ cell.cellToNode[ ( idx + 1 ) % 4 ] ].y;

    //        double v1Length = sqrt( v1.x*v1.x + v1.y*v1.y );
    //        double v2Length = sqrt( v2.x*v2.x + v2.y*v2.y );

    //        double v1DotV2  = v1.x*v2.x + v1.y*v2.y;

    //        bool isConcave = v1.x*v2.y - v1.y*v2.x < 0.0;

    //        double argument = v1DotV2 / ( v1Length * v2Length );

    //        bool isCollinear = false;
    //        
    //        // the acos function does not like the argument 1.0, hence the special treatment
    //        if( fabs( v1DotV2 / ( v1Length * v2Length ) - 1.0 ) < 1.0e-10 )
    //            isCollinear = true;
    //        else
    //            isCollinear = acos( v1DotV2 / ( v1Length * v2Length ) ) < M_PI / 4.0;

    //        //////////////////////////////////////////////////////////////////////////

    //        uint neighborIdx1 = cell.cellToCell[   idx           ];
    //        uint neighborIdx2 = cell.cellToCell[ ( idx + 1 ) % 4 ];

    //        bool neighborCellsExist = ( ( neighborIdx1 != INVALID_IDX && this->cells[ neighborIdx1 ].type != STOPPER_SOLID ) ||
    //                                    ( neighborIdx2 != INVALID_IDX && this->cells[ neighborIdx2 ].type != STOPPER_SOLID ) );

    //        //////////////////////////////////////////////////////////////////////////

    //        if( isConcave || isCollinear ){
    //            if( !neighborCellsExist ){
    //                cellToNodeTMP[idx] = INVALID_IDX;
    //            }
    //            else{
    //                nodes[ cell.cellToNode[ idx ] ].x -= this->displacement[ cell.cellToNode[ idx ] ].x / double( this->nDisplacements[ cell.cellToNode[ idx ] ] );
    //                nodes[ cell.cellToNode[ idx ] ].y -= this->displacement[ cell.cellToNode[ idx ] ].y / double( this->nDisplacements[ cell.cellToNode[ idx ] ] );
    //            }
    //        }
    //    }

    //    cell.cellToNode = cellToNodeTMP;
    //}
}

void GksMeshAdapter::computeCellGeometry()
{    
    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
        
        MeshCell& cell = this->cells[ cellIdx ];

        Vec3 cellCenter;

        for( uint node = 0; node < 8; node++ ){
            cellCenter = cellCenter + this->nodes[ cell.cellToNode[node] ];
        }

        cell.cellCenter.x = cellCenter.x / eight;
        cell.cellCenter.y = cellCenter.y / eight;
        cell.cellCenter.z = cellCenter.z / eight;

    //    MeshCell& cell = this->cells[ cellIdx ];

    //    uint nNodes = 0;
    //    for( uint idx : cell.cellToNode ) if( idx != INVALID_IDX ) nNodes++;

    //    if( nNodes == 3 ){
    //    
    //        arraytypes::uint3 cellToNode;
    //        for( uint idx = 0, counter = 0; idx < 4; idx++ )
    //            if( cell.cellToNode[idx] != INVALID_IDX ) 
    //                cellToNode[counter++] = cell.cellToNode[idx];
    //    
    //        cell.cellCenter.x = ( this->nodes[ cellToNode[0] ].x
    //                            + this->nodes[ cellToNode[1] ].x
    //                            + this->nodes[ cellToNode[2] ].x ) / 3.0;
    //        cell.cellCenter.y = ( this->nodes[ cellToNode[0] ].y
    //                            + this->nodes[ cellToNode[1] ].y
    //                            + this->nodes[ cellToNode[2] ].y ) / 3.0;

    //        cell.cellVolume = 0.5 * ( this->nodes[ cellToNode[0] ].x * ( this->nodes[ cellToNode[1] ].y - this->nodes[ cellToNode[2] ].y )
    //                                + this->nodes[ cellToNode[1] ].x * ( this->nodes[ cellToNode[2] ].y - this->nodes[ cellToNode[0] ].y )
    //                                + this->nodes[ cellToNode[2] ].x * ( this->nodes[ cellToNode[0] ].y - this->nodes[ cellToNode[1] ].y ) );

    //    }
    //    else if( nNodes == 4 ){
    //        
    //        arraytypes::Vec2 triCenter[2];
    //        triCenter[0].x =  (this->nodes[ cell.cellToNode[0] ].x + this->nodes[ cell.cellToNode[1] ].x +                                       this->nodes[ cell.cellToNode[3] ].x) / 3.0;
    //        triCenter[0].y =  (this->nodes[ cell.cellToNode[0] ].y + this->nodes[ cell.cellToNode[1] ].y +                                       this->nodes[ cell.cellToNode[3] ].y) / 3.0;
    //        triCenter[1].y =  (                                      this->nodes[ cell.cellToNode[1] ].y + this->nodes[ cell.cellToNode[2] ].y + this->nodes[ cell.cellToNode[3] ].y) / 3.0;
    //        triCenter[1].x =  (                                      this->nodes[ cell.cellToNode[1] ].x + this->nodes[ cell.cellToNode[2] ].x + this->nodes[ cell.cellToNode[3] ].x) / 3.0;

    //        double triVolume[2];
    //        triVolume[0] = 0.5 * fabs( this->nodes[ cell.cellToNode[0] ].x * ( this->nodes[ cell.cellToNode[1] ].y - this->nodes[ cell.cellToNode[3] ].y )
    //                                 + this->nodes[ cell.cellToNode[1] ].x * ( this->nodes[ cell.cellToNode[3] ].y - this->nodes[ cell.cellToNode[0] ].y )
    //                                 + this->nodes[ cell.cellToNode[3] ].x * ( this->nodes[ cell.cellToNode[0] ].y - this->nodes[ cell.cellToNode[1] ].y ) );
    //        triVolume[1] = 0.5 * fabs( this->nodes[ cell.cellToNode[2] ].x * ( this->nodes[ cell.cellToNode[3] ].y - this->nodes[ cell.cellToNode[1] ].y )
    //                                 + this->nodes[ cell.cellToNode[3] ].x * ( this->nodes[ cell.cellToNode[1] ].y - this->nodes[ cell.cellToNode[2] ].y )
    //                                 + this->nodes[ cell.cellToNode[1] ].x * ( this->nodes[ cell.cellToNode[2] ].y - this->nodes[ cell.cellToNode[3] ].y ) );

    //        cell.cellVolume   = triVolume[0] + triVolume[1];
    //        cell.cellCenter.x = ( triCenter[0].x * triVolume[0] + triCenter[1].x * triVolume[1] ) / cell.cellVolume;
    //        cell.cellCenter.y = ( triCenter[0].y * triVolume[0] + triCenter[1].y * triVolume[1] ) / cell.cellVolume;
    //    }
    }
}

void GksMeshAdapter::generateSolidGhostCells()
{
    //std::vector<MeshCell> newCells;

    //for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    MeshCell& cell = this->cells[ cellIdx ];

    //    if( cell.type != BC_SOLID ) continue;
    //    
    //    for( uint idx = 0; idx < 4; idx++ ){
    //        
    //        if( cell.cellToCell[ idx ] == INVALID_IDX ){
    //            
    //            cell.cellToCell[ idx ] = this->cells.size() + newCells.size();

    //            if( cell.cellToNode[ idx ] == INVALID_IDX )
    //                cell.cellToCell[ ( idx + 1 ) % 4 ] = this->cells.size() + newCells.size();

    //            if( cell.cellToNode[ ( idx - 1 ) % 4 ] == INVALID_IDX )
    //                cell.cellToCell[ ( idx - 1 ) % 4 ] = this->cells.size() + newCells.size();

    //            uint neighborIdx = cell.cellToCell[ ( idx + 1 ) % 4 ];
    //            if( neighborIdx != INVALID_IDX && neighborIdx < this->cells.size() )
    //                this->cells[ neighborIdx ].cellToCell[ 4 + ( idx - 1 ) % 4 ] = this->cells.size() + newCells.size();

    //            neighborIdx = cell.cellToCell[ ( idx - 1 ) % 4 ];
    //            if( neighborIdx != INVALID_IDX && neighborIdx < this->cells.size() )
    //                this->cells[ neighborIdx ].cellToCell[ 4 + ( idx + 0 ) % 4 ] = this->cells.size() + newCells.size();

    //            newCells.push_back( MeshCell() );
    //            
    //            newCells.back().isGhostCell = true;
    //            newCells.back().type        = STOPPER_SOLID;
    //            newCells.back().level       = cell.level;
    //        }
    //    }
    //}

    //this->cells.insert( this->cells.end(), newCells.begin(), newCells.end() );
}

void GksMeshAdapter::computeGhostCellGeometry()
{
    //for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    MeshCell& cell = this->cells[ cellIdx ];

    //    if( cell.type != BC_SOLID ) continue;

    //    for( uint idx = 0; idx < 4; idx++ ){
    //        
    //        MeshCell& neighborCell = this->cells[ cell.cellToCell[ idx ] ];

    //        if( neighborCell.type != STOPPER_SOLID ) continue;
    //    
    //        arraytypes::Vec2 faceCenter;

    //        uint node0 = cell.cellToNode[ ( idx - 1 ) % 4 ];
    //        uint node1 = cell.cellToNode[   idx           ];

    //        if( node0 == INVALID_IDX ) continue;

    //        if( node1 == INVALID_IDX ) node1 = cell.cellToNode[ ( idx + 1 ) % 4 ];

    //        faceCenter.x = 0.5 * ( this->nodes[ node0 ].x + this->nodes[ node1 ].x );
    //        faceCenter.y = 0.5 * ( this->nodes[ node0 ].y + this->nodes[ node1 ].y );

    //        arraytypes::Vec2 newNode;

    //        //newNode.x = 2.0 * faceCenter.x - cell.cellCenter.x;
    //        //newNode.y = 2.0 * faceCenter.y - cell.cellCenter.y;

    //        //////////////////////////////////////////////////////////////////////////

    //        arraytypes::Vec2 v01 ( this->nodes[ node1 ].x - this->nodes[ node0 ].x,
    //                               this->nodes[ node1 ].y - this->nodes[ node0 ].y );

    //        arraytypes::Vec2 vD0 ( this->nodes[ node0 ].x - cell.cellCenter.x,
    //                               this->nodes[ node0 ].y - cell.cellCenter.y );

    //        double alpha = - ( vD0.x * v01.x + vD0.y * v01.y ) / ( v01.x * v01.x + v01.y * v01.y );

    //        arraytypes::Vec2 vDP( this->nodes[ node0 ].x + alpha * v01.x - cell.cellCenter.x,
    //                              this->nodes[ node0 ].y + alpha * v01.y - cell.cellCenter.y );

    //        newNode.x = cell.cellCenter.x + 2.0 * vDP.x;
    //        newNode.y = cell.cellCenter.y + 2.0 * vDP.y;

    //        //////////////////////////////////////////////////////////////////////////

    //        neighborCell.cellCenter = newNode;
    //        
    //        neighborCell.cellToNode[0] = node1;
    //        neighborCell.cellToNode[1] = node0;
    //        neighborCell.cellToNode[2] = this->nodes.size();

    //        // dummy value to identify it as ghost cell
    //        neighborCell.cellVolume = -1.0;

    //        this->nodes.push_back( newNode );
    //    }
    //}
}

void GksMeshAdapter::generateFaces(SPtr<MultipleGridBuilder> gridBuilder)
{
    std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    this->faces.reserve( 2 * this->cells.size() );

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        if( cell.type == BC_SOLID || cell.type == STOPPER_SOLID ) continue;

        // generate faces in positive direction
        for( uint neighborIdx = 0; neighborIdx < 6; neighborIdx += 2 ){

            if( cell.faceExists[ neighborIdx ] ) continue;

            if( cell.cellToCell[ neighborIdx ] == INVALID_INDEX ) continue;

            uint neighborCellIdx = cell.cellToCell[ neighborIdx ];

            MeshCell& neighborCell = this->cells[ neighborCellIdx ];

            if( cell.isGhostCell && neighborCell.isGhostCell ) continue;

            if( cell.isCoarseGhostCell() || neighborCell.isCoarseGhostCell() ) continue;

            //////////////////////////////////////////////////////////////////////////

            MeshFace newFace;

            newFace.level = cell.level;

            if( neighborIdx == 0 )
            {
                newFace.faceToNode[ 0 ] = cell.cellToNode[ 3 ];
                newFace.faceToNode[ 1 ] = cell.cellToNode[ 1 ];
                newFace.faceToNode[ 2 ] = cell.cellToNode[ 0 ];
                newFace.faceToNode[ 3 ] = cell.cellToNode[ 2 ];
                newFace.orientation = 'x';
            }
            if( neighborIdx == 2 )
            {
                newFace.faceToNode[ 0 ] = cell.cellToNode[ 5 ];
                newFace.faceToNode[ 1 ] = cell.cellToNode[ 4 ];
                newFace.faceToNode[ 2 ] = cell.cellToNode[ 0 ];
                newFace.faceToNode[ 3 ] = cell.cellToNode[ 1 ];
                newFace.orientation = 'y';
            }
            if( neighborIdx == 4 )
            {
                newFace.faceToNode[ 0 ] = cell.cellToNode[ 6 ];
                newFace.faceToNode[ 1 ] = cell.cellToNode[ 2 ];
                newFace.faceToNode[ 2 ] = cell.cellToNode[ 0 ];
                newFace.faceToNode[ 3 ] = cell.cellToNode[ 4 ];
                newFace.orientation = 'z';
            }

            //////////////////////////////////////////////////////////////////////////

            cell.faceExists[ neighborIdx ] = true;

            // register face at neighbor
            for( uint idx = 0; idx < 6; idx++ ){
                if( neighborCell.cellToCell[ idx ] == cellIdx ){
                    neighborCell.faceExists[ idx ] = true;
                    break;
                }
            }

            //////////////////////////////////////////////////////////////////////////

            newFace.negCell = cellIdx;
            newFace.posCell = neighborCellIdx;

            //////////////////////////////////////////////////////////////////////////

            if ( cell.type == FLUID_CFF && neighborCell.type == FLUID_FCF ) newFace.negCellCoarse = cell.parent;
            if ( cell.type == FLUID_FCF && neighborCell.type == FLUID_CFF ) newFace.posCellCoarse = neighborCell.parent;

            //////////////////////////////////////////////////////////////////////////
            
            Vec3 faceCenter;

            for( uint node = 0; node < 4; node++ ){
                faceCenter = faceCenter + this->nodes[ newFace.faceToNode[node] ];
            }

            newFace.faceCenter.x = faceCenter.x / four;
            newFace.faceCenter.y = faceCenter.y / four;
            newFace.faceCenter.z = faceCenter.z / four;

            this->faces.push_back( newFace );
        }
    }
}

void GksMeshAdapter::generateSolidFaces(SPtr<MultipleGridBuilder> gridBuilder)
{

    //std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    //this->faces.reserve( 2 * this->cells.size() );

    //for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    //
    //    MeshCell& cell = this->cells[ cellIdx ];

    //    if( ! ( cell.type == BC_SOLID ) ) continue;

    //    for( uint neighborIdx = 0; neighborIdx < 4; neighborIdx++ ){

    //        if( cell.faceExists[ neighborIdx ] ) continue;

    //        if( cell.cellToCell[ neighborIdx ] == INVALID_IDX ) continue;

    //        uint neighborCellIdx = cell.cellToCell[ neighborIdx ];

    //        MeshCell& neighborCell = this->cells[ neighborCellIdx ];

    //        if( cell.isGhostCell && neighborCell.isGhostCell ) continue;

    //        if( cell.isCoarseGhostCell() || neighborCell.isCoarseGhostCell() ) continue;

    //        //////////////////////////////////////////////////////////////////////////

    //        MeshFace newFace;

    //        newFace.level = cell.level;

    //        uint node0 = cell.cellToNode[ ( neighborIdx - 1 ) % 4 ];
    //        uint node1 = cell.cellToNode[   neighborIdx           ];

    //        if( node0 == INVALID_IDX ) continue;

    //        if( node1 == INVALID_IDX ) node1 = cell.cellToNode[ ( neighborIdx + 1 ) % 4 ];

    //        newFace.faceToNode[ 0 ] = node0;
    //        newFace.faceToNode[ 1 ] = node1;

    //        //////////////////////////////////////////////////////////////////////////

    //        newFace.orientation = 'a';

    //        //////////////////////////////////////////////////////////////////////////

    //        cell.faceExists[ neighborIdx ] = true;

    //        // register face at neighbor
    //        for( uint idx = 0; idx < 4; idx++ ){
    //            if( neighborCell.cellToCell[ idx ] == cellIdx ){
    //                neighborCell.faceExists[ idx ] = true;
    //                break;
    //            }
    //        }

    //        //////////////////////////////////////////////////////////////////////////

    //        newFace.faceToCell[ 0 ] = cellIdx;
    //        newFace.faceToCell[ 1 ] = neighborCellIdx;
    //        
    //        //////////////////////////////////////////////////////////////////////////
    //        //////////////////////////////////////////////////////////////////////////

    //        if( false )
    //        {
    //            if (neighborCell.type != STOPPER_SOLID) {
    //                newFace.faceToCell[2] = cell.cellToCell[(neighborIdx + 1) % 4];
    //                newFace.faceToCell[3] = cell.cellToCell[(neighborIdx - 1) % 4];

    //                newFace.faceToCell[4] = cell.cellToCell[neighborIdx + 4];
    //                newFace.faceToCell[5] = cell.cellToCell[neighborIdx + 3];
    //            }
    //            else {
    //                newFace.faceToCell[2] = cell.cellToCell[(neighborIdx + 1) % 4];
    //                newFace.faceToCell[3] = cell.cellToCell[(neighborIdx - 1) % 4];

    //                if (newFace.faceToCell[2] == newFace.faceToCell[1])
    //                    newFace.faceToCell[2] = cell.cellToCell[(neighborIdx + 2) % 4];

    //                if (newFace.faceToCell[3] == newFace.faceToCell[1])
    //                    newFace.faceToCell[3] = cell.cellToCell[(neighborIdx - 2) % 4];
    //            }

    //        }
    //        //////////////////////////////////////////////////////////////////////////

    //        if( true )
    //        {
    //            // search for suitable cells

    //            std::vector<uint> candidates;

    //            std::function<void(uint)> findCandidates = [&](uint candidateIdx)
    //            {
    //                MeshCell candidateCell = this->cells[candidateIdx];

    //                bool sharesNode = false;

    //                for (auto faceNodeIdx : newFace.faceToNode)
    //                    for (auto candidateNodeIdx : candidateCell.cellToNode)
    //                        if (faceNodeIdx == candidateNodeIdx)
    //                            sharesNode = true;

    //                if (!sharesNode) return;

    //                candidates.push_back(candidateIdx);

    //                for (auto candidateCellIdx : candidateCell.cellToCell)
    //                    if (candidateCellIdx != INVALID_IDX)
    //                        if (std::find(candidates.begin(), candidates.end(), candidateCellIdx) == candidates.end())
    //                            findCandidates(candidateCellIdx);
    //            };

    //            findCandidates(cellIdx);

    //            // https://stackoverflow.com/questions/3385229/c-erase-vector-element-by-value-rather-than-by-position
    //            candidates.erase(std::remove(candidates.begin(), candidates.end(), cellIdx), candidates.end());
    //            candidates.erase(std::remove(candidates.begin(), candidates.end(), neighborCellIdx), candidates.end());

    //            //////////////////////////////////////////////////////////////////////////

    //            newFace.faceCenter.x = 0.5 * (this->nodes[newFace.faceToNode[1]].x + this->nodes[newFace.faceToNode[0]].x);
    //            newFace.faceCenter.y = 0.5 * (this->nodes[newFace.faceToNode[1]].y + this->nodes[newFace.faceToNode[0]].y);

    //            // choose cells

    //            if (candidates.size() <= 4) {

    //                uint j = 2;
    //                for (uint i = 0; i < candidates.size(); i++) {
    //                    newFace.faceToCell[j++] = candidates[i];
    //                }
    //                while (j < 6) {
    //                    newFace.faceToCell[j++] = INVALID_IDX;
    //                }
    //            }
    //            else
    //            {
    //                uint numberOfPermutations = pow(candidates.size(), (candidates.size() - 4));

    //                uint   minimalPermutation = INVALID_IDX;
    //                double minimalDeviation = 1.0e99;

    //                for (uint permutation = 0; permutation < numberOfPermutations; permutation++)
    //                {
    //                    std::vector<uint> jumpList;
    //                    uint tmp = permutation;
    //                    do
    //                    {
    //                        jumpList.push_back(tmp % candidates.size());
    //                        tmp = tmp / candidates.size();
    //                    } while (tmp > 0);

    //                    std::sort(jumpList.begin(), jumpList.end());
    //                    jumpList.erase(std::unique(jumpList.begin(), jumpList.end()), jumpList.end());

    //                    if (jumpList.size() != (candidates.size() - 4)) continue;

    //                    arraytypes::Vec2 centroid;

    //                    centroid.x += 1.0 / 6.0 * cell.cellCenter.x;
    //                    centroid.y += 1.0 / 6.0 * cell.cellCenter.y;
    //                    centroid.x += 1.0 / 6.0 * neighborCell.cellCenter.x;
    //                    centroid.y += 1.0 / 6.0 * neighborCell.cellCenter.y;

    //                    for (uint i = 0; i < candidates.size(); i++)
    //                    {
    //                        if (std::find(jumpList.begin(), jumpList.end(), i) == jumpList.end())
    //                        {
    //                            centroid.x += 1.0 / 6.0 * this->cells[candidates[i]].cellCenter.x;
    //                            centroid.y += 1.0 / 6.0 * this->cells[candidates[i]].cellCenter.y;
    //                        }
    //                    }

    //                    double deviation = sqrt((newFace.faceCenter.x - centroid.x)
    //                        * (newFace.faceCenter.x - centroid.x)
    //                        + (newFace.faceCenter.y - centroid.y)
    //                        * (newFace.faceCenter.y - centroid.y));

    //                    if (deviation < minimalDeviation)
    //                    {
    //                        minimalPermutation = permutation;
    //                        minimalDeviation = deviation;
    //                    }
    //                }

    //                std::vector<uint> jumpList;
    //                uint tmp = minimalPermutation;
    //                do
    //                {
    //                    jumpList.push_back(tmp % candidates.size());
    //                    tmp = tmp / candidates.size();
    //                } while (tmp > 0);

    //                uint j = 2;
    //                for (uint i = 0; i < candidates.size(); i++)
    //                {
    //                    if (std::find(jumpList.begin(), jumpList.end(), i) == jumpList.end())
    //                    {
    //                        newFace.faceToCell[j++] = candidates[i];
    //                    }
    //                }
    //            }
    //        }

    //        //////////////////////////////////////////////////////////////////////////

    //        newFace.negCell = cellIdx;
    //        newFace.posCell = neighborCellIdx;

    //        //////////////////////////////////////////////////////////////////////////

    //        if ( cell.type == FLUID_CFF && neighborCell.type == FLUID_FCF ) newFace.negCell = cell.parent;
    //        if ( cell.type == FLUID_FCF && neighborCell.type == FLUID_CFF ) newFace.posCell = neighborCell.parent;

    //        //////////////////////////////////////////////////////////////////////////

    //        arraytypes::Vec2 faceVec( this->nodes[ newFace.faceToNode[ 1 ] ].x - this->nodes[ newFace.faceToNode[ 0 ] ].x,
    //                                  this->nodes[ newFace.faceToNode[ 1 ] ].y - this->nodes[ newFace.faceToNode[ 0 ] ].y );

    //        newFace.faceArea = sqrt( faceVec.x * faceVec.x + faceVec.y * faceVec.y );

    //        newFace.faceCenter.x = 0.5 * ( this->nodes[ newFace.faceToNode[ 1 ] ].x + this->nodes[ newFace.faceToNode[ 0 ] ].x );
    //        newFace.faceCenter.y = 0.5 * ( this->nodes[ newFace.faceToNode[ 1 ] ].y + this->nodes[ newFace.faceToNode[ 0 ] ].y );

    //        newFace.faceNormal.x =   faceVec.y / newFace.faceArea;
    //        newFace.faceNormal.y = - faceVec.x / newFace.faceArea;

    //        //////////////////////////////////////////////////////////////////////////

    //        this->faces.push_back( newFace );
    //    }
    //}
}

void GksMeshAdapter::sortFaces()
{
    real xMax = ( *std::max_element( this->faces.begin(), this->faces.end(), [this]( MeshFace lhs, MeshFace rhs ){ return lhs.faceCenter.x < rhs.faceCenter.x; } ) ).faceCenter.x;
    real yMax = ( *std::max_element( this->faces.begin(), this->faces.end(), [this]( MeshFace lhs, MeshFace rhs ){ return lhs.faceCenter.y < rhs.faceCenter.y; } ) ).faceCenter.y;
    real zMax = ( *std::max_element( this->faces.begin(), this->faces.end(), [this]( MeshFace lhs, MeshFace rhs ){ return lhs.faceCenter.z < rhs.faceCenter.z; } ) ).faceCenter.z;
    
    real xMin = ( *std::min_element( this->faces.begin(), this->faces.end(), [this]( MeshFace lhs, MeshFace rhs ){ return lhs.faceCenter.x < rhs.faceCenter.x; } ) ).faceCenter.x;
    real yMin = ( *std::min_element( this->faces.begin(), this->faces.end(), [this]( MeshFace lhs, MeshFace rhs ){ return lhs.faceCenter.y < rhs.faceCenter.y; } ) ).faceCenter.y;
    real zMin = ( *std::min_element( this->faces.begin(), this->faces.end(), [this]( MeshFace lhs, MeshFace rhs ){ return lhs.faceCenter.z < rhs.faceCenter.z; } ) ).faceCenter.z;

    real xRange = xMax - xMin;
    real yRange = yMax - yMin;
    real zRange = zMax - zMin;

    std::sort(  this->faces.begin(), this->faces.end(),
                [&,this]( MeshFace lhs, MeshFace rhs ){

                    if( lhs.level != rhs.level ) return lhs.level < rhs.level;

                    if( lhs.orientation != rhs.orientation ){
                        if      ( lhs.orientation == 'x' && rhs.orientation == 'y' ) return true;
                        else if ( lhs.orientation == 'y' && rhs.orientation == 'z' ) return true;
                        else if ( lhs.orientation == 'x' && rhs.orientation == 'z' ) return true;
                        else                                                         return false;
                    }

                    return ( ten * ten * lhs.faceCenter.z / zRange + ten * lhs.faceCenter.y / yRange + lhs.faceCenter.x / xRange ) < 
                           ( ten * ten * rhs.faceCenter.z / zRange + ten * rhs.faceCenter.y / yRange + rhs.faceCenter.x / xRange );
                }
             );
}

void GksMeshAdapter::countFaces()
{
    this->numberOfFacesPerLevelXYZ.resize( 3 * this->numberOfLevels );
    this->startOfFacesPerLevelXYZ.resize ( 3 * this->numberOfLevels );

    for( auto& i : this->numberOfFacesPerLevelXYZ ) i = 0;
    for( auto& i : this->startOfFacesPerLevelXYZ  ) i = 0;

    for( auto& face : this->faces ){
        if      ( face.orientation == 'x' ) this->numberOfFacesPerLevelXYZ[ 3 * face.level     ]++;
        else if ( face.orientation == 'y' ) this->numberOfFacesPerLevelXYZ[ 3 * face.level + 1 ]++;
        else if ( face.orientation == 'z' ) this->numberOfFacesPerLevelXYZ[ 3 * face.level + 2 ]++;
    }

    this->startOfFacesPerLevelXYZ[0] = 0;

    for( uint level = 1; level < 3 * this->numberOfLevels; level++ ){
        
        this->startOfFacesPerLevelXYZ[level] = this->startOfFacesPerLevelXYZ [level - 1]
                                             + this->numberOfFacesPerLevelXYZ[level - 1];
    }
}

void GksMeshAdapter::findSolidFaces()
{
    //for( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

    //    MeshFace& face = this->faces[ faceIdx ];

    //    if( this->cells[ face.posCell ].type == STOPPER_SOLID )
    //        this->solidFaces.push_back( faceIdx );
    //}
}

void GksMeshAdapter::generateInterfaceConnectivity()
{
    this->numberOfFineToCoarsePerLevel.resize( this->numberOfLevels );
    this->startOfFineToCoarsePerLevel.resize ( this->numberOfLevels );
    this->numberOfCoarseToFinePerLevel.resize( this->numberOfLevels );
    this->startOfCoarseToFinePerLevel.resize ( this->numberOfLevels );

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

        MeshCell& cell = this->cells[ cellIdx ];

        if( cell.type == FLUID_FCC ){

            uint_9 connectivity;

            connectivity[ 0 ] = cellIdx;
            connectivity[ 1 ] = cell.children[ 0 ];
            connectivity[ 2 ] = cell.children[ 1 ];
            connectivity[ 3 ] = cell.children[ 2 ];
            connectivity[ 4 ] = cell.children[ 3 ];
            connectivity[ 5 ] = cell.children[ 4 ];
            connectivity[ 6 ] = cell.children[ 5 ];
            connectivity[ 7 ] = cell.children[ 6 ];
            connectivity[ 8 ] = cell.children[ 7 ];

            this->fineToCoarse.push_back( connectivity );

            this->numberOfFineToCoarsePerLevel[ cell.level ]++;
        }

        if( cell.type == FLUID_CFC ){
            
            uint_15 connectivity;

            connectivity[  0 ] = cellIdx;

            connectivity[  1 ] = cell.cellToCell[ 0 ];
            connectivity[  2 ] = cell.cellToCell[ 1 ];
            connectivity[  3 ] = cell.cellToCell[ 2 ];
            connectivity[  4 ] = cell.cellToCell[ 3 ];
            connectivity[  5 ] = cell.cellToCell[ 4 ];
            connectivity[  6 ] = cell.cellToCell[ 5 ];

            connectivity[  7 ] = cell.children[ 0 ];
            connectivity[  8 ] = cell.children[ 1 ];
            connectivity[  9 ] = cell.children[ 2 ];
            connectivity[ 10 ] = cell.children[ 3 ];
            connectivity[ 11 ] = cell.children[ 4 ];
            connectivity[ 12 ] = cell.children[ 5 ];
            connectivity[ 13 ] = cell.children[ 6 ];
            connectivity[ 14 ] = cell.children[ 7 ];

            this->coarseToFine.push_back( connectivity );

            numberOfCoarseToFinePerLevel[ cell.level ]++;
        }
    }

        std::cout << numberOfCoarseToFinePerLevel[ 0 ] << " " << this->numberOfFineToCoarsePerLevel[ 0 ] << std::endl;
    
    this->startOfFineToCoarsePerLevel[0] = 0;
    this->startOfCoarseToFinePerLevel[0] = 0;

    for( uint level = 1; level < this->numberOfLevels; level++ ){
        
        this->startOfFineToCoarsePerLevel[level] = this->startOfFineToCoarsePerLevel [level - 1]
                                                 + this->numberOfFineToCoarsePerLevel[level - 1];
        
        this->startOfCoarseToFinePerLevel[level] = this->startOfCoarseToFinePerLevel [level - 1]
                                                 + this->numberOfCoarseToFinePerLevel[level - 1];
    }
}

void GksMeshAdapter::writeMeshVTK(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "writeMeshVTK( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "by MeshGenerator\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << nodes.size() << " float" << std::endl;

    for (auto node : nodes){
        file << node.x << " " << node.y << " " << node.z << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELLS " << this->cells.size() << " " << this->cells.size() * 9 << std::endl;

    for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        uint_8 nodes;
        for( auto& i : nodes ) i = INVALID_INDEX;

        nodes[0] = cell.cellToNode[7];//[ 6 ];
        nodes[1] = cell.cellToNode[3];//[ 5 ];
        nodes[2] = cell.cellToNode[1];//[ 2 ];
        nodes[3] = cell.cellToNode[5];//[ 1 ];
        nodes[4] = cell.cellToNode[6];//[ 4 ];
        nodes[5] = cell.cellToNode[2];//[ 7 ];
        nodes[6] = cell.cellToNode[0];//[ 0 ];
        nodes[7] = cell.cellToNode[4];//[ 3 ];

        file << 8 << " ";

        for( uint i = 0; i < 8; i++ ) file << nodes[i] << " ";

        file << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELL_TYPES " << this->cells.size() << std::endl;

    for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
        file << 12 << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////

    file << "\nCELL_DATA " << this->cells.size() << std::endl;

    file << "FIELD Label " << 4 << std::endl;

    //////////////////////////////////////////////////////////////////////////

    file << "CellIdx 1 " << this->cells.size() << " int" << std::endl;

    for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

        file << cellIdx << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "level 1 " << this->cells.size() << " int" << std::endl;

    for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        file << cell.level << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "type 1 " << this->cells.size() << " int" << std::endl;

    for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        file << (uint) cell.type << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "isGhostCell 1 " << this->cells.size() << " int" << std::endl;

    for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        file << (uint) cell.isGhostCell << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file.close();
}

void GksMeshAdapter::writeMeshFaceVTK(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "writeMeshFaceVTK( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "by MeshGenerator\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << nodes.size() << " float" << std::endl;

    for (auto node : nodes){
        file << node.x << " " << node.y << " " << node.z << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELLS " << this->faces.size() << " " << 5 * this->faces.size() << std::endl;

    for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

        file << "4 ";

        file << this->faces[ faceIdx ].faceToNode[0] << " ";
        file << this->faces[ faceIdx ].faceToNode[1] << " ";
        file << this->faces[ faceIdx ].faceToNode[2] << " ";
        file << this->faces[ faceIdx ].faceToNode[3] << " ";

        file << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELL_TYPES " << this->faces.size() << std::endl;

    for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){
        file << "9" << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////

    file << "\nCELL_DATA " << this->faces.size() << std::endl;

    file << "FIELD Label " << 3 << std::endl;

    //////////////////////////////////////////////////////////////////////////

    file << "FaceIdx 1 " << this->faces.size() << " int" << std::endl;

    for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

        file << faceIdx << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "level 1 " << this->faces.size() << " int" << std::endl;

    for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

        file << this->faces[ faceIdx ].level << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "orientation 1 " << this->faces.size() << " int" << std::endl;

    for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

        file << (int)this->faces[ faceIdx ].orientation << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    file << "VECTORS posCell double" << std::endl;

	for ( auto face : this->faces )
    {
        MeshCell& cell = this->cells[ face.posCell ];

        Vec3 vec = cell.cellCenter - face.faceCenter;
            
		file << vec.x << " ";
		file << vec.y << " ";
		file << vec.z << std::endl;
    }

    file << "VECTORS negCell double" << std::endl;

	for ( auto face : this->faces )
    {
        MeshCell& cell = this->cells[ face.negCell ];

        Vec3 vec = cell.cellCenter - face.faceCenter;
            
		file << vec.x << " ";
		file << vec.y << " ";
		file << vec.z << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file.close();
}

void GksMeshAdapter::writeMeshCellToCellVTK(std::string filename)
{
    //std::ofstream file;

    //file.open(filename);

    //file << "# vtk DataFile Version 3.0\n";
    //file << "by MeshGenerator\n";
    //file << "ASCII\n";
    //file << "DATASET UNSTRUCTURED_GRID\n";

    //file << "POINTS " << this->cells.size() << " float" << std::endl;

    //for (auto cell : cells){
    //    file << cell.cellCenter.x << " " << cell.cellCenter.y << " 0.0" << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file << "CELLS " << 8 * this->cells.size() << " " << 3 * 8 * this->cells.size() << std::endl;

    //for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    for( uint i = 0; i < 8; i++ )
    //        if(  this->cells[ cellIdx ].cellToCell[ i ] != INVALID_INDEX )
    //            file << "2 " << cellIdx << " " << this->cells[ cellIdx ].cellToCell[ i ] << " " << std::endl;
    //        else
    //            file << "2 " << cellIdx << " " << cellIdx << " " << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file << "CELL_TYPES " << 8 * this->cells.size() << std::endl;

    //for ( uint i = 0; i < 8 * this->cells.size(); i++ ){
    //    file << "3" << std::endl;
    //}
    ////////////////////////////////////////////////////////////////////////////

    //file << "\nCELL_DATA " << 8 * this->cells.size() << std::endl;

    //file << "FIELD Label " << 2 << std::endl;

    ////////////////////////////////////////////////////////////////////////////

    //file << "CellIdx 1 " << 8 * this->cells.size() << " int" << std::endl;

    //for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    for( uint i = 0; i < 8; i++ )
    //        file << cellIdx << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file << "CellToCell 1 " << 8 * this->cells.size() << " int" << std::endl;

    //for ( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){

    //    for( uint i = 0; i < 8; i++ )
    //        file << i << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file.close();
}

void GksMeshAdapter::writeMeshFaceToCellVTK(std::string filename)
{
    //std::ofstream file;

    //file.open(filename);

    //file << "# vtk DataFile Version 3.0\n";
    //file << "by MeshGenerator\n";
    //file << "ASCII\n";
    //file << "DATASET UNSTRUCTURED_GRID\n";

    //file << "POINTS " << this->cells.size() + this->faces.size() << " float" << std::endl;

    //for (auto cell : cells){
    //    file << cell.cellCenter.x << " " << cell.cellCenter.y << " 0.0" << std::endl;
    //}

    //for (auto face : faces){
    //    file << face.faceCenter.x << " " << face.faceCenter.y << " 0.0" << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file << "CELLS " << 6 * this->faces.size() << " " << 3 * 6 * this->faces.size() << std::endl;

    //for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

    //    for( uint i = 0; i < 6; i++ )
    //        if(  this->faces[ faceIdx ].faceToCell[ i ] != INVALID_INDEX )
    //            file << "2 " << this->cells.size() + faceIdx << " " << this->faces[ faceIdx ].faceToCell[ i ] << " " << std::endl;
    //        else
    //            file << "2 " << this->cells.size() + faceIdx << " " << this->cells.size() + faceIdx << " " << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file << "CELL_TYPES " << 6 * this->faces.size() << std::endl;

    //for ( uint i = 0; i < 6 * this->faces.size(); i++ ){
    //    file << "3" << std::endl;
    //}
    ////////////////////////////////////////////////////////////////////////////

    //file << "\nCELL_DATA " << 6 * this->faces.size() << std::endl;

    //file << "FIELD Label " << 2 << std::endl;

    ////////////////////////////////////////////////////////////////////////////

    //file << "FaceIdx 1 " << 6 * this->faces.size() << " int" << std::endl;

    //for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

    //    for( uint i = 0; i < 6; i++ )
    //        file << faceIdx << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file << "FaceToCell 1 " << 6 * this->faces.size() << " int" << std::endl;

    //for ( uint faceIdx = 0; faceIdx < this->faces.size(); faceIdx++ ){

    //    for( uint i = 0; i < 6; i++ )
    //        file << i << std::endl;
    //}

    ////////////////////////////////////////////////////////////////////////////

    //file.close();
}

double GksMeshAdapter::getDx(uint level)
{
    return dxCoarse / pow( 2.0, level );
}
