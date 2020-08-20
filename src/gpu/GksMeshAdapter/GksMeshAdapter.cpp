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
//! \file GksMeshAdapter.cpp
//! \ingroup GksMeshAdapter
//! \author Stephan Lenz
//=======================================================================================
#include "GksMeshAdapter.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <fstream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Core/Logger/Logger.h"

#include "GridGenerator/grid/distributions/D3Q27.h"
#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "GridGenerator/grid/NodeValues.h"
#include "GridGenerator/utilities/math/Math.h"

#include "MeshCell.h"
#include "MeshFace.h"

GksMeshAdapter::GksMeshAdapter(SPtr<MultipleGridBuilder> gridBuilder)
    : gridBuilder(gridBuilder)
{}

void GksMeshAdapter::inputGrid()
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "inputGrid()" << "\n";

    this->numberOfLevels = this->gridBuilder->getNumberOfGridLevels();

    std::vector< SPtr<Grid> > grids = this->gridBuilder->getGrids();

    this->dxCoarse = grids[0]->getDelta();

    //////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Allocate gridToMesh[][]" << "\n";

    this->gridToMesh.resize( this->gridBuilder->getNumberOfGridLevels() );

    for( uint level = 0; level < this->gridBuilder->getNumberOfGridLevels(); level++ ){
        this->gridToMesh[level].resize( grids[level]->getSize() );

        for( auto& cellIdx : this->gridToMesh[level] ) cellIdx = INVALID_INDEX;
    }

    //////////////////////////////////////////////////////////////////////////
    //
    //    I d e n t i f y    C e l l s    i n    L B - G r i d
    //
    //////////////////////////////////////////////////////////////////////////

    uint numberOfCells = 0;

    for( uint level = 0; level < this->gridBuilder->getNumberOfGridLevels(); level++ ){
        for( uint gridIdx = 0; gridIdx < grids[level]->getSize(); gridIdx++ ){
            if (grids[level]->getFieldEntry(gridIdx)  != STOPPER_COARSE_UNDER_FINE &&
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

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Allocate " << numberOfCells << " cells" << "\n";

    this->cells.resize( numberOfCells );

    for( uint level = 0; level < this->gridBuilder->getNumberOfGridLevels(); level++ ){
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

    this->findCellToCellConnectivity();
    this->countCells();
    this->partitionCells();
    this->generateNodes();
    this->computeCellGeometry();

    this->generateFaces();
    this->sortFaces();
    this->countFaces();

    //////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "inputGrid() finished!" << "\n";
}

void GksMeshAdapter::findCellToCellConnectivity()
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "findCellToCellConnectivity()" << "\n";

    std::vector< SPtr<Grid> > grids = this->gridBuilder->getGrids();

    Distribution dirs = DistributionHelper::getDistribution27();

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        real x, y, z;
        grids[cell.level]->transIndexToCoords(cell.gridIdx, x, y, z);

        real d = grids[cell.level]->getDelta();

        for( uint idx = 0; idx < 27; idx++ )
        {
            if( idx == DIR_27_ZERO ) continue;

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
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "countCells()" << "\n";

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
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "partitionCells()" << "\n";

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
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "refreshCellConnectivity()" << "\n";

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

void GksMeshAdapter::generateNodes()
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "generateNodes()" << "\n";

    std::vector< SPtr<Grid> > grids = gridBuilder->getGrids();

    nodes.reserve( 2 * this->cells.size() );

    Distribution dirs = DistributionHelper::getDistribution27();

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        //if( cell.type == STOPPER_SOLID ) continue;

        real x, y, z;
        grids[cell.level]->transIndexToCoords(cell.gridIdx, x, y, z);

        real d = 0.5 * grids[cell.level]->getDelta();

        std::array<Vec3,8> dir;

        for( uint idx = 0; idx < 8; idx++ )
        {
            if( cell.cellToNode[idx] == INVALID_INDEX )
            {

                real dx = dirs.directions[idx + 19][0] * d;
                real dy = dirs.directions[idx + 19][1] * d;
                real dz = dirs.directions[idx + 19][2] * d;

                nodes.push_back( Vec3( x + dx, y + dy, z + dz ) );

                cell.cellToNode[idx] = nodes.size()-1;

                //// register new node at neighbor cells on same level
                for (uint idx = 0; idx < 8; idx++)
                {
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
    }
}

void GksMeshAdapter::computeCellGeometry()
{    
    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
        
        MeshCell& cell = this->cells[ cellIdx ];

        Vec3 cellCenter;

        for( uint node = 0; node < 8; node++ ){
            cellCenter = cellCenter + this->nodes[ cell.cellToNode[node] ];
        }

        cell.cellCenter.x = cellCenter.x / c8o1;
        cell.cellCenter.y = cellCenter.y / c8o1;
        cell.cellCenter.z = cellCenter.z / c8o1;
    }
}

void GksMeshAdapter::generateFaces()
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "generateFaces()" << "\n";

    std::vector< SPtr<Grid> > grids = this->gridBuilder->getGrids();

    this->faces.reserve( 2 * this->cells.size() );

    for( uint cellIdx = 0; cellIdx < this->cells.size(); cellIdx++ ){
    
        MeshCell& cell = this->cells[ cellIdx ];

        //if( cell.type == BC_SOLID || cell.type == STOPPER_SOLID ) continue;

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
            
            Vec3 faceCenter;

            for( uint node = 0; node < 4; node++ ){
                faceCenter = faceCenter + this->nodes[ newFace.faceToNode[node] ];
            }

            newFace.faceCenter.x = faceCenter.x / c4o1;
            newFace.faceCenter.y = faceCenter.y / c4o1;
            newFace.faceCenter.z = faceCenter.z / c4o1;

            this->faces.push_back( newFace );
        }
    }
}

void GksMeshAdapter::sortFaces()
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // sort by level and orientation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "sortFaces()" << "\n";

    std::stable_sort(this->faces.begin(), this->faces.end(),
            [&, this](MeshFace lhs, MeshFace rhs)
            {
                if( lhs.level != rhs.level ) return lhs.level < rhs.level;

                if (lhs.orientation != rhs.orientation) {
                    if      (lhs.orientation == 'x' && rhs.orientation == 'y') return true;
                    else if (lhs.orientation == 'y' && rhs.orientation == 'z') return true;
                    else if (lhs.orientation == 'x' && rhs.orientation == 'z') return true;
                    else                                                       return false;
                }

                return false;
            }
    );

    this->countFaces();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // sort into blocks
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::array<char, 3> orientations = {'x', 'y', 'z'};

    for( uint level = 0; level < this->gridBuilder->getNumberOfLevels(); level++ )
    {
        for( uint idx = 0; idx < 3; idx++ )
        {
            uint start =         this->startOfFacesPerLevelXYZ [ 3 * level + idx];
            uint end   = start + this->numberOfFacesPerLevelXYZ[ 3 * level + idx];

            real xMax = (*std::max_element(this->faces.begin() + start, this->faces.begin() + end, [this](MeshFace lhs, MeshFace rhs) { return lhs.faceCenter.x < rhs.faceCenter.x; })).faceCenter.x;
            real yMax = (*std::max_element(this->faces.begin() + start, this->faces.begin() + end, [this](MeshFace lhs, MeshFace rhs) { return lhs.faceCenter.y < rhs.faceCenter.y; })).faceCenter.y;
            real zMax = (*std::max_element(this->faces.begin() + start, this->faces.begin() + end, [this](MeshFace lhs, MeshFace rhs) { return lhs.faceCenter.z < rhs.faceCenter.z; })).faceCenter.z;

            real xMin = (*std::min_element(this->faces.begin() + start, this->faces.begin() + end, [this](MeshFace lhs, MeshFace rhs) { return lhs.faceCenter.x < rhs.faceCenter.x; })).faceCenter.x;
            real yMin = (*std::min_element(this->faces.begin() + start, this->faces.begin() + end, [this](MeshFace lhs, MeshFace rhs) { return lhs.faceCenter.y < rhs.faceCenter.y; })).faceCenter.y;
            real zMin = (*std::min_element(this->faces.begin() + start, this->faces.begin() + end, [this](MeshFace lhs, MeshFace rhs) { return lhs.faceCenter.z < rhs.faceCenter.z; })).faceCenter.z;

            real xRange = xMax - xMin;
            real yRange = yMax - yMin;
            real zRange = zMax - zMin;

            uint blockDim = 8;

            real dx = this->gridBuilder->getGrid(level)->getDelta();

            std::sort(this->faces.begin() + start, this->faces.begin() + end,
                [&, this](MeshFace lhs, MeshFace rhs)
            {
                uint xIdxLhs = lround((lhs.faceCenter.x - xMin) / dx);
                uint yIdxLhs = lround((lhs.faceCenter.y - yMin) / dx);
                uint zIdxLhs = lround((lhs.faceCenter.z - zMin) / dx);

                uint xIdxRhs = lround((rhs.faceCenter.x - xMin) / dx);
                uint yIdxRhs = lround((rhs.faceCenter.y - yMin) / dx);
                uint zIdxRhs = lround((rhs.faceCenter.z - zMin) / dx);

                real xBlockLhs = xIdxLhs / blockDim;
                real yBlockLhs = yIdxLhs / blockDim;
                real zBlockLhs = zIdxLhs / blockDim;

                real xBlockRhs = xIdxRhs / blockDim;
                real yBlockRhs = yIdxRhs / blockDim;
                real zBlockRhs = zIdxRhs / blockDim;

                if (zBlockLhs < zBlockRhs) return true;
                if (zBlockLhs > zBlockRhs) return false;
                if (yBlockLhs < yBlockRhs) return true;
                if (yBlockLhs > yBlockRhs) return false;
                if (xBlockLhs < xBlockRhs) return true;
                if (xBlockLhs > xBlockRhs) return false;

                if (zIdxLhs < zIdxRhs) return true;
                if (zIdxLhs > zIdxRhs) return false;
                if (yIdxLhs < yIdxRhs) return true;
                if (yIdxLhs > yIdxRhs) return false;
                if (xIdxLhs < xIdxRhs) return true;
                if (xIdxLhs > xIdxRhs) return false;

                return true;
            }
            );
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // partition by inner and out for communication hiding
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    this->numberOfInnerFacesPerLevel.resize( this->numberOfLevels );

    for( uint level = 0; level < this->gridBuilder->getNumberOfLevels(); level++ )
    {
        auto bound =
        std::stable_partition(  this->faces.begin() + this->startOfFacesPerLevelXYZ [3 * level], 
                                this->faces.begin() + this->startOfFacesPerLevelXYZ [3 * level] 
                                                    + this->numberOfFacesPerLevelXYZ[3 * level + 0] 
                                                    + this->numberOfFacesPerLevelXYZ[3 * level + 1] 
                                                    + this->numberOfFacesPerLevelXYZ[3 * level + 2], 
                                    [this](MeshFace& lhs)
                                    {
                                        for( uint neighborIndex = 0; neighborIndex < 6; neighborIndex++ )
                                        {
                                            uint neighborCellIndex = this->cells[ lhs.posCell ].cellToCell[ neighborIndex ];
                                            if( neighborCellIndex != INVALID_INDEX && this->cells[ neighborCellIndex ].isRecvCell )
                                            {
                                                return false;
                                            }
                                        }
                                        for( uint neighborIndex = 0; neighborIndex < 6; neighborIndex++ )
                                        {
                                            uint neighborCellIndex = this->cells[ lhs.negCell ].cellToCell[ neighborIndex ];
                                            if( neighborCellIndex != INVALID_INDEX && this->cells[ neighborCellIndex ].isRecvCell )
                                            {
                                                return false;
                                            }
                                        }

                                        return true;
                                    }
                                 );

        this->numberOfInnerFacesPerLevel[ level ] = 0;
        for( auto it = this->faces.begin() + this->startOfFacesPerLevelXYZ [3 * level]; it != bound; it++ )
        {
            this->numberOfInnerFacesPerLevel[ level ]++;
        }

        *logging::out << logging::Logger::INFO_LOW << "    Level " << level << ": " << this->numberOfFacesPerLevelXYZ[ 3 * level + 0 ]
                                                                                     + this->numberOfFacesPerLevelXYZ[ 3 * level + 1 ]
                                                                                     + this->numberOfFacesPerLevelXYZ[ 3 * level + 2 ]
                                                                                    << " faces" << "\n";
        *logging::out << logging::Logger::INFO_LOW << "    Level " << level << ": " << this->numberOfInnerFacesPerLevel[ level ] << " inner faces" << "\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

void GksMeshAdapter::findPeriodicBoundaryNeighbors()
{
    for( uint level = 0; level < this->numberOfLevels; level++ )
    {
        SPtr<Grid> grid = this->gridBuilder->getGrid(level);

        if( !grid->getPeriodicityX() && !grid->getPeriodicityY() && !grid->getPeriodicityZ() )
            throw std::runtime_error( "GksMeshAdapter::findPeriodicBoundaryNeighbors() failed, because no periodic direction is set!" );

        uint startIdx = startOfCellsPerLevel[ level ] + numberOfBulkCellsPerLevel[ level ];

        uint endIdx   = startOfCellsPerLevel[ level ] + numberOfCellsPerLevel[ level ];

        for( uint cellIdx = startIdx; cellIdx < endIdx; cellIdx++ )
        {
            MeshCell cell = this->cells[ cellIdx ];

            if( cell.type != STOPPER_OUT_OF_GRID && cell.type != STOPPER_OUT_OF_GRID_BOUNDARY && cell.type != STOPPER_SOLID ) continue;

            Vec3 gridStart ( grid->getStartX() + c1o2 * grid->getDelta(),
                             grid->getStartY() + c1o2 * grid->getDelta(),
                             grid->getStartZ() + c1o2 * grid->getDelta() );

            Vec3 gridEnd   ( grid->getEndX()   - c1o2 * grid->getDelta(),
                             grid->getEndY()   - c1o2 * grid->getDelta(),
                             grid->getEndZ()   - c1o2 * grid->getDelta() );

            Vec3 size = gridEnd - gridStart;

            Vec3 delta;

            if( grid->getPeriodicityX() && cell.cellCenter.x < gridStart.x ) delta.x =   size.x;
            if( grid->getPeriodicityX() && cell.cellCenter.x > gridEnd.x   ) delta.x = - size.x;

            if( grid->getPeriodicityY() && cell.cellCenter.y < gridStart.y ) delta.y =   size.y;
            if( grid->getPeriodicityY() && cell.cellCenter.y > gridEnd.y   ) delta.y = - size.y;

            if( grid->getPeriodicityZ() && cell.cellCenter.z < gridStart.z ) delta.z =   size.z;
            if( grid->getPeriodicityZ() && cell.cellCenter.z > gridEnd.z   ) delta.z = - size.z;

            uint neighborGridIdx = grid->transCoordToIndex( cell.cellCenter.x + delta.x,
                                                            cell.cellCenter.y + delta.y,
                                                            cell.cellCenter.z + delta.z );
            
            if( neighborGridIdx == INVALID_INDEX ) throw std::runtime_error( std::string("No periodic cell found!") );

            uint neighborIdx = this->gridToMesh[ level ][ neighborGridIdx ];

            if( neighborIdx == INVALID_INDEX )
            {
                std::stringstream s;

                s << "No periodic cell found: ";
                s << "( " << cell.cellCenter.x           << ", " << cell.cellCenter.y           << ", " << cell.cellCenter.z           << " )";
                s << "( " << cell.cellCenter.x + delta.x << ", " << cell.cellCenter.y + delta.y << ", " << cell.cellCenter.z + delta.z << " )";

                s << "( " << gridStart.x << ", " << gridStart.y << ", " << gridStart.z << " )";
                s << "( " << gridEnd.x   << ", " << gridEnd.y   << ", " << gridEnd.z   << " )";

                throw std::runtime_error( s.str() );
            }
            
            this->periodicBoundaryNeighbors.push_back( {cellIdx, neighborIdx} );
        }
    }
}

double GksMeshAdapter::getDx(uint level)
{
    return dxCoarse / pow( 2.0, level );
}
