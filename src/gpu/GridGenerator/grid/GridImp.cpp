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
//! \file GridImp.cpp
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz, Martin Schoenherr
//=======================================================================================
#include "GridImp.h"

#include <iostream>
#include <omp.h>
#include <sstream>
# include <algorithm>
#include <cmath>

#include "global.h"

#include "geometries/Object.h"
#include "geometries/Vertex/Vertex.h"
#include "geometries/Triangle/Triangle.h"
#include "geometries/TriangularMesh/TriangularMesh.h"
#include "geometries/TriangularMesh/TriangularMeshStrategy.h"
#include "geometries/BoundingBox/BoundingBox.h"

#include "grid/distributions/Distribution.h"
#include "grid/Field.h"
#include "grid/GridInterface.h"
#include "grid/NodeValues.h"

#include "io/GridVTKWriter/GridVTKWriter.h"

#include "utilities/communication.h"
#include "utilities/math/Math.h"

int DIRECTIONS[DIR_END_MAX][DIMENSION];

using namespace vf::gpu;

GridImp::GridImp(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, Distribution distribution, uint level) 
            : object(object), 
    startX(startX),
    startY(startY),
    startZ(startZ),
    endX(endX),
    endY(endY),
    endZ(endZ),
    delta(delta),
    distribution(distribution),
    level(level),
    periodicityX(false),
    periodicityY(false),
    periodicityZ(false),
    enableFixRefinementIntoTheWall(false),
    gridInterface(nullptr),
    neighborIndexX(nullptr),
    neighborIndexY(nullptr),
    neighborIndexZ(nullptr),
    neighborIndexNegative(nullptr),
    sparseIndices(nullptr),
    qIndices(nullptr),
    qValues(nullptr),
    qPatches(nullptr),
    innerRegionFromFinerGrid(false),
    numberOfLayers(0),
    qComputationStage(qComputationStageType::ComputeQs)
{
    initalNumberOfNodesAndSize();
}

SPtr<GridImp> GridImp::makeShared(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, std::string d3Qxx, uint level)
{
    Distribution distribution = DistributionHelper::getDistribution(d3Qxx);
    SPtr<GridImp> grid(new GridImp(object, startX, startY, startZ, endX, endY, endZ, delta, distribution, level));
    return grid;
}


void GridImp::initalNumberOfNodesAndSize()
{
    const real length = endX - startX;
    const real width = endY - startY;
    const real height = endZ - startZ;

    nx = std::lround((length + delta) / delta);
    ny = std::lround((width + delta) / delta);
    nz = std::lround((height + delta) / delta);

    this->size = nx * ny * nz;
    this->sparseSize = size;
    distribution.setSize(size);
}

void GridImp::inital(const SPtr<Grid> fineGrid, uint numberOfLayers)
{
    field = Field(size);
    field.allocateMemory();

    this->neighborIndexX        = new int[this->size];
    this->neighborIndexY        = new int[this->size];
    this->neighborIndexZ        = new int[this->size];
    this->neighborIndexNegative = new int[this->size];

    this->sparseIndices = new int[this->size];

    this->qIndices = new uint[this->size];
    for (uint i = 0; i < this->size; i++)
        this->qIndices[i] = INVALID_INDEX;

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start initalNodesToOutOfGrid()\n";
#pragma omp parallel for
    for (int index = 0; index < (int)this->size; index++)
        this->initalNodeToOutOfGrid(index);
    
    if( this->innerRegionFromFinerGrid ){
        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start setInnerBasedOnFinerGrid()\n";
        this->setInnerBasedOnFinerGrid(fineGrid);
    }
    else{
        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start findInnerNodes()\n";
        this->object->findInnerNodes( shared_from_this() );
    }

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start addOverlap()\n";
    this->addOverlap();
    
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start fixOddCells()\n";
#pragma omp parallel for
    for (int index = 0; index < (int)this->size; index++)
        this->fixOddCell(index);
    
    if( enableFixRefinementIntoTheWall )
    {
        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start fixRefinementIntoWall()\n";
#pragma omp parallel for
        for (int xIdx = 0; xIdx < (int)this->nx; xIdx++) {
            for (uint yIdx = 0; yIdx < this->ny; yIdx++) {
                this->fixRefinementIntoWall( xIdx, yIdx, 0           ,  3 );
                this->fixRefinementIntoWall( xIdx, yIdx, this->nz - 1, -3 );
            }
        }

#pragma omp parallel for
        for (int xIdx = 0; xIdx < (int)this->nx; xIdx++) {
            for (uint zIdx = 0; zIdx < this->nz; zIdx++) {
                this->fixRefinementIntoWall( xIdx, 0           , zIdx,  2 );
                this->fixRefinementIntoWall( xIdx, this->ny - 1, zIdx, -2 );
            }
        }

#pragma omp parallel for
        for (int yIdx = 0; yIdx < (int)this->ny; yIdx++) {
            for (uint zIdx = 0; zIdx < this->nz; zIdx++) {
                this->fixRefinementIntoWall( 0           , yIdx, zIdx,  1 );
                this->fixRefinementIntoWall( this->nx - 1, yIdx, zIdx, -1 );
            }
        }
    }
    
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start findEndOfGridStopperNodes()\n";
#pragma omp parallel for
    for (int index = 0; index < (int)this->size; index++)
        this->findEndOfGridStopperNode(index);
    
    *logging::out << logging::Logger::INFO_INTERMEDIATE
        << "Grid created: " << "from (" << this->startX << ", " << this->startY << ", " << this->startZ << ") to (" << this->endX << ", " << this->endY << ", " << this->endZ << ")\n"
        << "nodes: " << this->nx << " x " << this->ny << " x " << this->nz << " = " << this->size << "\n";
}

void GridImp::setOddStart(bool xOddStart, bool yOddStart, bool zOddStart)
{
    this->xOddStart = xOddStart;
    this->yOddStart = yOddStart;
    this->zOddStart = zOddStart;
}

void GridImp::initalNodeToOutOfGrid(uint index) {
    this->field.setFieldEntryToInvalidOutOfGrid(index);
}

void GridImp::freeMemory()
{
    if( this->neighborIndexX        != nullptr ) { delete[] this->neighborIndexX;        this->neighborIndexX        = nullptr; }
    if( this->neighborIndexY        != nullptr ) { delete[] this->neighborIndexY;        this->neighborIndexY        = nullptr; }
    if( this->neighborIndexZ        != nullptr ) { delete[] this->neighborIndexZ;        this->neighborIndexZ        = nullptr; }
    if( this->neighborIndexNegative != nullptr ) { delete[] this->neighborIndexNegative; this->neighborIndexNegative = nullptr; }
    if( this->sparseIndices         != nullptr ) { delete[] this->sparseIndices;         this->sparseIndices         = nullptr; }
	if( this->qIndices              != nullptr ) { delete[] this->qIndices;              this->qIndices              = nullptr; }
	if( this->qValues               != nullptr ) { delete[] this->qValues;               this->qValues               = nullptr; }
	if( this->qPatches              != nullptr ) { delete[] this->qPatches;              this->qPatches              = nullptr; }

    field.freeMemory();
}

void GridImp::findInnerNodes()
{
#pragma omp parallel for
    for (int index = 0; index < (int)this->size; index++)
        this->findInnerNode(index);
}

void GridImp::findInnerNode(uint index)
{
    this->sparseIndices[index] = index;

    if( this->level != 0 ){
        const Cell cell = getOddCellFromIndex(index);
        if (isInside(cell))
            this->field.setFieldEntryToFluid(index);
    }
    else{
        real x, y, z;
        this->transIndexToCoords(index, x, y, z);
        const uint xIndex = getXIndex(x);
        const uint yIndex = getYIndex(y);
        const uint zIndex = getZIndex(z);

        if( xIndex != 0 && xIndex != this->nx-1 &&
            yIndex != 0 && yIndex != this->ny-1 &&
            zIndex != 0 && zIndex != this->nz-1 )
            this->field.setFieldEntryToFluid(index);
    }
}

void GridImp::discretize(Object* solidObject, char innerType, char outerType)
{
#pragma omp parallel for
    for (int index = 0; index < (int)this->size; index++)
    {
        this->sparseIndices[index] = index;

        if( this->getFieldEntry(index) == innerType ) continue;
        
        real x, y, z;
        this->transIndexToCoords(index, x, y, z);

        if( solidObject->isPointInObject(x, y, z, 0.0, 0.0) )
            this->setFieldEntry(index, innerType);
        //else
        //    this->setFieldEntry(index, outerType);
    }
}

bool GridImp::isInside(const Cell& cell) const
{
    return object->isCellInObject(cell);
}

////TODO: check where the fine grid starts (0.25 or 0.75) and if even or odd-cell is needed
// Cell numbering:
//       even start                            odd start
//    +---------+                           +---------+
//    |       +-----+-----+-----+           | +-----+-----+-----+
//    |       | |   |     |     |           | |     | |   |     |
//    |       +-----+-----+-----+           | +-----+-----+-----+
//    +---------+                           +---------+
//               0     1     2                   0     1     2
//              even      even                        even     
//                   odd                        odd         odd
//
Cell GridImp::getOddCellFromIndex(uint index) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    const uint xIndex = getXIndex(x);
    const uint yIndex = getYIndex(y);
    const uint zIndex = getZIndex(z);

    real xCellStart;
    if( this->xOddStart ) xCellStart = xIndex % 2 != 0 ? x - this->delta : x;
    else                  xCellStart = xIndex % 2 != 0 ? x               : x - this->delta;

    real yCellStart;
    if( this->yOddStart ) yCellStart = yIndex % 2 != 0 ? y - this->delta : y;
    else                  yCellStart = yIndex % 2 != 0 ? y               : y - this->delta;

    real zCellStart;
    if( this->zOddStart ) zCellStart = zIndex % 2 != 0 ? z - this->delta : z;
    else                  zCellStart = zIndex % 2 != 0 ? z               : z - this->delta;

    return Cell(xCellStart, yCellStart, zCellStart, delta);
}

void GridImp::setInnerBasedOnFinerGrid(const SPtr<Grid> fineGrid)
{
    for( uint index = 0; index < this->size; index++ ){

        real x, y, z;
        this->transIndexToCoords(index, x, y, z);

        uint childIndex[8];

        childIndex[0] = fineGrid->transCoordToIndex( x + 0.25 * this->delta, y + 0.25 * this->delta, z + 0.25 * this->delta );
        childIndex[1] = fineGrid->transCoordToIndex( x + 0.25 * this->delta, y + 0.25 * this->delta, z - 0.25 * this->delta );
        childIndex[2] = fineGrid->transCoordToIndex( x + 0.25 * this->delta, y - 0.25 * this->delta, z + 0.25 * this->delta );
        childIndex[3] = fineGrid->transCoordToIndex( x + 0.25 * this->delta, y - 0.25 * this->delta, z - 0.25 * this->delta );
        childIndex[4] = fineGrid->transCoordToIndex( x - 0.25 * this->delta, y + 0.25 * this->delta, z + 0.25 * this->delta );
        childIndex[5] = fineGrid->transCoordToIndex( x - 0.25 * this->delta, y + 0.25 * this->delta, z - 0.25 * this->delta );
        childIndex[6] = fineGrid->transCoordToIndex( x - 0.25 * this->delta, y - 0.25 * this->delta, z + 0.25 * this->delta );
        childIndex[7] = fineGrid->transCoordToIndex( x - 0.25 * this->delta, y - 0.25 * this->delta, z - 0.25 * this->delta );

        for( uint i = 0; i < 8; i++ ){
            if( childIndex[i] != INVALID_INDEX && fineGrid->getFieldEntry( childIndex[i] ) == FLUID ){
                this->setFieldEntry(index, FLUID);
                break;
            }
        }
    }
}

void GridImp::addOverlap()
{
    for( uint layer = 0; layer < this->numberOfLayers; layer++ ){
#pragma omp parallel for
        for (int index = 0; index < (int)this->size; index++)
            this->setOverlapTmp(index);

#pragma omp parallel for
        for (int index = 0; index < (int)this->size; index++)
            this->setOverlapFluid(index);
    }
}

void GridImp::setOverlapTmp( uint index )
{
    if( this->field.is( index, INVALID_OUT_OF_GRID ) ){
        
        if( this->hasNeighborOfType(index, FLUID) ){
            this->field.setFieldEntry( index, OVERLAP_TMP );
        }
    }
}

void GridImp::setOverlapFluid( uint index )
{
    if( this->field.is( index, OVERLAP_TMP ) ){
        this->field.setFieldEntry( index, FLUID );
    }
}

void GridImp::fixRefinementIntoWall(uint xIndex, uint yIndex, uint zIndex, int dir)
{

    real x = this->startX + this->delta * xIndex;
    real y = this->startY + this->delta * yIndex;
    real z = this->startZ + this->delta * zIndex;

    uint index = this->transCoordToIndex(x, y, z);

    if( !this->xOddStart && ( dir == 1 || dir == -1 ) && ( xIndex % 2 == 1 || xIndex == 0 ) ) return;
    if( !this->yOddStart && ( dir == 2 || dir == -2 ) && ( yIndex % 2 == 1 || yIndex == 0 ) ) return;
    if( !this->zOddStart && ( dir == 3 || dir == -3 ) && ( zIndex % 2 == 1 || zIndex == 0 ) ) return;

    // Dont do this if inside of the domain
    if(  this->xOddStart && ( dir == 1 || dir == -1 ) && ( xIndex % 2 == 0 && xIndex != 0 ) ) return;
    if(  this->yOddStart && ( dir == 2 || dir == -2 ) && ( yIndex % 2 == 0 && yIndex != 0 ) ) return;
    if(  this->zOddStart && ( dir == 3 || dir == -3 ) && ( zIndex % 2 == 0 && zIndex != 0 ) ) return;
    
    //////////////////////////////////////////////////////////////////////////

    real dx{ 0.0 }, dy{ 0.0 }, dz{ 0.0 };

    if      ( dir ==  1 ){ dx =   this->delta; dy = 0.0;           dz = 0.0;           }
    else if ( dir == -1 ){ dx = - this->delta; dy = 0.0;           dz = 0.0;           }
    else if ( dir ==  2 ){ dx = 0.0;           dy =   this->delta; dz = 0.0;           }
    else if ( dir == -2 ){ dx = 0.0;           dy = - this->delta; dz = 0.0;           }
    else if ( dir ==  3 ){ dx = 0.0;           dy = 0.0;           dz =   this->delta; }
    else if ( dir == -3 ){ dx = 0.0;           dy = 0.0;           dz = - this->delta; }

    //////////////////////////////////////////////////////////////////////////

    char type = this->field.getFieldEntry(index);

    char type2    = ( type == FLUID ) ? ( INVALID_OUT_OF_GRID ) : ( FLUID );
    uint distance = ( type == FLUID ) ? ( 9                   ) : ( 5     );

    bool allTypesAreTheSame = true;

    for( uint i = 1; i <= distance; i++ ){
        uint neighborIndex = this->transCoordToIndex(x + i * dx, y + i * dy, z + i * dz);

        if( neighborIndex != INVALID_INDEX && !this->field.is( neighborIndex, type ) )
            allTypesAreTheSame = false;
    }

    //////////////////////////////////////////////////////////////////////////

    if( allTypesAreTheSame )
        return;

    this->setFieldEntry(index, type2);

    for( uint i = 1; i <= distance; i++ ){
        uint neighborIndex = this->transCoordToIndex(x + i * dx, y + i * dy, z + i * dz);

        this->setFieldEntry(neighborIndex, type2);
    }
}

void GridImp::findStopperNode(uint index) // deprecated
{
    if(isValidEndOfGridStopper(index))
        this->field.setFieldEntryToStopperOutOfGrid(index);

    if (isValidSolidStopper(index))
        this->field.setFieldEntry(index, STOPPER_SOLID);
}

void GridImp::findEndOfGridStopperNode(uint index)
{
	if (isValidEndOfGridStopper(index)){
        if( this->level != 0 )
		    this->field.setFieldEntryToStopperOutOfGrid(index);
        else
            this->field.setFieldEntryToStopperOutOfGridBoundary(index);
    }
    
	if (isValidEndOfGridBoundaryStopper(index))
		this->field.setFieldEntryToStopperOutOfGridBoundary(index);
}

void GridImp::findSolidStopperNode(uint index)
{
	if (isValidSolidStopper(index))
		this->field.setFieldEntry(index, STOPPER_SOLID);
}

void GridImp::findBoundarySolidNode(uint index)
{
	if (shouldBeBoundarySolidNode(index)) 
	{
		this->field.setFieldEntry(index, BC_SOLID);
		this->qIndices[index] = this->numberOfSolidBoundaryNodes++;
		//grid->setNumberOfSolidBoundaryNodes(grid->getNumberOfSolidBoundaryNodes() + 1);
	}
}

void GridImp::fixOddCell(uint index)
{
    Cell cell = getOddCellFromIndex(index);
    if (isOutSideOfGrid(cell))
        return;
    if (contains(cell, FLUID))
        setNodeTo(cell, FLUID);
}

bool GridImp::isOutSideOfGrid(Cell &cell) const
{
    for (const auto point : cell) {
        if (point.x < startX || point.x > endX
            || point.y < startY || point.y > endY
            || point.z < startZ || point.z > endZ)
            return true;
    }
    return false;
}

bool GridImp::contains(Cell &cell, char type) const
{
    for (const auto point : cell) {
		uint index = transCoordToIndex(point.x, point.y, point.z);
		if (index == INVALID_INDEX)
			continue;
        if (field.is(index, type))
            return true;
    }
    return false;
}

bool GridImp::cellContainsOnly(Cell &cell, char type) const
{
    for (const auto point : cell) {
		uint index = transCoordToIndex(point.x, point.y, point.z);
		if (index == INVALID_INDEX)
            return false;
        if (!field.is(index, type))
            return false;
    }
    return true;
}

bool GridImp::cellContainsOnly(Cell &cell, char typeA, char typeB) const
{
    for (const auto point : cell) {
		uint index = transCoordToIndex(point.x, point.y, point.z);
		if (index == INVALID_INDEX)
            return false;
        if (!field.is(index, typeA) && !field.is(index, typeB))
            return false;
    }
    return true;
}

const Object * GridImp::getObject() const
{
    return this->object;
}

void GridImp::setNodeTo(Cell &cell, char type)
{
    for (const auto point : cell) {
		uint index = transCoordToIndex(point.x, point.y, point.z);
		if (index == INVALID_INDEX)
			continue;
		field.setFieldEntry(index, type);
    }
}

void GridImp::setNodeTo(uint index, char type)
{
	if( index != INVALID_INDEX )
		field.setFieldEntry(index, type);
}

bool GridImp::isNode(uint index, char type) const
{
    if( index != INVALID_INDEX )
		return field.is(index, type);

    throw std::runtime_error("GridImp::isNode() -> index == INVALID_INDEX not supported.");
}

bool GridImp::isValidEndOfGridStopper(uint index) const
{
	// Lenz: also includes corner stopper nodes
	if (!this->field.is(index, INVALID_OUT_OF_GRID))
		return false;

	return hasNeighborOfType(index, FLUID);
}

bool GridImp::isValidEndOfGridBoundaryStopper(uint index) const
{
	// Lenz: also includes corner stopper nodes
	if (!this->field.is(index, FLUID))
		return false;

	return ! hasAllNeighbors(index);
}

bool GridImp::isValidSolidStopper(uint index) const
{
	// Lenz: also includes corner stopper nodes
	if (!this->field.is(index, INVALID_SOLID))
		return false;

	return hasNeighborOfType(index, FLUID);
}

bool GridImp::shouldBeBoundarySolidNode(uint index) const
{
	if (!this->field.is(index, FLUID))
		return false;

	return hasNeighborOfType(index, STOPPER_SOLID);
}

bool GridImp::hasAllNeighbors(uint index) const
{
	// new version by Lenz, utilizes the range based for loop for all directions
	real x, y, z;
	this->transIndexToCoords(index, x, y, z);
	for (const auto dir : this->distribution) {
		const uint neighborIndex = this->transCoordToIndex(x + dir[0] * this->getDelta(), y + dir[1] * this->getDelta(), z + dir[2] * this->getDelta());

		if (neighborIndex == INVALID_INDEX) return false;
	}

	return true;
}

bool GridImp::hasNeighborOfType(uint index, char type) const
{
	// new version by Lenz, utilizes the range based for loop for all directions
	real x, y, z;
	this->transIndexToCoords(index, x, y, z);
	for (const auto dir : this->distribution) {
		const uint neighborIndex = this->transCoordToIndex(x + dir[0] * this->getDelta(), y + dir[1] * this->getDelta(), z + dir[2] * this->getDelta());

		if (neighborIndex == INVALID_INDEX) continue;

		if (this->field.is(neighborIndex, type))
			return true;
	}

	return false;
}

bool GridImp::nodeInNextCellIs(int index, char type) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    const real neighborX = x + this->delta > endX ? endX : x + this->delta;
    const real neighborY = y + this->delta > endY ? endY : y + this->delta;
    const real neighborZ = z + this->delta > endZ ? endZ : z + this->delta;

    const uint indexX = transCoordToIndex(neighborX, y, z);
    const uint indexY = transCoordToIndex(x, neighborY, z);
    const uint indexZ = transCoordToIndex(x, y, neighborZ);

    const uint indexXY = transCoordToIndex(neighborX, neighborY, z);
    const uint indexYZ = transCoordToIndex(x, neighborY, neighborZ);
    const uint indexXZ = transCoordToIndex(neighborX, y, neighborZ);

    const uint indexXYZ = transCoordToIndex(neighborX, neighborY, neighborZ);

	const bool typeX   = indexX   == INVALID_INDEX ? false : this->field.is(indexX, type);
	const bool typeY   = indexY   == INVALID_INDEX ? false : this->field.is(indexY, type);
	const bool typeXY  = indexXY  == INVALID_INDEX ? false : this->field.is(indexXY, type);
	const bool typeZ   = indexZ   == INVALID_INDEX ? false : this->field.is(indexZ, type);
	const bool typeYZ  = indexYZ  == INVALID_INDEX ? false : this->field.is(indexYZ, type);
	const bool typeXZ  = indexXZ  == INVALID_INDEX ? false : this->field.is(indexXZ, type);
	const bool typeXYZ = indexXYZ == INVALID_INDEX ? false : this->field.is(indexXYZ, type);

    return typeX || typeY || typeXY || typeZ || typeYZ
        || typeXZ || typeXYZ;
}

bool GridImp::nodeInPreviousCellIs(int index, char type) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    const real neighborX = x - this->delta < startX ? startX : x - this->delta;
    const real neighborY = y - this->delta < startY ? startY : y - this->delta;
    const real neighborZ = z - this->delta < startZ ? startZ : z - this->delta;

    const uint indexX = transCoordToIndex(neighborX, y, z);
    const uint indexY = transCoordToIndex(x, neighborY, z);
    const uint indexZ = transCoordToIndex(x, y, neighborZ);

    const uint indexXY = transCoordToIndex(neighborX, neighborY, z);
    const uint indexYZ = transCoordToIndex(x, neighborY, neighborZ);
    const uint indexXZ = transCoordToIndex(neighborX, y, neighborZ);

    const uint indexXYZ = transCoordToIndex(neighborX, neighborY, neighborZ);

	const bool typeX   = indexX   == INVALID_INDEX ? false : this->field.is(indexX  , type);
	const bool typeY   = indexY   == INVALID_INDEX ? false : this->field.is(indexY  , type);
	const bool typeXY  = indexXY  == INVALID_INDEX ? false : this->field.is(indexXY , type);
	const bool typeZ   = indexZ   == INVALID_INDEX ? false : this->field.is(indexZ  , type);
	const bool typeYZ  = indexYZ  == INVALID_INDEX ? false : this->field.is(indexYZ , type);
	const bool typeXZ  = indexXZ  == INVALID_INDEX ? false : this->field.is(indexXZ , type);
	const bool typeXYZ = indexXYZ == INVALID_INDEX ? false : this->field.is(indexXYZ, type);

    return typeX || typeY || typeXY || typeZ || typeYZ
        || typeXZ || typeXYZ;
}

bool GridImp::nodeInCellIs(Cell& cell, char type) const
{
    for (const auto node : cell)
    {
        const uint index = transCoordToIndex(node.x, node.y, node.z);
		if (index == INVALID_INDEX)
			continue;
        if (field.is(index, type))
            return true;
    }
    return false;
}


void GridImp::setCellTo(uint index, char type)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    Cell cell(x, y, z, this->delta);
    for (const auto node : cell)
    {
        const uint nodeIndex = transCoordToIndex(node.x, node.y, node.z);
		if (nodeIndex == INVALID_INDEX)
			continue;
		this->field.setFieldEntry(nodeIndex, type);
    }
}


void GridImp::setNonStopperOutOfGridCellTo(uint index, char type)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    Cell cell(x, y, z, this->delta);
    for (const auto node : cell)
    {
        const uint nodeIndex = transCoordToIndex(node.x, node.y, node.z);
		if (nodeIndex == INVALID_INDEX)
			continue;

        if( this->getFieldEntry( nodeIndex ) != STOPPER_OUT_OF_GRID && 
            this->getFieldEntry( nodeIndex ) != STOPPER_OUT_OF_GRID_BOUNDARY )
            this->field.setFieldEntry(nodeIndex, type);
    }
}


void GridImp::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ)
{
    this->periodicityX = periodicityX;
    this->periodicityY = periodicityY;
    this->periodicityZ = periodicityZ;
}

void GridImp::setPeriodicityX(bool periodicity)
{
    this->periodicityX = periodicity;
}

void GridImp::setPeriodicityY(bool periodicity)
{
    this->periodicityY = periodicity;
}

void GridImp::setPeriodicityZ(bool periodicity)
{
    this->periodicityZ = periodicity;
}

bool GridImp::getPeriodicityX()
{
    return this->periodicityX;
}

bool GridImp::getPeriodicityY()
{
    return this->periodicityY;
}

bool GridImp::getPeriodicityZ()
{
    return this->periodicityZ;
}

void GridImp::setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall)
{
    this->enableFixRefinementIntoTheWall = enableFixRefinementIntoTheWall;
}

uint GridImp::transCoordToIndex(const real &x, const real &y, const real &z) const
{
    const uint xIndex = getXIndex(x);
    const uint yIndex = getYIndex(y);
    const uint zIndex = getZIndex(z);

	if (xIndex >= nx || yIndex >= ny || zIndex >= nz)
        return INVALID_INDEX;

    return xIndex + nx * (yIndex + ny * zIndex);
}

void GridImp::transIndexToCoords(uint index, real &x, real &y, real &z) const
{
    if (index == INVALID_INDEX)
        printf("Function: transIndexToCoords. GridImp Index: %d, size: %d. Exit Program!\n", index, size);

    x = (real)(index % nx);
    y = (real)((index / nx) % ny);
    z = (real)(((index / nx) / ny) % nz);

    x = (x * delta) + startX;
    y = (y * delta) + startY;
    z = (z * delta) + startZ;
}

uint GridImp::getLevel(real startDelta) const
{
    uint level = 0;
    real delta = this->delta;
    while(!vf::Math::equal(delta, startDelta))
    {
        delta *= 2;
        level++;
    }
    return level;
}

uint GridImp::getLevel() const
{
    return this->level;
}

void GridImp::setTriangularMeshDiscretizationStrategy(TriangularMeshDiscretizationStrategy* triangularMeshDiscretizationStrategy)
{
    this->triangularMeshDiscretizationStrategy = triangularMeshDiscretizationStrategy;
}

TriangularMeshDiscretizationStrategy * GridImp::getTriangularMeshDiscretizationStrategy()
{
    return this->triangularMeshDiscretizationStrategy;
}

uint GridImp::getNumberOfSolidBoundaryNodes() const
{
	return this->numberOfSolidBoundaryNodes;
}

void GridImp::setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes)
{
	if (numberOfSolidBoundaryNodes < INVALID_INDEX)
		this->numberOfSolidBoundaryNodes = numberOfSolidBoundaryNodes;
}

real GridImp::getQValue(const uint index, const uint dir) const
{
	const int qIndex = dir * this->numberOfSolidBoundaryNodes + this->qIndices[index];

	return this->qValues[qIndex];
}

uint GridImp::getQPatch(const uint index) const
{
    return this->qPatches[ this->qIndices[index] ];
}

void GridImp::setInnerRegionFromFinerGrid(bool innerRegionFromFinerGrid)
{
   this->innerRegionFromFinerGrid = innerRegionFromFinerGrid;
}

void GridImp::setNumberOfLayers(uint numberOfLayers)
{
    this->numberOfLayers = numberOfLayers;
}

// --------------------------------------------------------- //
//                  Set Sparse Indices                       //
// --------------------------------------------------------- //

void GridImp::findSparseIndices(SPtr<Grid> finerGrid)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Find sparse indices...";
    auto fineGrid = std::static_pointer_cast<GridImp>(finerGrid);
    
    this->updateSparseIndices();

#pragma omp parallel for
    for (int index = 0; index < (int)this->getSize(); index++)
        this->setNeighborIndices(index);

    if (fineGrid) {
        fineGrid->updateSparseIndices();
        this->findForGridInterfaceNewIndices(fineGrid);
    }

    const uint newGridSize = this->getSparseSize();
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "... done. new size: " << newGridSize
                  << ", delete nodes:" << this->getSize() - newGridSize << "\n";
}

void GridImp::findForGridInterfaceNewIndices(SPtr<GridImp> fineGrid)
{
#pragma omp parallel for
    for (int index = 0; index < (int)this->getNumberOfNodesCF(); index++)
        this->gridInterface->findForGridInterfaceSparseIndexCF(this, fineGrid.get(), index);

#pragma omp parallel for
    for (int index = 0; index < (int)this->getNumberOfNodesFC(); index++)
        this->gridInterface->findForGridInterfaceSparseIndexFC(this, fineGrid.get(), index);
}

void GridImp::updateSparseIndices()
{
    int removedNodes = 0;
    int newIndex = 0;
    for (uint index = 0; index < size; index++)
    {
        if (this->field.isInvalidCoarseUnderFine(index) || this->field.isInvalidOutOfGrid(index) || this->field.isInvalidSolid(index))
        {
            sparseIndices[index] = -1;
            removedNodes++;
        }
        else
        {
            sparseIndices[index] = newIndex;
            newIndex++;
        }
    }
    sparseSize = size - removedNodes;
}

void GridImp::findFluidNodeIndices(bool splitDomain) 
{
    // find sparse index of all fluid nodes
    this->fluidNodeIndices.clear();
    for (uint index = 0; index < this->size; index++) {
        int sparseIndex = this->getSparseIndex(index);
        if (sparseIndex == -1)
            continue;
        if (this->field.isFluid(index))
            this->fluidNodeIndices.push_back((uint)sparseIndex+1); // + 1 for numbering shift between GridGenerator and VF_GPU
    }

    // If splitDomain: find fluidNodeIndicesBorder and remove all indices in fluidNodeIndicesBorder from fluidNodeIndices
    if (splitDomain) {
        findFluidNodeIndicesBorder();
        std::sort(this->fluidNodeIndices.begin(), this->fluidNodeIndices.end());
        auto iterator = std::set_difference(this->fluidNodeIndices.begin(), this->fluidNodeIndices.end(),
                            this->fluidNodeIndicesBorder.begin(), this->fluidNodeIndicesBorder.end(),
                            this->fluidNodeIndices.begin());
        this->fluidNodeIndices.resize(iterator - this->fluidNodeIndices.begin());
    }
}

void GridImp::findFluidNodeIndicesBorder() {
    this->fluidNodeIndicesBorder.clear();

    // resize fluidNodeIndicesBorder (for better performance in copy operation)
    size_t newSize = 0;
    for (CommunicationIndices& ci : this->communicationIndices)
        newSize += ci.sendIndices.size();    
    this->fluidNodeIndicesBorder.reserve(newSize);

    // copy all send indices to fluidNodeIndicesBorder
    for (CommunicationIndices& ci : this->communicationIndices)
        std::copy(ci.sendIndices.begin(), ci.sendIndices.end(), std::back_inserter(this->fluidNodeIndicesBorder));

    // remove duplicate elements
    std::sort(this->fluidNodeIndicesBorder.begin(), this->fluidNodeIndicesBorder.end());
    this->fluidNodeIndicesBorder.erase(
        std::unique(this->fluidNodeIndicesBorder.begin(), this->fluidNodeIndicesBorder.end()),
        this->fluidNodeIndicesBorder.end());

    // + 1 for numbering shift between GridGenerator and VF_GPU
    for (size_t i = 0; i < this->fluidNodeIndicesBorder.size(); i++)
        this->fluidNodeIndicesBorder[i] = this->getSparseIndex(this->fluidNodeIndicesBorder[i])+1;
}

void GridImp::setNeighborIndices(uint index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    this->neighborIndexX[index]        = -1;
    this->neighborIndexY[index]        = -1;
    this->neighborIndexZ[index]        = -1;
    this->neighborIndexNegative[index] = -1;

    if (this->field.isStopper(index) || this->field.is(index, STOPPER_OUT_OF_GRID_BOUNDARY))
    {
        this->setStopperNeighborCoords(index);
        return;
    }
     
    if (this->sparseIndices[index] == -1)
        return;


    real neighborXCoord, neighborYCoord, neighborZCoord;
    this->getNeighborCoords(neighborXCoord, neighborYCoord, neighborZCoord, x, y, z);
    const int neighborX = getSparseIndex(neighborXCoord, y, z);
    const int neighborY = getSparseIndex(x, neighborYCoord, z);
    const int neighborZ = getSparseIndex(x, y, neighborZCoord);

    this->getNegativeNeighborCoords(neighborXCoord, neighborYCoord, neighborZCoord, x, y, z);
    const int neighborNegative = getSparseIndex(neighborXCoord, neighborYCoord, neighborZCoord);

    this->neighborIndexX[index]        = neighborX;
    this->neighborIndexY[index]        = neighborY;
    this->neighborIndexZ[index]        = neighborZ;
    this->neighborIndexNegative[index] = neighborNegative;
}

void GridImp::setStopperNeighborCoords(uint index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    if (vf::Math::lessEqual(x + delta, endX + (0.5 * delta)) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x + delta, y, z)))
        neighborIndexX[index] = getSparseIndex(x + delta, y, z);

    if (vf::Math::lessEqual(y + delta, endY + (0.5 * delta)) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x, y + delta, z)))
        neighborIndexY[index] = getSparseIndex(x, y + delta, z);

    if (vf::Math::lessEqual(z + delta, endZ + (0.5 * delta)) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x, y, z + delta)))
        neighborIndexZ[index] = getSparseIndex(x, y, z + delta);

    if (vf::Math::greaterEqual(x - delta, endX) && 
        vf::Math::greaterEqual(y - delta, endY) && 
        vf::Math::greaterEqual(z - delta, endZ) && 
        !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x - delta, y - delta, z - delta)))
    {
        neighborIndexNegative[index] = getSparseIndex(x - delta, y - delta, z - delta);
    }
}

void GridImp::getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const
{
    real coords[3] = { x, y, z };
    neighborX = getNeighborCoord(periodicityX, startX, coords, 0);
    neighborY = getNeighborCoord(periodicityY, startY, coords, 1);
    neighborZ = getNeighborCoord(periodicityZ, startZ, coords, 2);
}

real GridImp::getNeighborCoord(bool periodicity, real startCoord, real coords[3], int direction) const
{
    if (periodicity)
    {
        real neighborCoords[3] = {coords[0], coords[1] , coords[2] };
        neighborCoords[direction] = neighborCoords[direction] + delta;
        const int neighborIndex = this->transCoordToIndex(neighborCoords[0], neighborCoords[1], neighborCoords[2]);

        //////////////////////////////////////////////////////////////////////////

        if( field.is(neighborIndex, STOPPER_OUT_OF_GRID_BOUNDARY) )
            return getFirstFluidNode(coords, direction, startCoord);
        else
            return coords[direction] + delta;

    }
    
    return coords[direction] + delta;
}

void GridImp::getNegativeNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const
{
    real coords[3] = { x, y, z };

    neighborX = getNegativeNeighborCoord(periodicityX, endX, coords, 0);
    neighborY = getNegativeNeighborCoord(periodicityY, endY, coords, 1);
    neighborZ = getNegativeNeighborCoord(periodicityZ, endZ, coords, 2);
}

real GridImp::getNegativeNeighborCoord(bool periodicity, real startCoord, real coords[3], int direction) const
{
    if (periodicity)
    {
        real neighborCoords[3] = {coords[0], coords[1] , coords[2] };
        neighborCoords[direction] = neighborCoords[direction] - delta;
        const uint neighborIndex = this->transCoordToIndex(neighborCoords[0], neighborCoords[1], neighborCoords[2]);

        if(neighborIndex != INVALID_INDEX && !field.isStopperOutOfGrid(neighborIndex) && !field.is(neighborIndex, STOPPER_OUT_OF_GRID_BOUNDARY) )
            return coords[direction] - delta;

        return getLastFluidNode(coords, direction, startCoord);
    }
    
    return coords[direction] - delta;
}


real GridImp::getLastFluidNode(real coords[3], int direction, real startCoord) const
{
    coords[direction] = startCoord;
    uint index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    while (index != INVALID_INDEX && !field.isFluid(index))
    {
        coords[direction] -= delta;
        index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    }
    return coords[direction];
}

real GridImp::getFirstFluidNode(real coords[3], int direction, real startCoord) const
{
    coords[direction] = startCoord;
    uint index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    while (!field.isFluid(index))
    {
        coords[direction] += delta;
        index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    }
    return coords[direction];
}


int GridImp::getSparseIndex(const real &x, const real &y, const real &z) const
{
    const int matrixIndex = transCoordToIndex(x, y, z);
    return sparseIndices[matrixIndex];
}

// --------------------------------------------------------- //
//                    Find Interface                         //
// --------------------------------------------------------- //
void GridImp::findGridInterface(SPtr<Grid> finerGrid, LbmOrGks lbmOrGks)
{
    auto fineGrid          = std::static_pointer_cast<GridImp>(finerGrid);
    const auto coarseLevel = this->getLevel();
    const auto fineLevel   = fineGrid->getLevel();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "find interface level " << coarseLevel << " -> "
                  << fineLevel;

    this->gridInterface = new GridInterface();
    // TODO: this is stupid! concave refinements can easily have many more interface cells
    const uint sizeCF = 10 * (fineGrid->nx * fineGrid->ny + fineGrid->ny * fineGrid->nz + fineGrid->nx * fineGrid->nz);
    this->gridInterface->cf.coarse = new uint[sizeCF];
    this->gridInterface->cf.fine   = new uint[sizeCF];
    this->gridInterface->cf.offset = new uint[sizeCF];
    this->gridInterface->fc.coarse = new uint[sizeCF];
    this->gridInterface->fc.fine   = new uint[sizeCF];
    this->gridInterface->fc.offset = new uint[sizeCF];

    for (uint index = 0; index < this->getSize(); index++)
        this->findGridInterfaceCF(index, *fineGrid, lbmOrGks);

    for (uint index = 0; index < this->getSize(); index++)
        this->findGridInterfaceFC(index, *fineGrid);

    for (uint index = 0; index < this->getSize(); index++)
        this->findOverlapStopper(index, *fineGrid);

    if (lbmOrGks == GKS) {
        for (uint index = 0; index < this->getSize(); index++)
            this->findInvalidBoundaryNodes(index);
    }

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "  ... done. \n";
}

void GridImp::repairGridInterfaceOnMultiGPU(SPtr<Grid> fineGrid)
{
    this->gridInterface->repairGridInterfaceOnMultiGPU( shared_from_this(), std::static_pointer_cast<GridImp>(fineGrid) );
}

void GridImp::limitToSubDomain(SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks)
{
    for( uint index = 0; index < this->size; index++ ){

        real x, y, z;
        this->transIndexToCoords( index, x, y, z );

        {
            BoundingBox tmpSubDomainBox = *subDomainBox;

            // one layer for receive nodes and one for stoppers
            if( lbmOrGks == LBM )
                tmpSubDomainBox.extend(this->delta);
            
            if (!tmpSubDomainBox.isInside(x, y, z) 
                && ( this->getFieldEntry(index) == FLUID ||
                     this->getFieldEntry(index) == FLUID_CFC ||
                     this->getFieldEntry(index) == FLUID_CFF ||
                     this->getFieldEntry(index) == FLUID_FCC ||
                     this->getFieldEntry(index) == FLUID_FCF ||
                     this->getFieldEntry(index) == BC_SOLID ) )
            {   
                this->setFieldEntry(index, STOPPER_OUT_OF_GRID_BOUNDARY);
            }
        }

        {
            BoundingBox tmpSubDomainBox = *subDomainBox;

            // one layer for receive nodes and one for stoppers
            if( lbmOrGks == LBM )
                tmpSubDomainBox.extend(2.0 * this->delta);
            else
                tmpSubDomainBox.extend(1.0 * this->delta);

            if (!tmpSubDomainBox.isInside(x, y, z))
                this->setFieldEntry(index, INVALID_OUT_OF_GRID);
        }
    }
}

void GridImp::findGridInterfaceCF(uint index, GridImp& finerGrid, LbmOrGks lbmOrGks)
{
	if (lbmOrGks == LBM)
	{
		gridInterface->findInterfaceCF            (index, this, &finerGrid);
		gridInterface->findBoundaryGridInterfaceCF(index, this, &finerGrid);
	}
	else if (lbmOrGks == GKS)
		gridInterface->findInterfaceCF_GKS(index, this, &finerGrid);
}

void GridImp::findGridInterfaceFC(uint index, GridImp& finerGrid)
{
    gridInterface->findInterfaceFC(index, this, &finerGrid);
}

void GridImp::findOverlapStopper(uint index, GridImp& finerGrid)
{
    gridInterface->findOverlapStopper(index, this, &finerGrid);
}

void GridImp::findInvalidBoundaryNodes(uint index)
{
    gridInterface->findInvalidBoundaryNodes(index, this);
}

// --------------------------------------------------------- //
//                    Mesh Triangle                          //
// --------------------------------------------------------- //
void GridImp::mesh(Object* object)
{
    TriangularMesh *triangularMesh = dynamic_cast<TriangularMesh *>(object);
    if (triangularMesh)
        triangularMeshDiscretizationStrategy->discretize(triangularMesh, this, INVALID_SOLID, FLUID);
    else
		//new method for geometric primitives (not cell based) to be implemented
        this->discretize(object, INVALID_SOLID, FLUID);

    this->closeNeedleCells();

	#pragma omp parallel for
    for (int index = 0; index < (int)this->size; index++)
        this->findSolidStopperNode(index);

	//#pragma omp parallel for
    for (int index = 0; index < (int)this->size; index++) {
        this->findBoundarySolidNode(index);
    }
}


void GridImp::mesh(TriangularMesh &triangularMesh)
{
    const clock_t begin = clock();

#pragma omp parallel for
    for (int i = 0; i < triangularMesh.size; i++)
        this->mesh(triangularMesh.triangles[i]);

    const clock_t end = clock();
    const real time = (real)(real(end - begin) / CLOCKS_PER_SEC);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "time grid generation: " << time << "s\n";
}

void GridImp::mesh(Triangle &triangle)
{
    auto box = this->getBoundingBoxOnNodes(triangle);
    triangle.initalLayerThickness(getDelta());

    for (real x = box.minX; x <= box.maxX; x += delta)
    {
        for (real y = box.minY; y <= box.maxY; y += delta)
        {
            for (real z = box.minZ; z <= box.maxZ; z += delta)
            {
                const uint index = this->transCoordToIndex(x, y, z);
                if (!field.isFluid(index))
                    continue;

                const Vertex point(x, y, z);
                const char value = triangle.isUnderFace(point);
                //setDebugPoint(index, value);

                if (value == Q_DEPRECATED)
                    calculateQs(point, triangle);
            }
        }
    }
}

void GridImp::closeNeedleCells()
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start closeNeedleCells()\n";

    uint numberOfClosedNeedleCells = 0;

    do{
        numberOfClosedNeedleCells = 0;
#pragma omp parallel for reduction(+ : numberOfClosedNeedleCells)
        for (int index = 0; index < (int)this->size; index++) {
            if (this->closeCellIfNeedle(index))
                numberOfClosedNeedleCells++;
        }

        *logging::out << logging::Logger::INFO_INTERMEDIATE << numberOfClosedNeedleCells << " cells closed!\n";
    }
    while( numberOfClosedNeedleCells > 0 );
}

bool GridImp::closeCellIfNeedle(uint index)
{
    if( !this->getField().is( index, FLUID ) ) return false;

    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    bool noValidNeighborInX = this->getField().is( this->transCoordToIndex( x + this->delta, y,               z               ) , INVALID_SOLID ) &&
                              this->getField().is( this->transCoordToIndex( x - this->delta, y,               z               ) , INVALID_SOLID );
    bool noValidNeighborInY = this->getField().is( this->transCoordToIndex( x,               y + this->delta, z               ) , INVALID_SOLID ) &&
                              this->getField().is( this->transCoordToIndex( x,               y - this->delta, z               ) , INVALID_SOLID );
    bool noValidNeighborInZ = this->getField().is( this->transCoordToIndex( x,               y,               z + this->delta ) , INVALID_SOLID ) &&
                              this->getField().is( this->transCoordToIndex( x,               y,               z - this->delta ) , INVALID_SOLID );

    if( noValidNeighborInX || noValidNeighborInY || noValidNeighborInZ ){
        this->setFieldEntry(index, INVALID_SOLID);
        return true;
    }

    return false;
}

void GridImp::closeNeedleCellsThinWall()
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start closeNeedleCellsThinWall()\n";

    uint numberOfClosedNeedleCells = 0;

    do{
        numberOfClosedNeedleCells = 0;
#pragma omp parallel for reduction(+ : numberOfClosedNeedleCells)
        for (int index = 0; index < (int)this->size; index++) {
            if (this->closeCellIfNeedleThinWall(index))
                numberOfClosedNeedleCells++;
        }

        *logging::out << logging::Logger::INFO_INTERMEDIATE << numberOfClosedNeedleCells << " cells closed!\n";
    }
    while( numberOfClosedNeedleCells > 0 );
}

bool GridImp::closeCellIfNeedleThinWall(uint index)
{
    if( !this->getField().is( index, BC_SOLID ) ) return false;

    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    if( !this->hasNeighborOfType(index, FLUID) ){
        this->setFieldEntry(index, STOPPER_SOLID);
        return true;
    }

    return false;
}



void GridImp::findQs(Object* object) //TODO: enable qs for primitive objects
{
    TriangularMesh* triangularMesh = dynamic_cast<TriangularMesh*>(object);
    if (triangularMesh)
        findQs(*triangularMesh);
    else
        findQsPrimitive(object);
}

void GridImp::allocateQs() 
{
    this->qPatches = new uint[this->getNumberOfSolidBoundaryNodes()];

    for (uint i = 0; i < this->getNumberOfSolidBoundaryNodes(); i++)
        this->qPatches[i] = INVALID_INDEX;

    const uint numberOfQs = this->getNumberOfSolidBoundaryNodes() * (this->distribution.dir_end + 1);
    this->qValues         = new real[numberOfQs];
#pragma omp parallel for
    for (int i = 0; i < (int)numberOfQs; i++)
        this->qValues[i] = -1.0;
}

void GridImp::findQs(TriangularMesh &triangularMesh)
{
    const clock_t begin = clock();

    if( this->qComputationStage == qComputationStageType::ComputeQs )
        allocateQs();
    
    
#pragma omp parallel for
    for (int i = 0; i < triangularMesh.size; i++)
        this->findQs(triangularMesh.triangles[i]);

    const clock_t end = clock();
    const real time = (real)(real(end - begin) / CLOCKS_PER_SEC);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "time finding qs: " << time << "s\n";
}

void GridImp::findQs(Triangle &triangle)
{
    auto box = this->getBoundingBoxOnNodes(triangle);
    triangle.initalLayerThickness(getDelta());

    for (real x = box.minX; x <= box.maxX; x += delta)
    {
        for (real y = box.minY; y <= box.maxY; y += delta)
        {
            for (real z = box.minZ; z <= box.maxZ; z += delta)
            {
                const uint index = this->transCoordToIndex(x, y, z);
                //if (!field.isFluid(index))
                //    continue;

				if( index == INVALID_INDEX ) continue;

                const Vertex point(x, y, z);

                if( this->qComputationStage == qComputationStageType::ComputeQs ){
                    if(this->field.is(index, BC_SOLID))
                    {
					    calculateQs(index, point, triangle);
				    }
                }
                else if( this->qComputationStage == qComputationStageType::FindSolidBoundaryNodes )
                {
                    //if( this->field.is(index, BC_SOLID) || this->field.is(index, STOPPER_SOLID ) ) continue;

                    if( !this->field.is(index, FLUID) ) continue;

                    if( checkIfAtLeastOneValidQ(index, point, triangle) )
                    {
                        // similar as in void GridImp::findBoundarySolidNode(uint index)
                        this->field.setFieldEntry( index, BC_SOLID );
                        this->qIndices[index] = this->numberOfSolidBoundaryNodes++;
                    }
                }
            }
        }
    }
}

void GridImp::findQsPrimitive(Object * object)
{

    if( this->qComputationStage == qComputationStageType::ComputeQs )
        allocateQs();


    for( int index = 0; index < (int)this->size; index++ )
    {

        if( this->qIndices[index] == INVALID_INDEX ) continue;

        real x,y,z;

        this->transIndexToCoords(index,x,y,z);
        
        const Vertex point(x, y, z);

        if( this->qComputationStage == qComputationStageType::ComputeQs ){
            if(this->field.is(index, BC_SOLID))
            {
				calculateQs(index, point, object);
			}
        }
        else if( this->qComputationStage == qComputationStageType::FindSolidBoundaryNodes )
        {
            if( !this->field.is(index, FLUID) ) continue;

            if( checkIfAtLeastOneValidQ(index, point, object) )
            {
                // similar as in void GridImp::findBoundarySolidNode(uint index)
                this->field.setFieldEntry( index, BC_SOLID );
                this->qIndices[index] = this->numberOfSolidBoundaryNodes++;
            }
        }

    }
}

void GridImp::calculateQs(const uint index, const Vertex &point, Object* object) const
{
    Vertex pointOnTriangle, direction;

	real subdistance;
	int error;
	for (int i = distribution.dir_start; i <= distribution.dir_end; i++)
	{
		direction = Vertex( real(distribution.dirs[i * DIMENSION + 0]), 
                            real(distribution.dirs[i * DIMENSION + 1]),
			                real(distribution.dirs[i * DIMENSION + 2]) );

		uint neighborIndex = this->transCoordToIndex(point.x + direction.x * this->delta,
													    point.y + direction.y * this->delta,
													    point.z + direction.z * this->delta);

		if (neighborIndex == INVALID_INDEX) continue;

		error = object->getIntersection(point, direction, pointOnTriangle, subdistance);

		subdistance /= this->delta;

		if (error == 0 && vf::Math::lessEqual(subdistance, 1.0) && vf::Math::greaterEqual(subdistance, 0.0))
		{
			if ( -0.5        > this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] ||
                    subdistance < this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] )
			{

				this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] = subdistance;
                    
                this->qPatches[ this->qIndices[index] ] = 0;

			}
		}
	}
}

bool GridImp::checkIfAtLeastOneValidQ(const uint index, const Vertex &point, Object* object) const
{
    Vertex pointOnTriangle, direction;

	real subdistance;
	int error;
	for (int i = distribution.dir_start; i <= distribution.dir_end; i++)
	{
		direction = Vertex( real(distribution.dirs[i * DIMENSION + 0]), 
                            real(distribution.dirs[i * DIMENSION + 1]),
			                real(distribution.dirs[i * DIMENSION + 2]) );

		uint neighborIndex = this->transCoordToIndex(point.x + direction.x * this->delta,
													 point.y + direction.y * this->delta,
													 point.z + direction.z * this->delta);

		if (neighborIndex == INVALID_INDEX) continue;

		error = object->getIntersection(point, direction, pointOnTriangle, subdistance);

		subdistance /= this->delta;

		if (error == 0 && vf::Math::lessEqual(subdistance, 1.0) && vf::Math::greaterEqual(subdistance, 0.0))
		{
			return true;
		}
	}
    return false;
}

void GridImp::setDebugPoint(uint index, int pointValue)
{
    if (field.isInvalidCoarseUnderFine(index) && pointValue == INVALID_SOLID)
        field.setFieldEntry(index, pointValue);

    if(!field.isInvalidSolid(index) && !field.isQ(index) && !field.isInvalidCoarseUnderFine(index) && pointValue != 3 && pointValue != 2)
        field.setFieldEntry(index, pointValue);
}

void GridImp::calculateQs(const Vertex &point, const Triangle &triangle) const   // NOT USED !!!!
{
    Vertex pointOnTriangle, direction;
    real subdistance;
    int error;
    for (int i = distribution.dir_start; i <= distribution.dir_end; i++)
    {
#if defined(__CUDA_ARCH__)
        direction = Vertex(DIRECTIONS[i][0], DIRECTIONS[i][1], DIRECTIONS[i][2]);
#else
        direction = Vertex(real(distribution.dirs[i * DIMENSION + 0]), real(distribution.dirs[i * DIMENSION + 1]),
            real(distribution.dirs[i * DIMENSION + 2]));
#endif

        error = triangle.getTriangleIntersection(point, direction, pointOnTriangle, subdistance);

		subdistance /= this->delta;

        if (error == 0 && subdistance < 1.0 && subdistance > 0.0)
        {
            distribution.f[i*size + transCoordToIndex(point.x, point.y, point.z)] = subdistance;
        }
    }
}


void GridImp::calculateQs(const uint index, const Vertex &point, const Triangle &triangle) const
{
	Vertex pointOnTriangle, direction;
	real subdistance;
	int error;
	for (int i = distribution.dir_start; i <= distribution.dir_end; i++)
	{
#if defined(__CUDA_ARCH__)
		direction = Vertex(DIRECTIONS[i][0], DIRECTIONS[i][1], DIRECTIONS[i][2]);
#else
		direction = Vertex( real(distribution.dirs[i * DIMENSION + 0]), 
                            real(distribution.dirs[i * DIMENSION + 1]),
			                real(distribution.dirs[i * DIMENSION + 2]) );
#endif

		uint neighborIndex = this->transCoordToIndex(point.x + direction.x * this->delta,
													 point.y + direction.y * this->delta,
													 point.z + direction.z * this->delta);

		if (neighborIndex == INVALID_INDEX) continue;

		error = triangle.getTriangleIntersection(point, direction, pointOnTriangle, subdistance);

		subdistance /= this->delta;

		if (error == 0 && vf::Math::lessEqual(subdistance, 1.0) && vf::Math::greaterEqual(subdistance, 0.0))
		{
			if ( -0.5        > this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] ||
                 subdistance < this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] )
			{
				this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] = subdistance;

                this->qPatches[ this->qIndices[index] ] = triangle.patchIndex;
			}
		}
	}
}

bool GridImp::checkIfAtLeastOneValidQ(const uint index, const Vertex & point, const Triangle & triangle) const
{
	Vertex pointOnTriangle, direction;
	real subdistance;
	int error;
	for (int i = distribution.dir_start; i <= distribution.dir_end; i++)
	{
#if defined(__CUDA_ARCH__)
		direction = Vertex(DIRECTIONS[i][0], DIRECTIONS[i][1], DIRECTIONS[i][2]);
#else
		direction = Vertex(real(distribution.dirs[i * DIMENSION + 0]), 
                           real(distribution.dirs[i * DIMENSION + 1]),
			               real(distribution.dirs[i * DIMENSION + 2]));
#endif

		uint neighborIndex = this->transCoordToIndex(point.x + direction.x * this->delta,
													 point.y + direction.y * this->delta,
													 point.z + direction.z * this->delta);
		if (neighborIndex == INVALID_INDEX) continue;

		error = triangle.getTriangleIntersection(point, direction, pointOnTriangle, subdistance);

		subdistance /= this->delta;

		if (error == 0 && vf::Math::lessEqual(subdistance, 1.0) && vf::Math::greaterEqual(subdistance, 0.0))
		{
			return true;
		}
	}
    return false;
}

void GridImp::findCommunicationIndices(int direction, SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks)
{
    for( uint index = 0; index < this->size; index++ ){
        
        real x, y, z;
        this->transIndexToCoords(index, x, y, z);

        if( this->getFieldEntry(index) == INVALID_OUT_OF_GRID ||
            this->getFieldEntry(index) == INVALID_SOLID ||
            this->getFieldEntry(index) == INVALID_COARSE_UNDER_FINE ||
            this->getFieldEntry(index) == STOPPER_OUT_OF_GRID ||
            this->getFieldEntry(index) == STOPPER_COARSE_UNDER_FINE ) continue;

        if( lbmOrGks == LBM && this->getFieldEntry(index) == STOPPER_OUT_OF_GRID_BOUNDARY ) continue;
        if( lbmOrGks == LBM && this->getFieldEntry(index) == STOPPER_SOLID ) continue;

        if( direction == CommunicationDirections::MX ) findCommunicationIndex( index, x, subDomainBox->minX, direction);
        if( direction == CommunicationDirections::PX ) findCommunicationIndex( index, x, subDomainBox->maxX, direction);
        if( direction == CommunicationDirections::MY ) findCommunicationIndex( index, y, subDomainBox->minY, direction);
        if( direction == CommunicationDirections::PY ) findCommunicationIndex( index, y, subDomainBox->maxY, direction);
        if( direction == CommunicationDirections::MZ ) findCommunicationIndex( index, z, subDomainBox->minZ, direction);
        if( direction == CommunicationDirections::PZ ) findCommunicationIndex( index, z, subDomainBox->maxZ, direction);
    }
    real xS, yS, zS;
        this->transIndexToCoords(this->communicationIndices[direction].sendIndices[0], xS, yS, zS);
    real xR, yR, zR;
        this->transIndexToCoords(this->communicationIndices[direction].receiveIndices[0], xR, yR, zR);
    std::cout << "Dir: " << direction << " Sender: " << xS << " " << yS<< " "  << zS<< " " 
                                      << this->communicationIndices[direction].sendIndices[0] << " "
                                      << (int)this->getFieldEntry(this->communicationIndices[direction].sendIndices[0]) 
                                      << " Receiver: "  << xR << " " << yR << " " << zR << " " 
                                      << this->communicationIndices[direction].receiveIndices[0] << " "
                                      << (int)this->getFieldEntry(this->communicationIndices[direction].receiveIndices[0])
                                      << std::endl<< std::endl;

}

void GridImp::findCommunicationIndex( uint index, real coordinate, real limit, int direction ){
    // negative direction get a negative sign
    real s = ( direction % 2 == 0 ) ? ( -1.0 ) : ( 1.0 );  
    bool send = false;
    bool rec = false;
	if (std::abs(coordinate - (limit + s * 0.5 * this->delta)) < 0.1 * this->delta) {
		this->communicationIndices[direction].receiveIndices.push_back(index);
        rec = true;
	}

	if (std::abs(coordinate - (limit - s * 0.5 * this->delta)) < 0.1 * this->delta) {
		this->communicationIndices[direction].sendIndices.push_back(index);
        send = true;
	}
    if( false && (send || rec ) ) 
    {
        std::cout << "Send Idx: " << index << " limit: " << limit << " coord: " << coordinate << " dir: " << direction << " send: " << send << " receive: " << rec << std::endl;
        std::cout << "IDX: " << (int)this->getFieldEntry(index)<< std::endl;
    }
}

bool GridImp::isSendNode(int index) const
{
    bool isSendNode = false;
    for (size_t direction = 0; direction < this->communicationIndices.size(); direction++)
        if (std::find(this->communicationIndices[direction].sendIndices.begin(),
                      this->communicationIndices[direction].sendIndices.end(), index) != this->communicationIndices[direction].sendIndices.end())
            isSendNode = true;
    return isSendNode;
}

bool GridImp::isReceiveNode(int index) const
{
    bool isReceiveNode = false;
    for (size_t direction = 0; direction < this->communicationIndices.size(); direction++)
        if (std::find(this->communicationIndices[direction].receiveIndices.begin(),
                      this->communicationIndices[direction].receiveIndices.end(),
                      index) != this->communicationIndices[direction].receiveIndices.end())
            isReceiveNode = true;
    return isReceiveNode;
}

uint GridImp::getNumberOfSendNodes(int direction)
{
    return (uint)this->communicationIndices[direction].sendIndices.size();
}

uint GridImp::getNumberOfReceiveNodes(int direction)
{
    return (uint)this->communicationIndices[direction].receiveIndices.size();
}

uint GridImp::getSendIndex(int direction, uint index)
{
    return this->communicationIndices[direction].sendIndices[ index ];
}

uint GridImp::getReceiveIndex(int direction, uint index)
{
    return this->communicationIndices[direction].receiveIndices[ index ];
}

void GridImp::repairCommunicationIndices(int direction)
{
    this->communicationIndices[direction].sendIndices.insert( this->communicationIndices[direction].sendIndices.end(), 
                                                              this->communicationIndices[direction+1].sendIndices.begin(), 
                                                              this->communicationIndices[direction+1].sendIndices.end() );



    this->communicationIndices[direction+1].receiveIndices.insert( this->communicationIndices[direction+1].receiveIndices.end(), 
                                                                 this->communicationIndices[direction].receiveIndices.begin(), 
                                                                 this->communicationIndices[direction].receiveIndices.end() );

    this->communicationIndices[direction].receiveIndices = this->communicationIndices[direction+1].receiveIndices;






    *logging::out << logging::Logger::INFO_INTERMEDIATE << "size send " << (int)this->communicationIndices[direction].sendIndices.size() << "\n";
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "recv send " << (int)this->communicationIndices[direction].receiveIndices.size() << "\n";
}


// --------------------------------------------------------- //
//                        Getter                             //
// --------------------------------------------------------- //
int GridImp::getSparseIndex(uint matrixIndex) const
{
    return this->sparseIndices[matrixIndex];
}

real* GridImp::getDistribution() const
{
    return this->distribution.f;
}

int* GridImp::getDirection() const
{
    return this->distribution.dirs;
}

int GridImp::getStartDirection() const
{
    return this->distribution.dir_start;
}

int GridImp::getEndDirection() const
{
    return this->distribution.dir_end;
}

BoundingBox GridImp::getBoundingBoxOnNodes(Triangle &triangle) const
{
    real minX, maxX, minY, maxY, minZ, maxZ;
    triangle.setMinMax(minX, maxX, minY, maxY, minZ, maxZ);

    int minXIndex = std::lround(floor((minX - this->startX) / this->delta)) - 1;
    int minYIndex = std::lround(floor((minY - this->startY) / this->delta)) - 1;
    int minZIndex = std::lround(floor((minZ - this->startZ) / this->delta)) - 1;

    int maxXIndex = std::lround(ceil((maxX - this->startX) / this->delta)) + 1;
    int maxYIndex = std::lround(ceil((maxY - this->startY) / this->delta)) + 1;
    int maxZIndex = std::lround(ceil((maxZ - this->startZ) / this->delta)) + 1;

    minX = this->startX + minXIndex * this->delta;
    minY = this->startY + minYIndex * this->delta;
    minZ = this->startZ + minZIndex * this->delta;

    maxX = this->startX + maxXIndex * this->delta;
    maxY = this->startY + maxYIndex * this->delta;
    maxZ = this->startZ + maxZIndex * this->delta;

    return BoundingBox(minX, maxX, minY, maxY, minZ, maxZ);
}

Vertex GridImp::getMinimumOnNode(Vertex exact) const // deprecated
{
    const real minX = getMinimumOnNodes(exact.x, vf::Math::getDecimalPart(startX), delta);
    const real minY = getMinimumOnNodes(exact.y, vf::Math::getDecimalPart(startY), delta);
    const real minZ = getMinimumOnNodes(exact.z, vf::Math::getDecimalPart(startZ), delta);
    return Vertex(minX, minY, minZ);
}

real GridImp::getMinimumOnNodes(const real &minExact, const real &decimalStart, const real &delta) // deprecated
{
    real minNode = ceil(minExact - 1.0);
    minNode += decimalStart;
    while (minNode > minExact)
        minNode -= delta;

    while (minNode + delta < minExact)
        minNode += delta;
    return minNode;
}

Vertex GridImp::getMaximumOnNode(Vertex exact) const // deprecated
{
    const real maxX = getMaximumOnNodes(exact.x, vf::Math::getDecimalPart(startX), delta);
    const real maxY = getMaximumOnNodes(exact.y, vf::Math::getDecimalPart(startY), delta);
    const real maxZ = getMaximumOnNodes(exact.z, vf::Math::getDecimalPart(startZ), delta);
    return Vertex(maxX, maxY, maxZ);
}

real GridImp::getMaximumOnNodes(const real &maxExact, const real &decimalStart, const real &delta) // deprecated
{
    real maxNode = ceil(maxExact - 1.0);
    maxNode += decimalStart;

    while (maxNode <= maxExact)
        maxNode += delta;
    return maxNode;
}

uint GridImp::getXIndex(real x) const 
{ 
    return std::lround((x - startX) / delta); 
}

uint GridImp::getYIndex(real y) const
{ 
    return std::lround((y - startY) / delta); 
}

uint GridImp::getZIndex(real z) const
{ 
    return std::lround((z - startZ) / delta); 
}

real GridImp::getDelta() const
{
    return delta;
}

uint GridImp::getSize() const
{
    return this->size;
}

uint GridImp::getSparseSize() const
{
    return this->sparseSize; 
}

uint GridImp::getNumberOfFluidNodes() const { 
    return (uint)this->fluidNodeIndices.size(); 
}

Field GridImp::getField() const
{
    return this->field;
}

char GridImp::getFieldEntry(uint index) const
{
    return this->field.getFieldEntry(index);
}

void GridImp::setFieldEntry(uint matrixIndex, char type)
{
    this->field.setFieldEntry(matrixIndex, type);
}


real GridImp::getStartX() const
{
    return startX;
}

real GridImp::getStartY() const
{
    return startY;
}

real GridImp::getStartZ() const
{
    return startZ;
}

real GridImp::getEndX() const
{
    return endX;
}

real GridImp::getEndY() const
{
    return endY;
}

real GridImp::getEndZ() const
{
    return endZ;
}

uint GridImp::getNumberOfNodesX() const
{
    return nx;
}

uint GridImp::getNumberOfNodesY() const
{
    return ny;
}

uint GridImp::getNumberOfNodesZ() const
{
    return nz;
}


int* GridImp::getNeighborsX() const
{
    return this->neighborIndexX;
}

int* GridImp::getNeighborsY() const
{
    return this->neighborIndexY;
}

int* GridImp::getNeighborsZ() const
{
    return this->neighborIndexZ;
}

int* GridImp::getNeighborsNegative() const
{
    return this->neighborIndexNegative;
}

uint GridImp::getNumberOfNodesCF() const
{
    if(this->gridInterface)
        return this->gridInterface->cf.numberOfEntries;
    return 0;
}

uint GridImp::getNumberOfNodesFC() const
{
    if (this->gridInterface)
        return this->gridInterface->fc.numberOfEntries;
    return 0;
}

uint* GridImp::getCF_coarse() const
{
    return this->gridInterface->cf.coarse;
}

uint* GridImp::getCF_fine() const
{
    return this->gridInterface->cf.fine;
}

uint * GridImp::getCF_offset() const
{
    return this->gridInterface->cf.offset;
}

uint* GridImp::getFC_coarse() const
{
    return this->gridInterface->fc.coarse;
}

uint* GridImp::getFC_fine() const
{
    return this->gridInterface->fc.fine;
}

uint * GridImp::getFC_offset() const
{
    return this->gridInterface->fc.offset;
}

void GridImp::getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const
{
    getGridInterface(iCellCfc, this->gridInterface->cf.coarse, this->gridInterface->cf.numberOfEntries);
    getGridInterface(iCellCff, this->gridInterface->cf.fine, this->gridInterface->cf.numberOfEntries);
    getGridInterface(iCellFcc, this->gridInterface->fc.coarse, this->gridInterface->fc.numberOfEntries);
    getGridInterface(iCellFcf, this->gridInterface->fc.fine, this->gridInterface->fc.numberOfEntries);
}

void GridImp::getGridInterface(uint* gridInterfaceList, const uint* oldGridInterfaceList, uint size)
{
    for (uint i = 0; i < size; i++)
        gridInterfaceList[i] = oldGridInterfaceList[i] + 1; // + 1 for numbering shift between GridGenerator and VF_GPU
}

bool GridImp::isSparseIndexInFluidNodeIndicesBorder(uint &sparseIndex) const
{
    return std::find(this->fluidNodeIndicesBorder.begin(), this->fluidNodeIndicesBorder.end(), sparseIndex) !=
           this->fluidNodeIndicesBorder.end();
}

#define GEOFLUID 19
#define GEOSOLID 16

void GridImp::getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, uint *geo) const
{
    xCoords[0] = 0;
    yCoords[0] = 0;
    zCoords[0] = 0;
    neighborX[0] = 0;
    neighborY[0] = 0;
    neighborZ[0] = 0;
    geo[0] = GEOSOLID;

    int nodeNumber = 0;
    for (uint i = 0; i < this->size; i++)
    {
        if (this->sparseIndices[i] == -1)
            continue;

        real x, y, z;
        this->transIndexToCoords(i, x, y, z);

        // + 1 for numbering shift between GridGenerator and VF_GPU
        const uint neighborXIndex        = uint(this->neighborIndexX[i] + 1);
        const uint neighborYIndex        = uint(this->neighborIndexY[i] + 1);
        const uint neighborZIndex        = uint(this->neighborIndexZ[i] + 1);
        const uint neighborNegativeIndex = uint(this->neighborIndexNegative[i] + 1);

        const uint type = uint(this->field.isFluid(i) ? GEOFLUID : GEOSOLID);

        xCoords[nodeNumber + 1] = x;
        yCoords[nodeNumber + 1] = y;
        zCoords[nodeNumber + 1] = z;

        neighborX       [nodeNumber + 1] = neighborXIndex;
        neighborY       [nodeNumber + 1] = neighborYIndex;
        neighborZ       [nodeNumber + 1] = neighborZIndex;
        neighborNegative[nodeNumber + 1] = neighborNegativeIndex;

        geo[nodeNumber + 1] = type;
        nodeNumber++;
    }
}

void GridImp::getFluidNodeIndices(uint *fluidNodeIndices) const 
{ 
    for (uint nodeNumber = 0; nodeNumber < (uint)this->fluidNodeIndices.size(); nodeNumber++)
        fluidNodeIndices[nodeNumber] = this->fluidNodeIndices[nodeNumber];
}

uint GridImp::getNumberOfFluidNodesBorder() const 
{ 
    return (uint)this->fluidNodeIndicesBorder.size(); 
}

void GridImp::getFluidNodeIndicesBorder(uint *fluidNodeIndicesBorder) const 
{
    for (uint nodeNumber = 0; nodeNumber < (uint)this->fluidNodeIndicesBorder.size(); nodeNumber++)
        fluidNodeIndicesBorder[nodeNumber] = this->fluidNodeIndicesBorder[nodeNumber];
}

void GridImp::addFluidNodeIndicesMacroVars(std::vector<uint> _fluidNodeIndicesMacroVars)
{
    size_t newSize = this->fluidNodeIndicesMacroVars.size()+_fluidNodeIndicesMacroVars.size();
    this->fluidNodeIndicesMacroVars.reserve(newSize);
    std::copy(_fluidNodeIndicesMacroVars.begin(), _fluidNodeIndicesMacroVars.end(), std::back_inserter(this->fluidNodeIndicesMacroVars));
}

void GridImp::addFluidNodeIndicesApplyBodyForce(std::vector<uint> _fluidNodeIndicesApplyBodyForce)
{    
    
    size_t newSize = this->fluidNodeIndicesApplyBodyForce.size()+_fluidNodeIndicesApplyBodyForce.size();
    this->fluidNodeIndicesApplyBodyForce.reserve(newSize);
    std::copy(_fluidNodeIndicesApplyBodyForce.begin(), _fluidNodeIndicesApplyBodyForce.end(), std::back_inserter(this->fluidNodeIndicesApplyBodyForce));
}

void GridImp::addFluidNodeIndicesAllFeatures(std::vector<uint> _fluidNodeIndicesAllFeatures) 
{

    size_t newSize = this->fluidNodeIndicesAllFeatures.size()+_fluidNodeIndicesAllFeatures.size();
    this->fluidNodeIndicesAllFeatures.reserve(newSize);
    std::copy(_fluidNodeIndicesAllFeatures.begin(), _fluidNodeIndicesAllFeatures.end(), std::back_inserter(this->fluidNodeIndicesAllFeatures));
}

void GridImp::sortFluidNodeIndicesMacroVars()
{
    if(this->fluidNodeIndicesMacroVars.size()>0)
    {
        sort(this->fluidNodeIndicesMacroVars.begin(), this->fluidNodeIndicesMacroVars.end());
        // Remove duplicates
        this->fluidNodeIndicesMacroVars.erase( unique( this->fluidNodeIndicesMacroVars.begin(), this->fluidNodeIndicesMacroVars.end() ), this->fluidNodeIndicesMacroVars.end() );

         // Remove indices of fluidNodeIndicesAllFeatures from fluidNodeIndicesMacroVars
        if(this->fluidNodeIndicesAllFeatures.size()>0)
        {
            this->fluidNodeIndicesMacroVars.erase(   std::remove_if(   this->fluidNodeIndicesMacroVars.begin(), this->fluidNodeIndicesMacroVars.end(), 
                                                        [&](auto x){return binary_search(fluidNodeIndicesAllFeatures.begin(),fluidNodeIndicesAllFeatures.end(),x);} ),
                                            this->fluidNodeIndicesMacroVars.end()
                                        );
        }

        // Remove indices of fluidNodeIndicesMacroVars from fluidNodeIndices
        this->fluidNodeIndices.erase(   std::remove_if(   this->fluidNodeIndices.begin(), this->fluidNodeIndices.end(), 
                                                        [&](auto x){return binary_search(fluidNodeIndicesMacroVars.begin(),fluidNodeIndicesMacroVars.end(),x);} ),
                                        this->fluidNodeIndices.end()
                                    );
    }
}

void GridImp::sortFluidNodeIndicesApplyBodyForce()
{
    if(this->fluidNodeIndicesApplyBodyForce.size()>0)
    {
        sort(this->fluidNodeIndicesApplyBodyForce.begin(), this->fluidNodeIndicesApplyBodyForce.end());
        // Remove duplicates
        this->fluidNodeIndicesApplyBodyForce.erase( unique( this->fluidNodeIndicesApplyBodyForce.begin(), this->fluidNodeIndicesApplyBodyForce.end() ), this->fluidNodeIndicesApplyBodyForce.end() );

         // Remove indices of fluidNodeIndicesAllFeatures from fluidNodeIndicesMacroVars
        if(this->fluidNodeIndicesAllFeatures.size()>0)
        {
            this->fluidNodeIndicesApplyBodyForce.erase(   std::remove_if(   this->fluidNodeIndicesApplyBodyForce.begin(), this->fluidNodeIndicesApplyBodyForce.end(), 
                                                        [&](auto x){return binary_search(fluidNodeIndicesAllFeatures.begin(),fluidNodeIndicesAllFeatures.end(),x);} ),
                                            this->fluidNodeIndicesApplyBodyForce.end()
                                        );
        }

        // Remove indices of fluidNodeIndicesMacroVars from fluidNodeIndices
        this->fluidNodeIndices.erase(   std::remove_if(   this->fluidNodeIndices.begin(), this->fluidNodeIndices.end(), 
                                                        [&](auto x){return binary_search(fluidNodeIndicesApplyBodyForce.begin(),fluidNodeIndicesApplyBodyForce.end(),x);} ),
                                        this->fluidNodeIndices.end()
                                    );
    }
}

void GridImp::sortFluidNodeIndicesAllFeatures()
{
    if(this->fluidNodeIndicesAllFeatures.size()>0)
    {
        sort(this->fluidNodeIndicesAllFeatures.begin(), this->fluidNodeIndicesAllFeatures.end());
        // Remove duplicates
        this->fluidNodeIndicesAllFeatures.erase( unique( this->fluidNodeIndicesAllFeatures.begin(), this->fluidNodeIndicesAllFeatures.end() ), this->fluidNodeIndicesAllFeatures.end() );
        // Remove indices of fluidNodeIndicesMacroVars from fluidNodeIndices
        this->fluidNodeIndices.erase(   std::remove_if(   this->fluidNodeIndices.begin(), this->fluidNodeIndices.end(), 
                                                        [&](auto x){return binary_search(fluidNodeIndicesAllFeatures.begin(),fluidNodeIndicesAllFeatures.end(),x);} ),
                                        this->fluidNodeIndices.end()
                                    );
    }
}

uint GridImp::getNumberOfFluidNodeIndicesMacroVars() const { 
    return (uint)this->fluidNodeIndicesMacroVars.size(); 
}

uint GridImp::getNumberOfFluidNodeIndicesApplyBodyForce() const { 
    return (uint)this->fluidNodeIndicesApplyBodyForce.size(); 
}

uint GridImp::getNumberOfFluidNodeIndicesAllFeatures() const { 
    return (uint)this->fluidNodeIndicesAllFeatures.size(); 
}

void GridImp::getFluidNodeIndicesMacroVars(uint *_fluidNodeIndicesMacroVars) const 
{
    std::copy(fluidNodeIndicesMacroVars.begin(), fluidNodeIndicesMacroVars.end(), _fluidNodeIndicesMacroVars);       
}
void GridImp::getFluidNodeIndicesApplyBodyForce(uint *_fluidNodeIndicesApplyBodyForce) const 
{
    std::copy(fluidNodeIndicesApplyBodyForce.begin(), fluidNodeIndicesApplyBodyForce.end(), _fluidNodeIndicesApplyBodyForce);
}
void GridImp::getFluidNodeIndicesAllFeatures(uint *_fluidNodeIndicesAllFeatures) const 
{
    std::copy(fluidNodeIndicesAllFeatures.begin(), fluidNodeIndicesAllFeatures.end(), _fluidNodeIndicesAllFeatures);
}






void GridImp::print() const
{
    printf("min: (%2.4f, %2.4f, %2.4f), max: (%2.4f, %2.4f, %2.4f), size: %d, delta: %2.4f\n", startX, startY, startZ,
           endX, endY, endZ, size, delta);
    if(this->gridInterface)
        this->gridInterface->print();
}
