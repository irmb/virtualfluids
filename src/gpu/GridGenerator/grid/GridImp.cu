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
//! \file GridImp.cu
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "GridImp.h"

#include <stdio.h>
#include <time.h>
#include <iostream>
#include <omp.h>
#include <sstream>

#include "global.h"

#include "geometries/Object.h"
#include "geometries/Vertex/Vertex.h"
#include "geometries/BoundingBox/BoundingBox.h"

#include "grid/GridStrategy/GridStrategy.h"
#include "grid/distributions/Distribution.h"
#include "grid/Field.h"
#include "grid/NodeValues.h"

#include "utilities/math/Math.h"

int DIRECTIONS[DIR_END_MAX][DIMENSION];


GridImp::GridImp(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution distribution, uint level) 
            : object(object), 
    startX(startX),
    startY(startY),
    startZ(startZ),
    endX(endX),
    endY(endY),
    endZ(endZ),
    delta(delta),
    gridStrategy(gridStrategy),
    distribution(distribution),
    level(level),
    periodicityX(false),
    periodicityY(false),
    periodicityZ(false),
    neighborIndexX(nullptr),
    neighborIndexY(nullptr),
    neighborIndexZ(nullptr),
    neighborIndexNegative(nullptr),
    sparseIndices(nullptr)
{
    initalNumberOfNodesAndSize();
}

SPtr<GridImp> GridImp::makeShared(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level)
{
    SPtr<GridImp> grid(new GridImp(object, startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, d, level));
    return grid;
}


void GridImp::initalNumberOfNodesAndSize()
{
    const real length = endX - startX;
    const real width = endY - startY;
    const real height = endZ - startZ;

    nx = lround((length + delta) / delta);
    ny = lround((width + delta) / delta);
    nz = lround((height + delta) / delta);

    this->size = nx * ny * nz;
    this->sparseSize = size;
    distribution.setSize(size);
}

void GridImp::inital(const SPtr<Grid> fineGrid, uint numberOfLayers)
{
    field = Field(gridStrategy, size);
    field.allocateMemory();
    gridStrategy->allocateGridMemory(shared_from_this());
    
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start initalNodesToOutOfGrid()\n";
    gridStrategy->initalNodesToOutOfGrid(shared_from_this());
    
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start findInnerNodes()\n";
    this->object->findInnerNodes( shared_from_this() );
    
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Start findEndOfGridStopperNodes()\n";
	gridStrategy->findEndOfGridStopperNodes(shared_from_this());

    *logging::out << logging::Logger::INFO_INTERMEDIATE
        << "Grid created: " << "from (" << this->startX << ", " << this->startY << ", " << this->startZ << ") to (" << this->endX << ", " << this->endY << ", " << this->endZ << ")\n"
        << "nodes: " << this->nx << " x " << this->ny << " x " << this->nz << " = " << this->size << "\n";
}

void GridImp::initalNodeToOutOfGrid(uint index)
{
    this->field.setFieldEntryToInvalidOutOfGrid(index);
}

void GridImp::freeMemory()
{
    gridStrategy->freeMemory(shared_from_this());
}

GridImp::GridImp()
{
}

GridImp::~GridImp()
{
}

void GridImp::findInnerNode(uint index)
{
    this->sparseIndices[index] = index;

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
    this->periodicityX = periodicityX;
}

void GridImp::setPeriodicityY(bool periodicity)
{
    this->periodicityY = periodicityY;
}

void GridImp::setPeriodicityZ(bool periodicity)
{
    this->periodicityZ = periodicityZ;
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

    x = index % nx;
    y = (index / nx) % ny;
    z = ((index / nx) / ny) % nz;

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

// --------------------------------------------------------- //
//                  Set Sparse Indices                       //
// --------------------------------------------------------- //

void GridImp::findSparseIndices(SPtr<Grid> fineGrid)
{
    this->gridStrategy->findSparseIndices(shared_from_this(), std::static_pointer_cast<GridImp>(fineGrid));
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

    if (vf::Math::lessEqual(x + delta, endX) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x + delta, y, z)))
        neighborIndexX[index] = getSparseIndex(x + delta, y, z);

    if (vf::Math::lessEqual(y + delta, endY) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x, y + delta, z)))
        neighborIndexY[index] = getSparseIndex(x, y + delta, z);

    if (vf::Math::lessEqual(z + delta, endZ) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x, y, z + delta)))
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
        const int neighborIndex = this->transCoordToIndex(neighborCoords[0], neighborCoords[1], neighborCoords[2]);

        if(neighborIndex != INVALID_INDEX && !field.isStopperOutOfGrid(neighborIndex) && !field.is(neighborIndex, STOPPER_OUT_OF_GRID_BOUNDARY) )
            return coords[direction] - delta;

        return getLastFluidNode(coords, direction, startCoord);
    }
    
    return coords[direction] - delta;
}


real GridImp::getLastFluidNode(real coords[3], int direction, real startCoord) const
{
    coords[direction] = startCoord;
    int index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
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

uint GridImp::getXIndex(real x) const
{
    return lround((x - startX) / delta);
}

uint GridImp::getYIndex(real y) const
{
    return lround((y - startY) / delta);
}

uint GridImp::getZIndex(real z) const
{
	return lround((z - startZ) / delta);
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

SPtr<GridStrategy> GridImp::getGridStrategy() const
{
    return gridStrategy;
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
    for (uint i = 0; i < this->getSize(); i++)
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

        const char type2 = this->field.getFieldEntry(i);

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
