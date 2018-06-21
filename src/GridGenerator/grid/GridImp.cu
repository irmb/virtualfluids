#include "GridImp.h"

#include "GridGenerator/global.h"
#include <stdio.h>
#include <time.h>

#include <sstream>


#include <GridGenerator/utilities/math/Math.h>
#include "distributions/Distribution.h"

#include <GridGenerator/geometries/Vertex/Vertex.h>
#include <GridGenerator/geometries/Triangle/Triangle.h>
#include <GridGenerator/geometries/TriangularMesh/TriangularMesh.h>
#include <GridGenerator/geometries/TriangularMesh/TriangularMeshStrategy.h>

#include <GridGenerator/geometries/BoundingBox/BoundingBox.h>

#include <GridGenerator/grid/NodeValues.h>
#include <GridGenerator/grid/distributions/Distribution.h>

#include <GridGenerator/grid/GridStrategy/GridStrategy.h>
#include <utilities/logger/Logger.h>
#include "GridInterface.h"

#include "geometries/Object.h"
#include "Field.h"

#include "io/GridVTKWriter/GridVTKWriter.h"



CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];


HOST GridImp::GridImp(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution distribution, uint level) 
: object(object), startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta), gridStrategy(gridStrategy), distribution(distribution), level(level)
{
    initalNumberOfNodesAndSize();
}

HOST SPtr<GridImp> GridImp::makeShared(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level)
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

	this->numberOfSolidBoundaryNodes = 0;
}


HOST void GridImp::inital()
{
    field = Field(gridStrategy, size);
    field.allocateMemory();
    gridStrategy->allocateGridMemory(shared_from_this());

    gridStrategy->initalNodesToOutOfGrid(shared_from_this());

    TriangularMesh* triangularMesh = dynamic_cast<TriangularMesh*>(object);
    if (triangularMesh)
        triangularMeshDiscretizationStrategy->discretize(triangularMesh, this, FLUID, INVALID_OUT_OF_GRID);
    else
        gridStrategy->findInnerNodes(shared_from_this());

    //gridStrategy->findStopperNodes(shared_from_this());
	gridStrategy->findEndOfGridStopperNodes(shared_from_this());

    *logging::out << logging::Logger::INFO_INTERMEDIATE
        << "Grid created: " << "from (" << this->startX << ", " << this->startY << ", " << this->startZ << ") to (" << this->endX << ", " << this->endY << ", " << this->endZ << ")\n"
        << "nodes: " << this->nx << " x " << this->ny << " x " << this->nz << " = " << this->size << "\n";
}

HOSTDEVICE void GridImp::initalNodeToOutOfGrid(uint index)
{
    this->field.setFieldEntryToInvalidOutOfGrid(index);
}

HOST void GridImp::freeMemory()
{
    gridStrategy->freeMemory(shared_from_this());
}

HOST GridImp::GridImp()
{
    //printf("Constructor\n");
    //this->print();
}

HOST GridImp::~GridImp()
{
    //printf("Destructor\n");
    //this->print();
}


HOSTDEVICE void GridImp::findInnerNode(uint index)
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

bool GridImp::isInside(const Cell& cell) const
{
    return object->isCellInObject(cell);
}

////TODO: check where the fine grid starts (0.25 or 0.75) and if even or odd-cell is needed
// *--*--*--*
// |  |  |  |
// *--*--*--*
//  0  1  2
HOSTDEVICE Cell GridImp::getOddCellFromIndex(uint index) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    const uint xIndex = getXIndex(x);
    const uint yIndex = getYIndex(y);
    const uint zIndex = getZIndex(z);

    const real xCellStart = xIndex % 2 != 0 ? x : x - this->delta;
    const real yCellStart = yIndex % 2 != 0 ? y : y - this->delta;
    const real zCellStart = zIndex % 2 != 0 ? z : z - this->delta;
    return Cell(xCellStart, yCellStart, zCellStart, delta);
}

HOSTDEVICE void GridImp::findStopperNode(uint index) // deprecated
{
    if(isValidEndOfGridStopper(index))
        this->field.setFieldEntryToStopperOutOfGrid(index);

    if (isValidSolidStopper(index))
        this->field.setFieldEntry(index, STOPPER_SOLID);
}

HOSTDEVICE void GridImp::findEndOfGridStopperNode(uint index)
{
	if (isValidEndOfGridStopper(index))
		this->field.setFieldEntryToStopperOutOfGrid(index);

	if (isValidEndOfGridBoundaryStopper(index))
		this->field.setFieldEntryToStopperOutOfGridBoundary(index);
}

HOSTDEVICE void GridImp::findSolidStopperNode(uint index)
{
	if (isValidSolidStopper(index))
		this->field.setFieldEntry(index, STOPPER_SOLID);
}

HOSTDEVICE void GridImp::findBoundarySolidNode(uint index)
{
	if (shouldBeBoundarySolidNode(index)) 
	{
		this->field.setFieldEntry(index, BC_SOLID);
		this->qIndices[index] = this->numberOfSolidBoundaryNodes++;
		//grid->setNumberOfSolidBoundaryNodes(grid->getNumberOfSolidBoundaryNodes() + 1);
	}
}

HOSTDEVICE void GridImp::removeOddBoundaryCellNode(uint index)
{
    Cell cell = getOddCellFromIndex(index);
    if (isOutSideOfGrid(cell))
        return;
    if (contains(cell, INVALID_OUT_OF_GRID))
        setNodeTo(cell, INVALID_OUT_OF_GRID);
}

HOSTDEVICE bool GridImp::isOutSideOfGrid(Cell &cell) const
{
    for (const auto point : cell) {
        if (point.x < startX || point.x > endX
            || point.y < startY || point.y > endY
            || point.z < startZ || point.z > endZ)
            return true;
    }
    return false;
}

HOSTDEVICE bool GridImp::contains(Cell &cell, char type) const
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

HOSTDEVICE void GridImp::setNodeTo(Cell &cell, char type)
{
    for (const auto point : cell) {
		uint index = transCoordToIndex(point.x, point.y, point.z);
		if (index == INVALID_INDEX)
			continue;
		field.setFieldEntry(index, type);
    }
}

HOSTDEVICE void GridImp::setNodeTo(uint index, char type)
{
	if( index != INVALID_INDEX )
		field.setFieldEntry(index, type);
}

HOSTDEVICE bool GridImp::isNode(uint index, char type) const
{
    if( index != INVALID_INDEX )
		return field.is(index, type);
}

HOSTDEVICE bool GridImp::isValidEndOfGridStopper(uint index) const
{
	// Lenz: also includes corner stopper nodes
	if (!this->field.is(index, INVALID_OUT_OF_GRID))
		return false;

	return hasNeighborOfType(index, FLUID);

	//previous version of Sören P.
    //return this->field.is(index, OUT_OF_GRID) && (nodeInNextCellIs(index, FLUID) || nodeInNextCellIs(index, FLUID_CFF))
    //    || this->field.is(index, OUT_OF_GRID) && (nodeInPreviousCellIs(index, FLUID) || nodeInPreviousCellIs(index, FLUID_CFF));
}

HOSTDEVICE bool GridImp::isValidEndOfGridBoundaryStopper(uint index) const
{
	// Lenz: also includes corner stopper nodes
	if (!this->field.is(index, FLUID))
		return false;

	return ! hasAllNeighbors(index);

	//previous version of Sören P.
    //return this->field.is(index, OUT_OF_GRID) && (nodeInNextCellIs(index, FLUID) || nodeInNextCellIs(index, FLUID_CFF))
    //    || this->field.is(index, OUT_OF_GRID) && (nodeInPreviousCellIs(index, FLUID) || nodeInPreviousCellIs(index, FLUID_CFF));
}

HOSTDEVICE bool GridImp::isValidSolidStopper(uint index) const
{
	// Lenz: also includes corner stopper nodes
	if (!this->field.is(index, INVALID_SOLID))
		return false;

	return hasNeighborOfType(index, FLUID);

	//previous version of Sören P.
	//return this->field.is(index, SOLID) && (nodeInNextCellIs(index, FLUID) || nodeInNextCellIs(index, FLUID_CFF))
 //       || this->field.is(index, SOLID) && (nodeInPreviousCellIs(index, FLUID) || nodeInPreviousCellIs(index, FLUID_CFF));
}

HOSTDEVICE bool GridImp::shouldBeBoundarySolidNode(uint index) const
{
	if (!this->field.is(index, FLUID))
		return false;

	return hasNeighborOfType(index, STOPPER_SOLID);
}

HOSTDEVICE bool GridImp::hasAllNeighbors(uint index) const
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

HOSTDEVICE bool GridImp::hasNeighborOfType(uint index, char type) const
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

    //real x, y, z;
    //this->transIndexToCoords(index, x, y, z);

    //const real neighborX = x + this->delta > endX ? endX : x + this->delta;
    //const real neighborY = y + this->delta > endY ? endY : y + this->delta;
    //const real neighborZ = z + this->delta > endZ ? endZ : z + this->delta;

    //const real neighborMinusX = x - this->delta < startX ? startX : x - this->delta;
    //const real neighborMinusY = y - this->delta < startY ? startY : y - this->delta;
    //const real neighborMinusZ = z - this->delta < startZ ? startZ : z - this->delta;


    //const uint indexMXY = transCoordToIndex(-neighborX, neighborY, z);
    //const uint indexMYZ = transCoordToIndex(x, -neighborY, neighborZ);
    //const uint indexMXZ = transCoordToIndex(-neighborX, y, neighborZ);

    //const uint indexXMY = transCoordToIndex(neighborX, -neighborY, z);
    //const uint indexYMZ = transCoordToIndex(x, neighborY, -neighborZ);
    //const uint indexXMZ = transCoordToIndex(neighborX, y, -neighborZ);

    //const uint indexMXYMZ = transCoordToIndex(-neighborX, neighborY, -neighborZ);
    //const uint indexMXYZ  = transCoordToIndex(-neighborX, neighborY, neighborZ);
    //const uint indexMXMYZ = transCoordToIndex(-neighborX, -neighborY, neighborZ);
    //const uint indexXMYMZ = transCoordToIndex(neighborX, -neighborY, -neighborZ);
    //const uint indexXMYZ  = transCoordToIndex(neighborX, -neighborY, neighborZ);
    //const uint indexXYMZ  = transCoordToIndex(neighborX, neighborY, -neighborZ);



    //return nodeInNextCellIs(index, type) || nodeInPreviousCellIs(index, type) || indexMXY || indexMYZ || indexMXZ || indexXMY || indexYMZ || indexXMZ 
    //    || indexMXYMZ
    //    || indexMXYZ
    //    || indexMXMYZ
    //    || indexXMYMZ
    //    || indexXMYZ
    //    || indexXYMZ;

}

HOSTDEVICE bool GridImp::nodeInNextCellIs(int index, char type) const
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

HOSTDEVICE bool GridImp::nodeInPreviousCellIs(int index, char type) const
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

HOSTDEVICE bool GridImp::nodeInCellIs(Cell& cell, char type) const
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


HOST void GridImp::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ)
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

HOSTDEVICE uint GridImp::transCoordToIndex(const real &x, const real &y, const real &z) const
{
    const uint xIndex = getXIndex(x);
    const uint yIndex = getYIndex(y);
    const uint zIndex = getZIndex(z);

	if (xIndex >= nx || yIndex >= ny || zIndex >= nz)
        return INVALID_INDEX;

    return xIndex + nx * (yIndex + ny * zIndex);
}

HOSTDEVICE void GridImp::transIndexToCoords(uint index, real &x, real &y, real &z) const
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

HOST uint GridImp::getLevel(real startDelta) const
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

HOST uint GridImp::getLevel() const
{
    return this->level;
}

HOST void GridImp::setTriangularMeshDiscretizationStrategy(TriangularMeshDiscretizationStrategy* triangularMeshDiscretizationStrategy)
{
    this->triangularMeshDiscretizationStrategy = triangularMeshDiscretizationStrategy;
}

HOST uint GridImp::getNumberOfSolidBoundaryNodes() const
{
	return this->numberOfSolidBoundaryNodes;
}

HOST void GridImp::setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes)
{
	if (numberOfSolidBoundaryNodes >= 0 && numberOfSolidBoundaryNodes < INVALID_INDEX)
		this->numberOfSolidBoundaryNodes = numberOfSolidBoundaryNodes;
}

HOST real GridImp::getQValue(const uint index, const uint dir) const
{
	const int qIndex = dir * this->numberOfSolidBoundaryNodes + this->qIndices[index];

	return this->qValues[qIndex];
}

// --------------------------------------------------------- //
//                  Set Sparse Indices                       //
// --------------------------------------------------------- //

HOST void GridImp::findSparseIndices(SPtr<Grid> fineGrid)
{
    this->gridStrategy->findSparseIndices(shared_from_this(), std::static_pointer_cast<GridImp>(fineGrid));
}


HOST void GridImp::updateSparseIndices()
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

HOSTDEVICE void GridImp::setNeighborIndices(uint index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    neighborIndexX[index] = -1;
    neighborIndexY[index] = -1;
    neighborIndexZ[index] = -1;

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

    this->neighborIndexX[index] = neighborX;
    this->neighborIndexY[index] = neighborY;
    this->neighborIndexZ[index] = neighborZ;
}

HOSTDEVICE void GridImp::setStopperNeighborCoords(uint index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    if (vf::Math::lessEqual(x + delta, endX) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x + delta, y, z)))
        neighborIndexX[index] = getSparseIndex(x + delta, y, z);

    if (vf::Math::lessEqual(y + delta, endY) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x, y + delta, z)))
        neighborIndexY[index] = getSparseIndex(x, y + delta, z);

    if (vf::Math::lessEqual(z + delta, endZ) && !this->field.isInvalidOutOfGrid(this->transCoordToIndex(x, y, z + delta)))
        neighborIndexZ[index] = getSparseIndex(x, y, z + delta);
}

HOSTDEVICE void GridImp::getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const
{
    real coords[3] = { x, y, z };
    neighborX = getNeighborCoord(periodicityX, startX, coords, 0);
    neighborY = getNeighborCoord(periodicityY, startY, coords, 1);
    neighborZ = getNeighborCoord(periodicityZ, startZ, coords, 2);
}

HOSTDEVICE real GridImp::getNeighborCoord(bool periodicity, real startCoord, real coords[3], int direction) const
{
    if (periodicity)
    {
        real neighborCoords[3] = {coords[0], coords[1] , coords[2] };
        neighborCoords[direction] = neighborCoords[direction] + delta;
        const int neighborIndex = this->transCoordToIndex(neighborCoords[0], neighborCoords[1], neighborCoords[2]);

        if(!field.isStopperOutOfGrid(neighborIndex) && !field.is(neighborIndex, STOPPER_OUT_OF_GRID_BOUNDARY) )
            return coords[direction] + delta;

        return getFirstFluidNode(coords, direction, startCoord);
    }
    
    return coords[direction] + delta;
}


HOSTDEVICE real GridImp::getLastFluidNode(real coords[3], int direction, real startCoord) const
{
    coords[direction] = startCoord;
    int index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    while (!field.isFluid(index))
    {
        coords[direction] -= delta;
        index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    }
    return coords[direction];
}

HOSTDEVICE real GridImp::getFirstFluidNode(real coords[3], int direction, real startCoord) const
{
    coords[direction] = startCoord;
    int index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    while (!field.isFluid(index))
    {
        coords[direction] += delta;
        index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    }
    return coords[direction];
}


HOSTDEVICE int GridImp::getSparseIndex(const real &x, const real &y, const real &z) const
{
    const int matrixIndex = transCoordToIndex(x, y, z);
    return sparseIndices[matrixIndex];
}


// --------------------------------------------------------- //
//                    Find Interface                         //
// --------------------------------------------------------- //
HOST void GridImp::findGridInterface(SPtr<Grid> finerGrid, LbmOrGks lbmOrGks)
{
    gridStrategy->findGridInterface(shared_from_this(), std::static_pointer_cast<GridImp>(finerGrid), lbmOrGks);
}

HOSTDEVICE void GridImp::findGridInterfaceCF(uint index, GridImp& finerGrid, LbmOrGks lbmOrGks)
{
	if(lbmOrGks == LBM)
		gridInterface->findInterfaceCF(index, this, &finerGrid);
	else if (lbmOrGks == GKS)
		gridInterface->findInterfaceCF_GKS(index, this, &finerGrid);
}

HOSTDEVICE void GridImp::findGridInterfaceFC(uint index, GridImp& finerGrid)
{
    gridInterface->findInterfaceFC(index, this, &finerGrid);
}

HOSTDEVICE void GridImp::findOverlapStopper(uint index, GridImp& finerGrid)
{
    gridInterface->findOverlapStopper(index, this, &finerGrid);
}

// --------------------------------------------------------- //
//                    Mesh Triangle                          //
// --------------------------------------------------------- //
HOST void GridImp::mesh(Object* object)
{
    TriangularMesh* triangularMesh = dynamic_cast<TriangularMesh*>(object);
    if (triangularMesh)
        triangularMeshDiscretizationStrategy->discretize(triangularMesh, this, INVALID_SOLID, FLUID);
    else
        gridStrategy->findInnerNodes(shared_from_this()); //TODO: adds INNERTYPE AND OUTERTYPE to findInnerNodes 
		//new method for geometric primitives (not cell based) to be implemented

	//gridStrategy->findStopperNodes(shared_from_this()); //deprecated
	gridStrategy->findSolidStopperNodes(shared_from_this());
	gridStrategy->findBoundarySolidNodes(shared_from_this());
}


HOST void GridImp::mesh(TriangularMesh &triangularMesh)
{
    const clock_t begin = clock();

    gridStrategy->mesh(shared_from_this(), triangularMesh);

    const clock_t end = clock();
    const real time = (real)(real(end - begin) / CLOCKS_PER_SEC);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "time grid generation: " << time << "s\n";
}

HOSTDEVICE void GridImp::mesh(Triangle &triangle)
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
                const int value = triangle.isUnderFace(point);
                //setDebugPoint(index, value);

                if (value == Q)
                    calculateQs(point, triangle);
            }
        }
    }
}



HOST void GridImp::findQs(Object* object) //TODO: enable qs for primitive objects
{
    TriangularMesh* triangularMesh = dynamic_cast<TriangularMesh*>(object);
    if (triangularMesh)
        findQs(*triangularMesh);
}

HOST void GridImp::findQs(TriangularMesh &triangularMesh)
{
    const clock_t begin = clock();

	gridStrategy->allocateQs(shared_from_this());
    gridStrategy->findQs(shared_from_this(), triangularMesh);

    const clock_t end = clock();
    const real time = (real)(real(end - begin) / CLOCKS_PER_SEC);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "time finding qs: " << time << "s\n";
}

HOSTDEVICE void GridImp::findQs(Triangle &triangle)
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

                if(this->field.is(index, BC_SOLID)) //TODO: not working with thin walls
                {
                    //calculateQs(point, triangle);
					calculateQs(index, point, triangle);
				}
            }
        }
    }
}

HOSTDEVICE void GridImp::setDebugPoint(uint index, int pointValue)
{
    if (field.isInvalidCoarseUnderFine(index) && pointValue == INVALID_SOLID)
        field.setFieldEntry(index, pointValue);

    if(!field.isInvalidSolid(index) && !field.isQ(index) && !field.isInvalidCoarseUnderFine(index) && pointValue != 3 && pointValue != 2)
        field.setFieldEntry(index, pointValue);
}

HOSTDEVICE void GridImp::calculateQs(const Vertex &point, const Triangle &triangle) const
{
    Vertex pointOnTriangle, direction;
    //VertexInteger solid_node;
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

		//real lengthDirection = sqrt(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z);
		subdistance /= /*lengthDirection **/ this->delta;

        if (error != 0 && subdistance < 1.0 && subdistance > 0.0)
        {
            //solid_node = VertexInteger(actualPoint.x + direction.x, actualPoint.y + direction.y, actualPoint.z + direction.z);
            distribution.f[i*size + transCoordToIndex(point.x, point.y, point.z)] = subdistance;
            //printf("Q%d %d: %2.8f \n", i, grid.transCoordToIndex(actualPoint), grid.d.f[index]);
        } /*else
            distribution.f[i*size + transCoordToIndex(point.x, point.y, point.z)] = -1.0;*/
    }
}


HOSTDEVICE void GridImp::calculateQs(const uint index, const Vertex &point, const Triangle &triangle) const
{
	Vertex pointOnTriangle, direction;
	//VertexInteger solid_node;
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

		if (error != 0 && vf::Math::lessEqual(subdistance, 1.0) && vf::Math::greaterEqual(subdistance, 0.0))
		{
			this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] = subdistance;
		}
		//else
		//{
		//	uint neighborIndex = this->transCoordToIndex( point.x + direction.x * this->delta,
		//		                                          point.y + direction.y * this->delta, 
		//		                                          point.z + direction.z * this->delta );
		//	if (neighborIndex == INVALID_INDEX) continue;

		//	if( this->field.is(neighborIndex, STOPPER_SOLID) )
		//		this->qValues[i*this->numberOfSolidBoundaryNodes + this->qIndices[index]] = 1.0;
		//}
	}
}


// --------------------------------------------------------- //
//                        Getter                             //
// --------------------------------------------------------- //
HOSTDEVICE int GridImp::getSparseIndex(uint matrixIndex) const
{
    return this->sparseIndices[matrixIndex];
}

HOST real* GridImp::getDistribution() const
{
    return this->distribution.f;
}

HOST int* GridImp::getDirection() const
{
    return this->distribution.dirs;
}

HOST int GridImp::getStartDirection() const
{
    return this->distribution.dir_start;
}

HOST int GridImp::getEndDirection() const
{
    return this->distribution.dir_end;
}

HOSTDEVICE BoundingBox GridImp::getBoundingBoxOnNodes(Triangle &triangle) const
{
    real minX, maxX, minY, maxY, minZ, maxZ;
    triangle.setMinMax(minX, maxX, minY, maxY, minZ, maxZ);
    const Vertex minOnNodes = getMinimumOnNode(Vertex(minX, minY, minZ));
    const Vertex maxOnNodes = getMaximumOnNode(Vertex(maxX, maxY, maxZ));

    return BoundingBox(minOnNodes.x, maxOnNodes.x, minOnNodes.y, maxOnNodes.y, minOnNodes.z, maxOnNodes.z);
}

HOSTDEVICE Vertex GridImp::getMinimumOnNode(Vertex exact) const
{
    const real minX = getMinimumOnNodes(exact.x, vf::Math::getDecimalPart(startX), delta);
    const real minY = getMinimumOnNodes(exact.y, vf::Math::getDecimalPart(startY), delta);
    const real minZ = getMinimumOnNodes(exact.z, vf::Math::getDecimalPart(startZ), delta);
    return Vertex(minX, minY, minZ);
}

HOSTDEVICE real GridImp::getMinimumOnNodes(const real& minExact, const real& decimalStart, const real& delta)
{
    real minNode = ceil(minExact - 1.0);
    minNode += decimalStart;
    while (minNode > minExact)
        minNode -= delta;

    while (minNode + delta < minExact)
        minNode += delta;
    return minNode;
}

HOSTDEVICE Vertex GridImp::getMaximumOnNode(Vertex exact) const
{
    const real maxX = getMaximumOnNodes(exact.x, vf::Math::getDecimalPart(startX), delta);
    const real maxY = getMaximumOnNodes(exact.y, vf::Math::getDecimalPart(startY), delta);
    const real maxZ = getMaximumOnNodes(exact.z, vf::Math::getDecimalPart(startZ), delta);
    return Vertex(maxX, maxY, maxZ);
}

HOSTDEVICE real GridImp::getMaximumOnNodes(const real& maxExact, const real& decimalStart, const real& delta)
{
    real maxNode = ceil(maxExact - 1.0);
    maxNode += decimalStart;

    while (maxNode <= maxExact)
        maxNode += delta;
    return maxNode;
}

HOSTDEVICE uint GridImp::getXIndex(real x) const
{
    return lround((x - startX) / delta);
	//return int((x - startX) / delta);
}

HOSTDEVICE uint GridImp::getYIndex(real y) const
{
    return lround((y - startY) / delta);
	//return int((y - startY) / delta);
}

HOSTDEVICE uint GridImp::getZIndex(real z) const
{
	return lround((z - startZ) / delta);
	//return int((z - startZ) / delta);
}

HOSTDEVICE real GridImp::getDelta() const
{
    return delta;
}

HOSTDEVICE uint GridImp::getSize() const
{
    return this->size;
}

HOSTDEVICE uint GridImp::getSparseSize() const
{
    return this->sparseSize;
}

HOSTDEVICE Field GridImp::getField() const
{
    return this->field;
}

char GridImp::getFieldEntry(uint index) const
{
    return this->field.getFieldEntry(index);
}

HOSTDEVICE void GridImp::setFieldEntry(uint matrixIndex, char type)
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

uint* GridImp::getFC_coarse() const
{
    return this->gridInterface->fc.coarse;
}

uint* GridImp::getFC_fine() const
{
    return this->gridInterface->fc.fine;
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
        gridInterfaceList[i] = oldGridInterfaceList[i] + 1;
}

#define GEOFLUID 19
#define GEOSOLID 16

HOST void GridImp::getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *geo) const
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

        const uint neighborXIndex = uint(this->neighborIndexX[i] + 1);
        const uint neighborYIndex = uint(this->neighborIndexY[i] + 1);
        const uint neighborZIndex = uint(this->neighborIndexZ[i] + 1);

        const char type2 = this->field.getFieldEntry(i);

        const uint type = uint(this->field.isFluid(i) ? GEOFLUID : GEOSOLID);

        xCoords[nodeNumber + 1] = x;
        yCoords[nodeNumber + 1] = y;
        zCoords[nodeNumber + 1] = z;

        neighborX[nodeNumber + 1] = neighborXIndex;
        neighborY[nodeNumber + 1] = neighborYIndex;
        neighborZ[nodeNumber + 1] = neighborZIndex;
        geo[nodeNumber + 1] = type;
        nodeNumber++;
    }
}

void GridImp::print() const
{
    printf("min: (%2.4f, %2.4f, %2.4f), max: (%2.4f, %2.4f, %2.4f), size: %d, delta: %2.4f\n", startX, startY, startZ,
           endX, endY, endZ, size, delta);
    if(this->gridInterface)
        this->gridInterface->print();
}
