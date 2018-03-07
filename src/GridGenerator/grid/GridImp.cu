#include "GridImp.h"

#include "GridGenerator/global.h"
#include <stdio.h>
#include <time.h>

#include <sstream>


#include <GridGenerator/utilities/math/CudaMath.cuh>
#include "distributions/Distribution.h"

#include <GridGenerator/geometries/Vertex/Vertex.cuh>
#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>

#include <GridGenerator/grid/NodeValues.h>
#include <GridGenerator/grid/distributions/Distribution.h>

#include <GridGenerator/grid/GridStrategy/GridStrategy.h>
#include <utilities/logger/Logger.h>
#include "GridInterface.h"

#include "geometries/Object.h"



CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

HOST GridImp::GridImp(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution distribution)
    : startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta), distribution(distribution),
    gridInterface(nullptr), neighborIndexX(nullptr), neighborIndexY(nullptr), neighborIndexZ(nullptr), sparseIndices(nullptr), gridStrategy(gridStrategy)
{
    initalNumberOfNodesAndSize();
}

HOST GridImp::GridImp(Object* object, real delta, SPtr<GridStrategy> gridStrategy, Distribution distribution) : object(object), delta(delta), gridStrategy(gridStrategy), distribution(distribution)
{
    initalBoundingBoxStartValues();
    initalNumberOfNodesAndSize();
}

HOST void GridImp::initalBoundingBoxStartValues()
{
    startX = object->getX1Minimum();
    startY = object->getX2Minimum();
    startZ = object->getX3Minimum();

    endX = object->getX1Maximum();
    endY = object->getX2Maximum();
    endZ = object->getX3Maximum();
}

void GridImp::initalNumberOfNodesAndSize()
{
    startX -= delta; // +1 stopper node
    startY -= delta; // +1 stopper node
    startZ -= delta; // +1 stopper node

    endX += delta; // +1 stopper node
    endY += delta; // +1 stopper node
    endZ += delta; // +1 stopper node

    const real length = endX - startX;
    const real width = endY - startY;
    const real height = endZ - startZ;

    nx = int((length + delta) / delta);
    ny = int((width + delta) / delta);
    nz = int((height + delta) / delta);

    this->size = nx * ny * nz;
    this->sparseSize = size;
    distribution.setSize(size);
}

HOST SPtr<GridImp> GridImp::makeShared(Object* object, real delta, SPtr<GridStrategy> gridStrategy, Distribution d)
{
    SPtr<GridImp> grid(new GridImp(object, delta, gridStrategy, d));
    return grid;
}

HOST SPtr<GridImp> GridImp::makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, 
                                        SPtr<GridStrategy> gridStrategy, Distribution d)
{
    SPtr<GridImp> grid(new GridImp(startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, d));
    return grid;
}

HOST void GridImp::inital()
{
    field = Field(gridStrategy, size);
    field.allocateMemory();
    gridStrategy->allocateGridMemory(shared_from_this());
    gridStrategy->initalNodes(shared_from_this());
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


HOSTDEVICE void GridImp::findInnerNode(uint index)
{
    this->sparseIndices[index] = index;

    const Cell cell = getOddCellFromIndex(index);
    if (object->isCellInObject(cell))
        this->field.setFieldEntryToFluid(index);
    else
        this->field.setFieldEntryToOutOfGrid(index);
}

//TODO: check where the fine grid starts (0.25 or 0.75) and if even or odd-cell is needed
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


uint GridImp::getXIndex(real x) const
{
   return uint((x - startX) / delta);
}

uint GridImp::getYIndex(real y) const
{
   return uint((y - startY) / delta);
}

uint GridImp::getZIndex(real z) const
{
    return uint((z - startZ) / delta);
}

HOSTDEVICE void GridImp::findStopperNode(uint index)
{
    if(isValidEndOfGridStopper(index))
        this->field.setFieldEntryToStopperEndOfGrid(index);

    if (isValidStartOfGridStopper(index))
        this->field.setFieldEntryToStopperEndOfGrid(index);
}

HOST void GridImp::mesh(Geometry &geometry)
{
    clock_t begin = clock();

    gridStrategy->mesh(shared_from_this(), geometry);

    clock_t end = clock();
    real time = (real)(real(end - begin) / CLOCKS_PER_SEC);

    *logging::out << logging::Logger::INTERMEDIATE << "time grid generation: " + SSTR(time) + "s\n";
}

HOST void GridImp::freeMemory()
{
    gridStrategy->freeMemory(shared_from_this());
}

HOST void GridImp::findGridInterface(SPtr<Grid> finerGrid)
{
    gridStrategy->findGridInterface(shared_from_this(), std::static_pointer_cast<GridImp>(finerGrid));
}

HOSTDEVICE void GridImp::findGridInterfaceCF(uint index, GridImp& finerGrid)
{
    gridInterface->findInterfaceCF(index, this, &finerGrid);
}

HOSTDEVICE void GridImp::findGridInterfaceFC(uint index, GridImp& finerGrid)
{
    gridInterface->findInterfaceFC(index, this, &finerGrid);
}

HOSTDEVICE void GridImp::findOverlapStopper(uint index, GridImp& finerGrid)
{
    gridInterface->findOverlapStopper(index, this, &finerGrid);
}


HOST void GridImp::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ)
{
    this->periodicityX = periodicityX;
    this->periodicityY = periodicityY;
    this->periodicityZ = periodicityZ;
}



HOSTDEVICE void GridImp::setFieldEntry(const Vertex &v, char val)
{
    const uint index = transCoordToIndex(v);
    this->field.setFieldEntry(index, val);
}

HOSTDEVICE char GridImp::getFieldEntry(const Vertex &v) const
{
    const uint index = transCoordToIndex(v);
    return this->field.getFieldEntry(index);
}

HOSTDEVICE int GridImp::transCoordToIndex(const real &x, const real &y, const real &z) const
{
    return transCoordToIndex(Vertex(x,y,z));
}

HOSTDEVICE int GridImp::transCoordToIndex(const Vertex &v) const
{
    const int x = int((v.x - startX) / delta);
    const int y = int((v.y - startY) / delta);
    const int z = int((v.z - startZ) / delta);

#ifdef DEBUG
    if (x < 0 || y < 0 || z < 0 || uint(x) > nx || uint(y) > ny || uint(z) > nz)
    {
        //printf( "Function: transCoordToIndex. Coordinates are out of range and cannot calculate the index. Exit Program!\n");
        /* exit(1);*/
        return -1;
    }
#endif

    return x + nx * (y + ny * z);
}

HOSTDEVICE void GridImp::transIndexToCoords(int index, real &x, real &y, real &z) const
{
#ifdef DEBUG
    if (index < 0 || index >= (int)size)
    {
        printf("Function: transIndexToCoords. GridImp Index: %d, size: %d. Exit Program!\n", index, size); /*exit(1);*/ 
    }
#endif
    x = index % nx;
    y = (index / nx) % ny;
    z = ((index / nx) / ny) % nz;

    x = (x * delta) + startX;
    y = (y * delta) + startY;
    z = (z * delta) + startZ;
}


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

HOSTDEVICE void GridImp::setDebugPoint(const Vertex &point, const int pointValue)
{
    if (getFieldEntry(point) == INVALID_NODE && pointValue == SOLID)
        setFieldEntry(point, pointValue);

    if (getFieldEntry(point) != SOLID && getFieldEntry(point) != Q && getFieldEntry(point) != INVALID_NODE && pointValue
        != 3 && pointValue != 2)
        setFieldEntry(point, pointValue);
}

HOSTDEVICE bool GridImp::isOutOfRange(const Vertex &v) const
{
    return v.x < startX || v.y < startY || v.z < startZ || v.x > endX || v.y > endY || v.z > endZ;
}

HOSTDEVICE void GridImp::meshTriangle(const Triangle &triangle)
{
    BoundingBox<real> box = BoundingBox<real>::makeRealNodeBox(triangle, delta);

    for (real x = box.minX; x <= box.maxX; x += delta)
    {
        for (real y = box.minY; y <= box.maxY; y += delta)
        {
            for (real z = box.minZ; z <= box.maxZ; z += delta)
            {
                Vertex point(x, y, z);
                if (isOutOfRange(point))
                    continue;
                const int value = triangle.isUnderFace(point);
                setDebugPoint(point, value);

                if (value == Q)
                    calculateQs(point, triangle);
            }
        }
    }
}

HOSTDEVICE void GridImp::calculateQs(const Vertex &point, const Triangle &triangle)
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

        if (error != 0 && subdistance <= 1.0f)
        {
            //solid_node = VertexInteger(actualPoint.x + direction.x, actualPoint.y + direction.y, actualPoint.z + direction.z);
            distribution.f[i*size + transCoordToIndex(point)] = subdistance;
            //printf("Q%d %d: %2.8f \n", i, grid.transCoordToIndex(actualPoint), grid.d.f[index]);
        }
    }
}

HOSTDEVICE void GridImp::setNeighborIndices(const int &index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    real neighborX, neighborY, neighborZ;
    this->getNeighborCoords(neighborX, neighborY, neighborZ, x, y, z);

    neighborIndexX[index] = uint(transCoordToIndex(neighborX, y, z));
    neighborIndexY[index] = uint(transCoordToIndex(x, neighborY, z));
    neighborIndexZ[index] = uint(transCoordToIndex(x, y, neighborZ));
}

HOSTDEVICE void GridImp::setInvalidNode(const int &index, bool &invalidNodeFound)
{
    if (field.isSolid(index))
        return;

    if (!field.isInvalid(index) && isNeighborInvalid(index))
    {
        field.setFieldEntryToInvalid(index);
        invalidNodeFound = true;
    }
}

HOSTDEVICE bool GridImp::isNeighborInvalid(const int &index) const
{
    return (field.isInvalid(neighborIndexX[index]) || field.isInvalid(neighborIndexY[index]) || field.isInvalid(neighborIndexZ[index]));
}

HOSTDEVICE void GridImp::findNeighborIndex(int index)
{
    neighborIndexX[index] = -1;
    neighborIndexY[index] = -1;
    neighborIndexZ[index] = -1;

    if (this->field.isStopperEndOfGrid(index) || this->field.isStopperOverlapGrid(index))
    {
        this->setStopperNeighborCoords(index);
        return;
    }

    if (this->sparseIndices[index] == -1)
        return;

    real x, y, z;
    this->transIndexToCoords(index, x, y, z);
    real neighborXCoord, neighborYCoord, neighborZCoord;
    getNeighborCoords(neighborXCoord, neighborYCoord, neighborZCoord, x, y, z);
    const int neighborX = getNeighborIndex(neighborXCoord, y, z);
    const int neighborY = getNeighborIndex(x, neighborYCoord, z);
    const int neighborZ = getNeighborIndex(x, y, neighborZCoord);

    this->neighborIndexX[index] = neighborX;
    this->neighborIndexY[index] = neighborY;
    this->neighborIndexZ[index] = neighborZ;
}

HOSTDEVICE void GridImp::setStopperNeighborCoords(int index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    if (CudaMath::lessEqual(x + delta, endX) && (this->field.isStopperEndOfGrid(this->transCoordToIndex(x + delta, y, z)) || this->field.isFluid(this->transCoordToIndex(x + delta, y, z))))
        neighborIndexX[index] = getNeighborIndex(x + delta, y, z);

    if (CudaMath::lessEqual(y + delta, endY) && (this->field.isStopperEndOfGrid(this->transCoordToIndex(x, y + delta, z)) || this->field.isFluid(this->transCoordToIndex(x, y + delta, z))))
        neighborIndexY[index] = getNeighborIndex(x, y + delta, z);

    if (CudaMath::lessEqual(z + delta, endZ) && (this->field.isStopperEndOfGrid(this->transCoordToIndex(x, y, z + delta)) || this->field.isFluid(this->transCoordToIndex(x, y, z + delta))))
        neighborIndexZ[index] = getNeighborIndex(x, y, z + delta);
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
        if(!field.isStopperEndOfGrid(neighborIndex))
            return coords[direction] + delta;

        return getFirstFluidNode(coords, direction, startCoord);
    }
    else
        return coords[direction] + delta;
}


HOSTDEVICE real GridImp::getFirstFluidNode(real coords[3], int direction, real startCoord) const
{
    coords[direction] = startCoord;
    int index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    while (field.isOutOfGrid(index))
    {
        coords[direction] += delta;
        index = this->transCoordToIndex(coords[0], coords[1], coords[2]);
    }
    return coords[direction];
}

HOSTDEVICE bool GridImp::isValidStartOfGridStopper(uint index) const
{
    return this->field.is(index, OUT_OF_GRID) && (nodeInNextCellIs(index, FLUID) || nodeInNextCellIs(index, FLUID_CFF));
}


HOSTDEVICE bool GridImp::isValidEndOfGridStopper(uint index) const
{
    return this->field.is(index, OUT_OF_GRID) && (nodeInPreviousCellIs(index, FLUID) || nodeInPreviousCellIs(index, FLUID_CFF));
}

HOSTDEVICE bool GridImp::isOverlapStopper(uint index) const
{
    return this->field.isFluid(index) && nodeInNextCellIs(index, INVALID_NODE);
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

    const bool isInvalidNeighborX = this->field.is(indexX, type);
    const bool isInvalidNeighborY = this->field.is(indexY, type);
    const bool isInvalidNeighborXY  = this->field.is(indexXY, type);
    const bool isInvalidNeighborZ   = this->field.is(indexZ, type);
    const bool isInvalidNeighborYZ = this->field.is(indexYZ, type);
    const bool isInvalidNeighborXZ = this->field.is(indexXZ, type);
    const bool isInvalidNeighborXYZ = this->field.is(indexXYZ, type);

    return isInvalidNeighborX || isInvalidNeighborY || isInvalidNeighborXY || isInvalidNeighborZ || isInvalidNeighborYZ
        || isInvalidNeighborXZ || isInvalidNeighborXYZ;
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

    const bool isInvalidNeighborX = this->field.is(indexX, type);
    const bool isInvalidNeighborY = this->field.is(indexY, type);
    const bool isInvalidNeighborXY = this->field.is(indexXY, type);
    const bool isInvalidNeighborZ = this->field.is(indexZ, type);
    const bool isInvalidNeighborYZ = this->field.is(indexYZ, type);
    const bool isInvalidNeighborXZ = this->field.is(indexXZ, type);
    const bool isInvalidNeighborXYZ = this->field.is(indexXYZ, type);

    return isInvalidNeighborX || isInvalidNeighborY || isInvalidNeighborXY || isInvalidNeighborZ || isInvalidNeighborYZ
        || isInvalidNeighborXZ || isInvalidNeighborXYZ;
}

HOSTDEVICE int GridImp::getNeighborIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const
{
    const int neighborIndex = transCoordToIndex(expectedX, expectedY, expectedZ);
    return sparseIndices[neighborIndex];
}

//void findForGridInterfaceNewIndex(GridInterface::Interface interface, uint interfaceIndex, int* indices)
//{
//    const uint oldIndex = interface.coarse[interfaceIndex];
//    const uint newIndex = indices[oldIndex];
//    interface.coarse[interfaceIndex] = newIndex;
//}

HOSTDEVICE void GridImp::findForGridInterfaceNewIndexCF(uint index)
{
    const uint oldIndex = gridInterface->cf.coarse[index];
    const uint newIndex = sparseIndices[oldIndex];
    gridInterface->cf.coarse[index] = newIndex;
}

HOSTDEVICE void GridImp::findForGridInterfaceNewIndexFC(uint index)
{
    const uint oldIndex = gridInterface->fc.coarse[index];
    const uint newIndex = sparseIndices[oldIndex];
    gridInterface->fc.coarse[index] = newIndex;
}


HOST void GridImp::removeInvalidNodes()
{
    int removedNodes = 0;
    int newIndex = 0;
    for (uint index = 0; index < size; index++)
    {
        if (this->field.isInvalid(index) || this->field.isOutOfGrid(index))
        {
            sparseIndices[index] = -1;
            removedNodes++;
        } else
        {
            sparseIndices[index] = newIndex;
            newIndex++;
        }
    }
    sparseSize = size - removedNodes;
    printf("new size nodes: %d , delete nodes: %d\n", sparseSize, removedNodes);
}


void GridImp::setCellTo(uint index, char type)
{
    this->field.setFieldEntry(index, type);

    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    const int neighborX = this->transCoordToIndex(x + this->delta, y, z);
    const int neighborY = this->transCoordToIndex(x, y + this->delta, z);
    const int neighborZ = this->transCoordToIndex(x, y, z + this->delta);

    this->field.setFieldEntry(neighborX, type);
    this->field.setFieldEntry(neighborY, type);
    this->field.setFieldEntry(neighborZ, type);

    const int neighborYZ = this->transCoordToIndex(x, y + this->delta, z + this->delta);
    const int neighborXZ = this->transCoordToIndex(x + this->delta, y, z + this->delta);
    const int neighborXY = this->transCoordToIndex(x + this->delta, y + this->delta, z);

    this->field.setFieldEntry(neighborYZ, type);
    this->field.setFieldEntry(neighborXZ, type);
    this->field.setFieldEntry(neighborXY, type);

    const int neighborXYZ = this->transCoordToIndex(x + this->delta, y + this->delta, z + this->delta);

    this->field.setFieldEntry(neighborXYZ, type);
}

char GridImp::getFieldEntry(uint index) const
{
    return this->field.getFieldEntry(index);
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

void GridImp::setStartX(real startX)
{
    this->startX = startX;
    initalNumberOfNodesAndSize();
}

void GridImp::setStartY(real startY)
{
    this->startY = startY;
    initalNumberOfNodesAndSize();
}

void GridImp::setStartZ(real startZ)
{
    this->startZ = startZ;
    initalNumberOfNodesAndSize();
}

void GridImp::setEndX(real endX)
{
    this->endX = endX;
    initalNumberOfNodesAndSize();
}

void GridImp::setEndY(real endY)
{
    this->endY = endY;
    initalNumberOfNodesAndSize();
}

void GridImp::setEndZ(real endZ)
{
    this->endZ = endZ;
    initalNumberOfNodesAndSize();
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
    setGridInterface(iCellCfc, this->gridInterface->cf.coarse, this->gridInterface->cf.numberOfEntries);
    setGridInterface(iCellCff, this->gridInterface->cf.fine, this->gridInterface->cf.numberOfEntries);
    setGridInterface(iCellFcc, this->gridInterface->fc.coarse, this->gridInterface->fc.numberOfEntries);
    setGridInterface(iCellFcf, this->gridInterface->fc.fine, this->gridInterface->fc.numberOfEntries);
}

void GridImp::setGridInterface(uint* gridInterfaceList, const uint* oldGridInterfaceList, uint size)
{
    for (uint i = 0; i < size; i++)
        gridInterfaceList[i] = oldGridInterfaceList[i] + 1;
}

#define GEOFLUID 19
#define GEOSOLID 16

HOST void GridImp::setNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *geo) const
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

        xCoords[nodeNumber + 1] = x;
        yCoords[nodeNumber + 1] = y;
        zCoords[nodeNumber + 1] = z;
        neighborX[nodeNumber + 1] = uint(this->neighborIndexX[i] + 1);
        neighborY[nodeNumber + 1] = uint(this->neighborIndexY[i] + 1);
        neighborZ[nodeNumber + 1] = uint(this->neighborIndexZ[i] + 1);
        geo[nodeNumber + 1] = uint(this->field.isFluid(i) ? GEOFLUID : GEOSOLID);
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

std::string GridImp::toString() const
{
    std::ostringstream oss;
    oss << 
        "min: (" << startX << ", " << startY << ", " << startZ << 
        "), max: " << endX << ", " << endY << ", " << endZ <<
        "), size: " << sparseSize << ", delta: " << delta << "\n";
    return oss.str();
}


