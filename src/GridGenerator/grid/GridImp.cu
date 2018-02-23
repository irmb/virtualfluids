#include "GridImp.cuh"

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
#include "GridInterface.cuh"

#include "geometries/Object.h"


CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

HOST GridImp::GridImp(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution distribution)
    : startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta), field(nullptr), distribution(distribution),
    gridInterface(nullptr), neighborIndexX(nullptr), neighborIndexY(nullptr), neighborIndexZ(nullptr), matrixIndex(nullptr), gridStrategy(gridStrategy)
{
    initalNumberOfNodesAndSize();
}

HOST GridImp::GridImp(Object* object, real delta, SPtr<GridStrategy> gridStrategy, Distribution distribution) : object(object), delta(delta), gridStrategy(gridStrategy), distribution(distribution)
{
    initalBoundingBoXStartValues();
    initalNumberOfNodesAndSize();
}

HOST void GridImp::initalBoundingBoXStartValues()
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
    const real length = endX - startX;
    const real width = endY - startY;
    const real height = endZ - startZ;

    nx = int((length + delta) / delta) + 1; // +1 stopper node
    ny = int((width + delta) / delta) + 1; // +1 stopper node
    nz = int((height + delta) / delta) + 1; // +1 stopper node

    this->size = nx * ny * nz;
    this->reducedSize = size;
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

void GridImp::allocateGridMemory()
{
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

HOSTDEVICE uint GridImp::getReducedSize() const
{
    return this->reducedSize;
}

HOSTDEVICE void GridImp::findInnerNode(uint index)
{
    this->matrixIndex[index] = index;

    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    if(object->isPointInObject(x, y, z, 0, 0))
    {
        this->setFieldEntryToFluid(index);
        this->setNeighborIndices(index);
    } else
        this->setFieldEntryToOutOfGrid(index);
}


HOSTDEVICE void GridImp::findStopperNode(uint index)
{
    if(isEndOfGridStopper(index))
        this->setFieldEntryToSolid(index);
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

HOST void GridImp::removeOverlapNodes(SPtr<Grid> finerGrid)
{
    gridStrategy->createGridInterface(shared_from_this(), std::static_pointer_cast<GridImp>(finerGrid));
}

HOSTDEVICE void GridImp::createGridInterface(uint index, const GridImp& finerGrid)
{
    gridInterface->findInterface(index, this, &finerGrid);
    //this->findGridInterface(index, finerGrid);
    //this->setOverlapNodeToInvalid(index, finerGrid);
}

HOSTDEVICE void GridImp::findGridInterface(uint index, const GridImp& finerGrid)
{
    gridInterface->findCF(index, this, &finerGrid);
    gridInterface->findFC(index, this, &finerGrid);
}

HOSTDEVICE void GridImp::setOverlapNodeToInvalid(uint index, const GridImp& finerGrid)
{
    if (this->isInside(index, finerGrid))
        this->setFieldEntryToInvalid(index);
}

HOSTDEVICE bool GridImp::isInside(uint index, const GridImp& finerGrid) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    const real overlapWithStopper = 3 * this->delta;
    const real overlap = 2 * this->delta;

    return finerGrid.isInside(x, y, z, overlapWithStopper, overlap);
}

HOSTDEVICE bool GridImp::isInside(const real& x, const real& y, const real& z, const real& minOffset, const real& maxOffset) const
{
    uint index1 = transCoordToIndex(x - minOffset, y - minOffset, z - minOffset);
    uint index2 = transCoordToIndex(x + minOffset, y + minOffset, z + maxOffset);
    return isFluid(index1) && isFluid(index2);

    //return this->object->isPointInObject(x, y, z, minOffset, maxOffset);
}

HOSTDEVICE bool GridImp::isOnInterface(const real& x, const real& y, const real& z, const real& minOffset, const real& maxOffset) const
{
    return this->object->isOnBoundary(x, y, z, minOffset, maxOffset);
}


HOST void GridImp::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ)
{
    this->periodicityX = periodicityX;
    this->periodicityY = periodicityY;
    this->periodicityZ = periodicityZ;
}

HOSTDEVICE bool GridImp::is(uint index, char type) const
{
    return field[index] == type;
}

HOSTDEVICE bool GridImp::isFluid(uint index) const
{
    return field[index] == FLUID;
}

HOSTDEVICE bool GridImp::isSolid(uint index) const
{
    return field[index] == SOLID;
}

HOSTDEVICE bool GridImp::isOutOfGrid(uint index) const
{
    return field[index] == OUT_OF_GRID;
}

HOSTDEVICE bool GridImp::isInvalid(uint index) const
{
    return field[index] == INVALID_NODE;
}

HOSTDEVICE bool GridImp::isQ(uint index) const
{
    return field[index] == Q;
}

HOSTDEVICE bool GridImp::isRb(uint index) const
{
    return field[index] == VELOCITY || field[index] == PRESSURE || field[index] == NOSLIP || field[index] == SOLID;
}

HOSTDEVICE void GridImp::setFieldEntryToFluid(uint index)
{
    this->field[index] = FLUID;
}

HOSTDEVICE void GridImp::setFieldEntryToSolid(uint index)
{
    this->field[index] = SOLID;
}

HOSTDEVICE void GridImp::setFieldEntryToInvalid(uint index)
{
    this->field[index] = INVALID_NODE;
}

HOSTDEVICE void GridImp::setFieldEntryToOutOfGrid(uint index)
{
    this->field[index] = OUT_OF_GRID;
}

HOSTDEVICE void GridImp::setFieldEntry(const Vertex &v, char val)
{
    this->field[transCoordToIndex(v)] = val;
}

HOSTDEVICE char GridImp::getFieldEntry(const Vertex &v) const
{
    return this->field[transCoordToIndex(v)];
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


HOSTDEVICE int GridImp::getIndex(uint matrixIndex) const
{
    return this->matrixIndex[matrixIndex];
}

HOSTDEVICE char GridImp::getFieldEntry(uint index) const
{
    return this->field[index];
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
    if (isSolid(index))
        return;

    if (field[index] != INVALID_NODE && isNeighborInvalid(index))
    {
        field[index] = INVALID_NODE;
        invalidNodeFound = true;
    }
}

HOSTDEVICE bool GridImp::isNeighborInvalid(const int &index)
{
    return (field[neighborIndexX[index]] == INVALID_NODE || field[neighborIndexY[index]] == INVALID_NODE || field[
        neighborIndexZ[index]] == INVALID_NODE);
}

HOSTDEVICE void GridImp::findNeighborIndex(int index)
{
    if (this->isOverlapStopper(index) || isEndOfGridStopper(index))
    {
        this->setFieldEntryToSolid(index);
        this->setStopperNeighborCoords(index);
        return;
    }

    if (this->matrixIndex[index] == -1)
    {
        this->neighborIndexX[index] = -1;
        this->neighborIndexY[index] = -1;
        this->neighborIndexZ[index] = -1;
        return;
    }

    real x, y, z;
    this->transIndexToCoords(index, x, y, z);
    real neighborXCoord, neighborYCoord, neighborZCoord;
    getNeighborCoords(neighborXCoord, neighborYCoord, neighborZCoord, x, y, z);
    this->neighborIndexX[index] = getNeighborIndex(neighborXCoord, y, z);
    this->neighborIndexY[index] = getNeighborIndex(x, neighborYCoord, z);
    this->neighborIndexZ[index] = getNeighborIndex(x, y, neighborZCoord);
}

HOSTDEVICE bool GridImp::isEndOfGridStopper(uint index) const
{
    return this->isFluid(index) && nodeInNextCellIs(index, OUT_OF_GRID);
}

//bool GridImp::isEndOfGridStopper(uint index) const
//{
//    real x, y, z;
//    this->transIndexToCoords(index, x, y, z);
//    return (x > this->endX || y > this->endY || z > this->endZ);
//}


HOSTDEVICE bool GridImp::isOverlapStopper(uint index) const
{
    return this->isFluid(index) && nodeInNextCellIs(index, INVALID_NODE);
}

HOSTDEVICE bool GridImp::nodeInNextCellIs(int index, char type) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    real neighborX, neighborY, neighborZ;
    this->getNeighborCoords(neighborX, neighborY, neighborZ, x, y, z);

    const uint indexX = transCoordToIndex(neighborX, y, z);
    const uint indexY = transCoordToIndex(x, neighborY, z);
    const uint indexZ = transCoordToIndex(x, y, neighborZ);
     
    const uint indexXY = transCoordToIndex(neighborX, neighborY, z);
    const uint indexYZ = transCoordToIndex(x, neighborY, neighborZ);
    const uint indexXZ = transCoordToIndex(neighborX, y, neighborZ);
     
    const uint indexXYZ = transCoordToIndex(neighborX, neighborY, neighborZ);

    const bool isInvalidNeighborX = this->is(indexX, type);
    const bool isInvalidNeighborY = this->is(indexY, type);
    const bool isInvalidNeighborXY  = this->is(indexXY, type);
    const bool isInvalidNeighborZ   = this->is(indexZ, type);
    const bool isInvalidNeighborYZ = this->is(indexYZ, type);
    const bool isInvalidNeighborXZ = this->is(indexXZ, type);
    const bool isInvalidNeighborXYZ = this->is(indexXYZ, type);

    return isInvalidNeighborX || isInvalidNeighborY || isInvalidNeighborXY || isInvalidNeighborZ || isInvalidNeighborYZ
        || isInvalidNeighborXZ || isInvalidNeighborXYZ;
}

HOSTDEVICE int GridImp::getNeighborIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const
{
    const int neighborIndex = transCoordToIndex(expectedX, expectedY, expectedZ);
    return matrixIndex[neighborIndex];
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
    const uint newIndex = matrixIndex[oldIndex];
    gridInterface->cf.coarse[index] = newIndex;
}

HOSTDEVICE void GridImp::findForGridInterfaceNewIndexFC(uint index)
{
    const uint oldIndex = gridInterface->fc.coarse[index];
    const uint newIndex = matrixIndex[oldIndex];
    gridInterface->fc.coarse[index] = newIndex;
}


HOST void GridImp::removeInvalidNodes()
{
    int removedNodes = 0;
    int newIndex = 0;
    for (uint index = 0; index < size; index++)
    {
        if (this->isInvalid(index) || this->isOutOfGrid(index))
        {
            matrixIndex[index] = -1;
            removedNodes++;
        } else
        {
            matrixIndex[index] = newIndex;
            newIndex++;
        }
    }
    reducedSize = size - removedNodes;
    printf("new size nodes: %d , delete nodes: %d\n", reducedSize, removedNodes);
}

HOSTDEVICE void GridImp::getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const
{
    neighborX = getNeighhborCoord(periodicityX, x, startX, endX);
    neighborY = getNeighhborCoord(periodicityY, y, startY, endY);
    neighborZ = getNeighhborCoord(periodicityZ, z, startZ, endZ);
}

HOSTDEVICE real GridImp::getNeighhborCoord(bool periodicity, real actualCoord, real startCoord, real endCoord) const
{
    if (periodicity)
        return CudaMath::lessEqual(actualCoord + delta, endCoord) ? actualCoord + delta : startCoord;
    else
        return actualCoord + delta;
}

HOSTDEVICE void GridImp::setStopperNeighborCoords(int index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    if (CudaMath::lessEqual(x + delta, endX + delta))
        neighborIndexX[index] = getNeighborIndex(x + delta, y, z);

    if (CudaMath::lessEqual(y + delta, endY + delta))
        neighborIndexY[index] = getNeighborIndex(x, y + delta, z);

    if (CudaMath::lessEqual(z + delta, endZ + delta))
        neighborIndexZ[index] = getNeighborIndex(x, y, z + delta);
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

void GridImp::setFieldEntry(uint index, char entry)
{
    this->field[index] = entry;
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

HOST void GridImp::setNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *neighborX,
                                 unsigned int *neighborY, unsigned int *neighborZ, unsigned int *geo) const
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
        if (this->matrixIndex[i] == -1)
            continue;

        real x, y, z;
        this->transIndexToCoords(i, x, y, z);

        xCoords[nodeNumber + 1] = x;
        yCoords[nodeNumber + 1] = y;
        zCoords[nodeNumber + 1] = z;
        neighborX[nodeNumber + 1] = uint(this->neighborIndexX[i] + 1);
        neighborY[nodeNumber + 1] = uint(this->neighborIndexY[i] + 1);
        neighborZ[nodeNumber + 1] = uint(this->neighborIndexZ[i] + 1);
        geo[nodeNumber + 1] = uint(this->isSolid(i) ? GEOSOLID : GEOFLUID);
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
        "), size: " << reducedSize << ", delta: " << delta << "\n";
    return oss.str();
}


