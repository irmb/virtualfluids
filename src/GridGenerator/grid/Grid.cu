#include "Grid.cuh"

#include "GridGenerator/global.h"
#include <stdio.h>
#include <time.h>

#include <sstream>
#include <algorithm>


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


CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

HOST Grid::Grid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, std::shared_ptr<GridStrategy> gridStrategy, Distribution distribution)
    : startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ), delta(delta), field(nullptr), distribution(distribution),
    gridInterface(nullptr), neighborIndexX(nullptr), neighborIndexY(nullptr), neighborIndexZ(nullptr), matrixIndex(nullptr), gridStrategy(gridStrategy)
{
    const real length = endX - startX;
    const real width = endY - startY;
    const real height = endZ - startZ;

    nx = int((length + delta) / delta) + 1; // +1 stopper node
    ny = int((width + delta) / delta) + 1; // +1 stopper node
    nz = int((height + delta) / delta) + 1; // +1 stopper node

    this->size = nx * ny * nz;
    this->reducedSize = size;
}

HOST SPtr<Grid> Grid::makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta,
                                 std::shared_ptr<GridStrategy> gridStrategy, Distribution d)
{
    SPtr<Grid> grid(new Grid(startX, startY, startZ, endX, endY, endZ, delta, gridStrategy, d));

    gridStrategy->allocateGridMemory(grid);
    *logging::out << logging::Logger::LOW << "-------------------------------------------\n";
    *logging::out << logging::Logger::LOW << "Initial field with fluid. \n";
    *logging::out << logging::Logger::LOW << "-------------------------------------------\n";
    time_t begin = clock();
    gridStrategy->initalNodes(grid);

    time_t end = clock();
    real time = (real)(real(end - begin) / CLOCKS_PER_SEC);
    *logging::out << logging::Logger::INTERMEDIATE << "Time initial field: " + SSTR(time / 1000) + "sec\n";
    *logging::out << logging::Logger::INTERMEDIATE << "-------------------------------------------\n";

    return grid;
}

HOST Grid::Grid()
{
    //printf("Constructor\n");
    //this->print();
};

HOST Grid::~Grid()
{
    //printf("Destructor\n");
    //this->print();
};

HOST void Grid::mesh(Geometry &geometry)
{
    clock_t begin = clock();

    gridStrategy->mesh(shared_from_this(), geometry);

    clock_t end = clock();
    real time = (real)(real(end - begin) / CLOCKS_PER_SEC);

    *logging::out << logging::Logger::INTERMEDIATE << "time grid generation: " + SSTR(time) + "s\n";
}

HOST void Grid::freeMemory()
{
    gridStrategy->freeMemory(shared_from_this());
}

HOST void Grid::removeOverlapNodes(SPtr<Grid> finerGrid)
{
    gridStrategy->removeOverlapNodes(shared_from_this(), finerGrid);
}

HOSTDEVICE void Grid::setOverlapNodeToInvalid(uint index, const Grid& finerGrid)
{
    gridInterface->findCF(index, this, &finerGrid);
    gridInterface->findFC(index, this, &finerGrid);

    if (this->isInside(index, finerGrid))
        this->setFieldEntryToInvalid(index);
}

HOSTDEVICE bool Grid::isInside(uint index, const Grid& finerGrid)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    const real overlapWithStopper = 3 * this->delta;
    const real overlap = 2 * this->delta;

    return 
        (x > finerGrid.startX + overlapWithStopper && x < finerGrid.endX - overlap) &&
        (y > finerGrid.startY + overlapWithStopper && y < finerGrid.endY - overlap) &&
        (z > finerGrid.startZ + overlapWithStopper && z < finerGrid.endZ - overlap);
}

HOST void Grid::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ)
{
    this->periodicityX = periodicityX;
    this->periodicityY = periodicityY;
    this->periodicityZ = periodicityZ;
}

HOSTDEVICE bool Grid::isFluid(uint index) const
{
    return field[index] == FLUID;
}

HOSTDEVICE bool Grid::isSolid(uint index) const
{
    return field[index] == SOLID;
}

HOSTDEVICE bool Grid::isInvalid(uint index) const
{
    return field[index] == INVALID_NODE;
}

HOSTDEVICE bool Grid::isQ(uint index) const
{
    return field[index] == Q;
}

HOSTDEVICE bool Grid::isRb(uint index) const
{
    return field[index] == VELOCITY || field[index] == PRESSURE || field[index] == NOSLIP || field[index] == SOLID;
}

HOSTDEVICE void Grid::setFieldEntryToFluid(uint index)
{
    this->field[index] = FLUID;
}

HOSTDEVICE void Grid::setFieldEntryToSolid(uint index)
{
    this->field[index] = SOLID;
}

HOST void Grid::setFieldEntryToInvalid(uint index)
{
    this->field[index] = INVALID_NODE;
}

HOSTDEVICE void Grid::setFieldEntry(const Vertex &v, char val)
{
    this->field[transCoordToIndex(v)] = val;
}

HOSTDEVICE char Grid::getFieldEntry(const Vertex &v) const
{
    return this->field[transCoordToIndex(v)];
}

HOSTDEVICE int Grid::transCoordToIndex(const real &x, const real &y, const real &z) const
{
    return transCoordToIndex(Vertex(x,y,z));
}

HOSTDEVICE int Grid::transCoordToIndex(const Vertex &v) const
{
    const int x = int((v.x - startX) / delta);
    const int y = int((v.y - startY) / delta);
    const int z = int((v.z - startZ) / delta);

#ifdef DEBUG
    if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz)
    {
        printf(
            "Function: transCoordToIndex. Coordinates are out of range and cannot calculate the index. Exit Program!\n");
        /* exit(1);*/
    }
#endif

    return x + nx * (y + ny * z);
}

HOSTDEVICE void Grid::transIndexToCoords(const int index, real &x, real &y, real &z) const
{
#ifdef DEBUG
    if (index < 0 || index >= (int)size)
    {
        printf("Function: transIndexToCoords. Grid Index: %d, size: %d. Exit Program!\n", index, size); /*exit(1);*/ 
    }
#endif
    x = index % nx;
    y = (index / nx) % ny;
    z = ((index / nx) / ny) % nz;

    x = (x * delta) + startX;
    y = (y * delta) + startY;
    z = (z * delta) + startZ;
}


HOSTDEVICE void Grid::setDebugPoint(const Vertex &point, const int pointValue)
{
    if (getFieldEntry(point) == INVALID_NODE && pointValue == SOLID)
        setFieldEntry(point, pointValue);

    if (getFieldEntry(point) != SOLID && getFieldEntry(point) != Q && getFieldEntry(point) != INVALID_NODE && pointValue
        != 3 && pointValue != 2)
        setFieldEntry(point, pointValue);
}

HOSTDEVICE bool Grid::isOutOfRange(const Vertex &v) const
{
    return v.x < startX || v.y < startY || v.z < startZ || v.x > endX || v.y > endY || v.z > endZ;
}

HOSTDEVICE void Grid::meshTriangle(const Triangle &triangle)
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

HOSTDEVICE void Grid::calculateQs(const Vertex &point, const Triangle &triangle)
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

HOSTDEVICE void Grid::setNeighborIndices(const int &index)
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    real neighborX, neighborY, neighborZ;
    this->getNeighborCoords(neighborX, neighborY, neighborZ, x, y, z);

    neighborIndexX[index] = uint(transCoordToIndex(neighborX, y, z));
    neighborIndexY[index] = uint(transCoordToIndex(x, neighborY, z));
    neighborIndexZ[index] = uint(transCoordToIndex(x, y, neighborZ));

    //if (grid.isRB(index)) {
    //if (neighborX == 0) neighborIndexX[index] = 0;
    //if (neighborY == 0) neighborIndexY[index] = 0;
    //if (neighborZ == 0) neighborIndexZ[index] = 0;
    //}
}

HOSTDEVICE void Grid::setInvalidNode(const int &index, bool &invalidNodeFound)
{
    if (isSolid(index))
        return;

    if (field[index] != INVALID_NODE && isNeighborInvalid(index))
    {
        field[index] = INVALID_NODE;
        invalidNodeFound = true;
    }
}

HOSTDEVICE bool Grid::isNeighborInvalid(const int &index)
{
    return (field[neighborIndexX[index]] == INVALID_NODE || field[neighborIndexY[index]] == INVALID_NODE || field[
        neighborIndexZ[index]] == INVALID_NODE);
}

HOSTDEVICE void Grid::findNeighborIndex(int index)
{
    if (this->matrixIndex[index] == -1)
    {
        this->neighborIndexX[index] = -1;
        this->neighborIndexY[index] = -1;
        this->neighborIndexZ[index] = -1;
        return;
    }

    if(this->isOverlapStopper(index) || isEndOfGridStopper(index))
    {
        this->setFieldEntryToSolid(index);
        this->setStopperNeighborCoords(index);
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

HOSTDEVICE bool Grid::isOverlapStopper(uint index) const
{
    return this->isFluid(index) && nodeInNextCellIsInvalid(index);
}


HOSTDEVICE bool Grid::nodeInNextCellIsInvalid(int index) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    real neighborX, neighborY, neighborZ;
    this->getNeighborCoords(neighborX, neighborY, neighborZ, x, y, z);

    const uint indexX = (uint)transCoordToIndex(neighborX, y, z);
    const uint indexY = (uint)transCoordToIndex(x, neighborY, z);
    const uint indexZ = (uint)transCoordToIndex(x, y, neighborZ);
     
    const uint indexXY = (uint)transCoordToIndex(neighborX, neighborY, z);
    const uint indexYZ = (uint)transCoordToIndex(x, neighborY, neighborZ);
    const uint indexXZ = (uint)transCoordToIndex(neighborX, y, neighborZ);
     
    const uint indexXYZ = (uint)transCoordToIndex(neighborX, neighborY, neighborZ);

    const bool isInvalidNeighborX = this->isInvalid(indexX);
    const bool isInvalidNeighborY = this->isInvalid(indexY);
    const bool isInvalidNeighborXY  = this->isInvalid(indexXY);
    const bool isInvalidNeighborZ   = this->isInvalid(indexZ);
    const bool isInvalidNeighborYZ  = this->isInvalid(indexYZ);
    const bool isInvalidNeighborXZ  = this->isInvalid(indexXZ);
    const bool isInvalidNeighborXYZ = this->isInvalid(indexXYZ);

    return isInvalidNeighborX || isInvalidNeighborY || isInvalidNeighborXY || isInvalidNeighborZ || isInvalidNeighborYZ
        || isInvalidNeighborXZ || isInvalidNeighborXYZ;
}

HOSTDEVICE int Grid::getNeighborIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const
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

HOST void Grid::findForGridInterfaceNewIndexCF(uint index)
{
    const uint oldIndex = gridInterface->cf.coarse[index];
    const uint newIndex = matrixIndex[oldIndex];
    gridInterface->cf.coarse[index] = newIndex;
}

HOST void Grid::findForGridInterfaceNewIndexFC(uint index)
{
    const uint oldIndex = gridInterface->fc.coarse[index];
    const uint newIndex = matrixIndex[oldIndex];
    gridInterface->fc.coarse[index] = newIndex;
}


HOST void Grid::removeInvalidNodes()
{
    int removedNodes = 0;
    int newIndex = 0;
    for (uint index = 0; index < size; index++)
    {
        if (this->isInvalid(index))
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
    printf("new size coords: %d , delete nodes: %d\n", reducedSize, removedNodes);
}

HOSTDEVICE void Grid::getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const
{
    neighborX = getNeighhborCoord(periodicityX, x, startX, endX);
    neighborY = getNeighhborCoord(periodicityY, y, startY, endY);
    neighborZ = getNeighhborCoord(periodicityZ, z, startZ, endZ);
}

HOSTDEVICE real Grid::getNeighhborCoord(bool periodicity, real actualCoord, real startCoord, real endCoord) const
{
    if (periodicity)
        return CudaMath::lessEqual(actualCoord + delta, endCoord) ? actualCoord + delta : startCoord;
    else
        return actualCoord + delta;
}

void Grid::setStopperNeighborCoords(int index)
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


HOST bool Grid::isEndOfGridStopper(uint index) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);
    return (x > this->endX || y > this->endY || z > this->endZ);
}


//HOSTDEVICE  bool Grid::isStopper(int index) const
//{
//    return isSolid(index) && previousCellHasFluid(index);
//}
//
//HOSTDEVICE  bool Grid::previousCellHasFluid(int index) const
//{
//    real x, y, z;
//    this->transIndexToCoords(index, x, y, z);
//
//    real previousX = x - delta >= 0 ? x - delta : this->endX;
//    real previousY = y - delta >= 0 ? y - delta : this->endY;
//    real previousZ = z - delta >= 0 ? z - delta : this->endZ;
//
//    int indexpreviousX   = this->transCoordToIndex(previousX, y, z);
//    int indexpreviousY   = this->transCoordToIndex(x, previousY, z);
//    int indexpreviousXY  = this->transCoordToIndex(previousX, previousY, z);
//    int indexpreviousZ   = this->transCoordToIndex(x, y, previousZ);
//    int indexpreviousZX  = this->transCoordToIndex(previousX, y, previousZ);
//    int indexpreviousZY  = this->transCoordToIndex(x, previousY, previousZ);
//    int indexpreviousZYX = this->transCoordToIndex(previousX, previousY, previousZ);
//
//    return (!this->isSolid(indexpreviousX) || !this->isSolid(indexpreviousY) || !this->isSolid(indexpreviousZ)
//        || !this->isSolid(indexpreviousZYX) || !this->isSolid(indexpreviousXY) || !this->isSolid(indexpreviousZY) || !
//        this->isSolid(indexpreviousZX));
//}


uint Grid::getNumberOfNodes() const
{
    return this->reducedSize;
}

uint Grid::getNumberOfNodesCF() const
{
    return this->gridInterface->cf.numberOfEntries;
}

uint Grid::getNumberOfNodesFC() const
{
    return this->gridInterface->fc.numberOfEntries;
}

void Grid::setCFC(uint* iCellCfc) const
{
    setGridInterface(iCellCfc, this->gridInterface->cf.coarse, this->gridInterface->cf.numberOfEntries);
}

void Grid::setCFF(uint* iCellCff) const
{
    setGridInterface(iCellCff, this->gridInterface->cf.fine, this->gridInterface->cf.numberOfEntries);
}

void Grid::setFCC(uint* iCellFcc) const
{
    setGridInterface(iCellFcc, this->gridInterface->fc.coarse, this->gridInterface->fc.numberOfEntries);
}

void Grid::setFCF(uint* iCellFcf) const
{
    setGridInterface(iCellFcf, this->gridInterface->fc.fine, this->gridInterface->fc.numberOfEntries);
}

void Grid::setGridInterface(uint* gridInterfaceList, const uint* oldGridInterfaceList, uint size)
{
    for (uint i = 0; i < size; i++)
        gridInterfaceList[i] = oldGridInterfaceList[i] + 1;
}



HOSTDEVICE void Grid::print() const
{
    printf("min: (%2.2f, %2.2f, %2.2f), max: (%2.2f, %2.2f, %2.2f), size: %d, delta: %2.2f\n", startX, startY, startZ,
           endX, endY, endZ, size, delta);
}

std::string Grid::toString() const
{
    std::ostringstream oss;
    oss << 
        "min: (" << startX << ", " << startY << ", " << startZ << 
        "), max: " << endX << ", " << endY << ", " << endZ <<
        "), size: " << reducedSize << ", delta: " << delta << "\n";
    return oss.str();
}


