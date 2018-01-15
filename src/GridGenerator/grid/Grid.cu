#include "Grid.cuh"

#include "GridGenerator/global.h"
#include <stdio.h>
#include <sstream>
#include <algorithm>

#include <GridGenerator/utilities/math/CudaMath.cuh>

#include <GridGenerator/geometries/Vertex/Vertex.cuh>
#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>

#include <GridGenerator/grid/NodeValues.h>
#include <GridGenerator/grid/distributions/Distribution.h>

__constant__ int DIRECTIONS[DIR_END_MAX][DIMENSION];

HOSTDEVICE Grid::Grid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, Distribution &d) : sX(startX), sY(startY), sZ(startZ), eX(endX), eY(endY), eZ(endZ), delta(delta), d(d)
{
    real length = endX - startX;
    real width = endY - startY;
    real height = endZ - startZ;

    nx = (int)((length + delta) / delta);
    ny = (int)((width  + delta) / delta);
    nz = (int)((height + delta) / delta);

    this->size = nx * ny * nz;
    this->reducedSize = size;

    const unsigned long distributionSize = size * (d.dir_end + 1);
    this->d.f = new real[distributionSize](); // automatic initialized with zeros
    this->field = new char[this->size];

    this->neighborIndexX = new int[this->size];
    this->neighborIndexY = new int[this->size];
    this->neighborIndexZ = new int[this->size];

    this->matrixIndex = new uint[this->size];
}

HOSTDEVICE Grid::Grid(){};
HOSTDEVICE Grid::Grid(char *field, int startX, int startY, int startZ, int nx, int ny, int nz, Distribution &d)
    : field(field), startX(startX), startY(startY), startZ(startZ), nx(nx), ny(ny), nz(nz), d(d)
{
    this->size = nx * ny * nz;
    this->reducedSize = size;
}

HOSTDEVICE bool Grid::isFluid(int index) const
{
    return field[index] == FLUID;
}

HOSTDEVICE bool Grid::isSolid(int index) const 
{
    return field[index] == SOLID;
}

HOSTDEVICE bool Grid::isQ(int index) const
{
    return field[index] == Q;
}

HOSTDEVICE  bool Grid::isRb(int index) const
{
    return field[index] == VELOCITY || field[index] == PRESSURE || field[index] == NOSLIP || field[index] == SOLID;
}

HOSTDEVICE void Grid::setFieldEntryToFluid(unsigned int index)
{
	this->field[index] = FLUID;
}

HOSTDEVICE void Grid::setFieldEntryToSolid(unsigned int index)
{
	this->field[index] = SOLID;
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
//#ifdef DEBUG
//	if (isOutOfRange(v))
//	{ printf("Function: transCoordToIndex. Coordinates are out of range and cannot calculate the index. Exit Program!\n");/* exit(1);*/ };
//#endif
    int x = (int)((v.x - sX) / delta);
    int y = (int)((v.y - sY) / delta);
    int z = (int)((v.z - sZ) / delta);

	return x + nx * (y + ny * z);
}

HOSTDEVICE void Grid::transIndexToCoords(const int index, real &x, real &y, real &z) const
{
//#ifdef DEBUG
//	if (index < 0 || index >= (int)size)
//	{
//        printf("Function: transIndexToCoords. Grid Index: %d, size: %d. Exit Program!\n", index, size); /*exit(1);*/ 
//    };
//#endif
    x = index % nx;
    y = (index / nx) % ny;
    z = ((index / nx) / ny) % nz;

    x = (x + sX) * delta;
    y = (y + sY) * delta;
    z = (z + sZ) * delta;
}

char* Grid::toString(const char* name) const
{
    std::stringstream ss;
    ss << "\n" << name << " " << nx << " " << ny << " " << nz;
    return strdup(ss.str().c_str());
}


HOSTDEVICE void Grid::print() const
{
    printf("Dimension: (%d, %d, %d), size: %d, offset: (%d, %d, %d)\n", nx, ny, nz, size, startX, startY, startZ);
}

HOSTDEVICE void Grid::setDebugPoint(const Vertex &point, const int pointValue)
{
    if (getFieldEntry(point) == INVALID_NODE && pointValue == SOLID)
        setFieldEntry(point, pointValue);

	if (getFieldEntry(point) != SOLID && getFieldEntry(point) != Q && getFieldEntry(point) != INVALID_NODE && pointValue != 3 && pointValue != 2)
		setFieldEntry(point, pointValue);
}

HOSTDEVICE bool Grid::isOutOfRange(const Vertex &v) const
{
	int gridX = (int)v.x - startX;
	int gridY = (int)v.y - startY;
	int gridZ = (int)v.z - startZ;
	return (gridX < 0 || gridY < 0 || gridZ < 0 || gridX >= (int)nx || gridY >= (int)ny || gridZ >= (int)nz);
}


HOSTDEVICE void Grid::meshTriangleExact(const Triangle &triangle)
{
    BoundingBox<real> box = BoundingBox<real>::makeRealNodeBox(triangle, delta);

    for (real x = box.minX; x <= box.maxX; x += delta) {
        for (real y = box.minY; y <= box.maxY; y += delta) {
            for (real z = box.minZ; z <= box.maxZ; z += delta) {
                Vertex point(x, y, z);
                //if (isOutOfRange(point))
                //    continue;
                const int value = triangle.isUnderFace(point);
                setDebugPoint(point, value);

                //if (value == Q)
                //    calculateQs(point, triangle);
            }
        }
    }
}

HOSTDEVICE void Grid::meshTriangle(const Triangle &triangle)
{
	int x, y, z;
	Vertex point;

	BoundingBox<int> box = BoundingBox<int>::makeNodeBox(triangle);

	for (x = box.minX; x <= box.maxX; x++) {
		for (y = box.minY; y <= box.maxY; y++) {
			for (z = box.minZ; z <= box.maxZ; z++) {
				point = Vertex((real)x, (real)y, (real)z);
				if (isOutOfRange(point))
					continue;
                int value = triangle.isUnderFace(point);
                setDebugPoint(point, value);

                if (value == Q)
                    calculateQs(point, triangle);
			}
		}
	}
}

HOSTDEVICE void Grid::calculateQs(const Vertex &point, const Triangle &triangle)
{
	Vertex pointOnTriangle, dir;
	//VertexInteger solid_node;
	real qVal;
	int err;
	for (int i = d.dir_start; i <= d.dir_end; i++) {
	#if defined(__CUDA_ARCH__)
		dir = Vertex(DIRECTIONS[i][0], DIRECTIONS[i][1], DIRECTIONS[i][2]);
	#else
		dir = Vertex((real)d.dirs[i * DIMENSION + 0], (real)d.dirs[i * DIMENSION + 1], (real)d.dirs[i * DIMENSION + 2]);
	#endif

		err = triangle.getTriangleIntersection(point, dir, pointOnTriangle, qVal);

		if (err != 0 && qVal <= 1.0f) {
			//solid_node = VertexInteger(actualPoint.x + dir.x, actualPoint.y + dir.y, actualPoint.z + dir.z);
			d.f[i*size + transCoordToIndex(point)] = qVal;
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

	neighborIndexX[index] = (unsigned int) transCoordToIndex(neighborX, y, z);
    neighborIndexY[index] = (unsigned int) transCoordToIndex(x, neighborY, z);
    neighborIndexZ[index] = (unsigned int) transCoordToIndex(x, y, neighborZ);

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
    return (field[neighborIndexX[index]] == INVALID_NODE || field[neighborIndexY[index]] == INVALID_NODE || field[neighborIndexZ[index]] == INVALID_NODE);
}

HOSTDEVICE void Grid::findNeighborIndex(int index)
{
    real x, y, z;
    int nodeIndex = this->matrixIndex[index];
    this->transIndexToCoords(this->matrixIndex[index], x, y, z);

    real neighborXCoord, neighborYCoord, neighborZCoord;
    getNeighborCoords(neighborXCoord, neighborYCoord, neighborZCoord, x, y, z);
    this->neighborIndexX[nodeIndex] = getNeighborIndex(index, this->neighborIndexX[nodeIndex], neighborXCoord, y, z);
    this->neighborIndexY[nodeIndex] = getNeighborIndex(index, this->neighborIndexY[nodeIndex], x, neighborYCoord, z);
    this->neighborIndexZ[nodeIndex] = getNeighborIndex(index, this->neighborIndexZ[nodeIndex], x, y, neighborZCoord);
}

HOSTDEVICE int Grid::getNeighborIndex(const int &nodeIndex, int &neighborIndex, const real &expectedX, const real &expectedY, const real &expectedZ)
{
    while (neighborIndex >= (int)this->reducedSize)
    {
        neighborIndex--;
        //printf("here\n");
    }
    if (neighborIndex >= 0) {
        real neighborX, neighborY, neighborZ;
        this->transIndexToCoords(this->matrixIndex[neighborIndex], neighborX, neighborY, neighborZ);
        while (!(neighborX == expectedX && neighborY == expectedY && neighborZ == expectedZ)) {
            neighborIndex--;
            //printf("expectedNeighborCoords:(%d, %d, %d), actualNeighborCoords:(%d,%d,%d), neighborIndex: %d\n", expectedX, expectedY, expectedZ,neighborX ,neighborY ,neighborZ, neighborIndex);
            this->transIndexToCoords(this->matrixIndex[neighborIndex], neighborX, neighborY, neighborZ);

            if (neighborIndex == nodeIndex) {
                neighborIndex = -1;
                break;
            }
        }
        return neighborIndex;
    }
    else {
        return -1;
    }
}


HOST void Grid::removeInvalidNodes()
{
    std::vector<unsigned int> stl_vector(size);
    stl_vector.assign(this->matrixIndex, this->matrixIndex + this->size);

    int oldsize = (int)stl_vector.size();
    printf("size coords: %d \n", oldsize);
    std::vector<unsigned int>::iterator end_vaild = std::remove_if(stl_vector.begin(), stl_vector.end(), [this](const unsigned int &index)
    {
        return this->field[index] == INVALID_NODE;
    });

    stl_vector.erase(end_vaild, stl_vector.end());
    printf("new size coords: %zd , delete nodes: %zd\n", stl_vector.size(), oldsize - stl_vector.size());

    unsigned int *indices_reduced = new  unsigned int[stl_vector.size()];
    for (size_t i = 0; i < stl_vector.size(); i++)
        indices_reduced[i] = stl_vector[i];
    
    this->reducedSize = (int)stl_vector.size();
    delete[]this->matrixIndex;
    this->matrixIndex = indices_reduced;
}

HOSTDEVICE void Grid::getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, const real x, const real y, const real z) const
{
	neighborX = x + delta < nx ? x + delta : startX;
	neighborY = y + delta < ny ? y + delta : startY;
	neighborZ = z + delta < nz ? z + delta : startZ;
}


HOSTDEVICE bool Grid::isStopper(int index) const
{
    return isSolid(index) && previousCellHasFluid(index);
}

HOSTDEVICE bool Grid::previousCellHasFluid(int index) const
{
    real x, y, z;
    this->transIndexToCoords(index, x, y, z);

    real previousX = x - delta >= 0 ? x - delta : this->eX;
    real previousY = y - delta >= 0 ? y - delta : this->eY;
    real previousZ = z - delta >= 0 ? z - delta : this->eZ;

    int indexpreviousX   = this->transCoordToIndex(previousX, y, z);
    int indexpreviousY   = this->transCoordToIndex(x, previousY, z);
    int indexpreviousXY  = this->transCoordToIndex(previousX, previousY, z);
    int indexpreviousZ   = this->transCoordToIndex(x, y, previousZ);
    int indexpreviousZX  = this->transCoordToIndex(previousX, y, previousZ);
    int indexpreviousZY  = this->transCoordToIndex(x, previousY, previousZ);
    int indexpreviousZYX = this->transCoordToIndex(previousX, previousY, previousZ);

    return (!this->isSolid(indexpreviousX) || !this->isSolid(indexpreviousY) || !this->isSolid(indexpreviousZ)
        || !this->isSolid(indexpreviousZYX) || !this->isSolid(indexpreviousXY) || !this->isSolid(indexpreviousZY) || !this->isSolid(indexpreviousZX));
}
