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

HOSTDEVICE int Grid::transCoordToIndex(const int &x, const int &y, const int &z) const
{
	return transCoordToIndex(Vertex((doubflo)x, (doubflo)y, (doubflo)z));
}

HOSTDEVICE int Grid::transCoordToIndex(const Vertex &v) const
{
#ifdef DEBUG
	if (isOutOfRange(v))
	{ printf("Function: transCoordToIndex. Coordinates are out of range and cannot calculate the index. Exit Program!\n");/* exit(1);*/ };
#endif
	return (int)(v.x + nx * (v.y + ny * v.z));
}

HOSTDEVICE void Grid::transIndexToCoords(const int index, unsigned int &x, unsigned int &y, unsigned int &z) const
{
#ifdef DEBUG
	if (index < 0 || index >= (int)size)
	{
        printf("Function: transIndexToCoords. Grid Index: %d, size: %d. Exit Program!\n", index, size); /*exit(1);*/ 
    };
#endif
    x = index % nx;
    y = (index / nx) % ny;
    z = ((index / nx) / ny) % nz;
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


HOSTDEVICE void Grid::meshTriangle(const Triangle &triangle)
{
	int x, y, z;
	Vertex point;

	BoundingBox<int> box = BoundingBox<int>::makeNodeBox(triangle);

	for (x = box.minX; x <= box.maxX; x++) {
		for (y = box.minY; y <= box.maxY; y++) {
			for (z = box.minZ; z <= box.maxZ; z++) {
				point = Vertex((doubflo)x, (doubflo)y, (doubflo)z);
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
	doubflo qVal;
	int err;
	for (int i = d.dir_start; i <= d.dir_end; i++) {
	#if defined(__CUDA_ARCH__)
		dir = Vertex(DIRECTIONS[i][0], DIRECTIONS[i][1], DIRECTIONS[i][2]);
	#else
		dir = Vertex((doubflo)d.dirs[i * DIMENSION + 0], (doubflo)d.dirs[i * DIMENSION + 1], (doubflo)d.dirs[i * DIMENSION + 2]);
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
    unsigned int x, y, z;
    this->transIndexToCoords(index, x, y, z);

	unsigned int neighborX, neighborY, neighborZ;
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
    unsigned int x, y, z;
    int nodeIndex = this->matrixIndex[index];
    this->transIndexToCoords(this->matrixIndex[index], x, y, z);

    unsigned int neighborXCoord, neighborYCoord, neighborZCoord;
    getNeighborCoords(neighborXCoord, neighborYCoord, neighborZCoord, x, y, z);
    this->neighborIndexX[nodeIndex] = getNeighborIndex(index, this->neighborIndexX[nodeIndex], neighborXCoord, y, z);
    this->neighborIndexY[nodeIndex] = getNeighborIndex(index, this->neighborIndexY[nodeIndex], x, neighborYCoord, z);
    this->neighborIndexZ[nodeIndex] = getNeighborIndex(index, this->neighborIndexZ[nodeIndex], x, y, neighborZCoord);
}

HOSTDEVICE int Grid::getNeighborIndex(const int &nodeIndex, int &neighborIndex, const int &expectedX, const int &expectedY, const int &expectedZ)
{
    while (neighborIndex >= (int)this->reducedSize)
    {
        neighborIndex--;
        //printf("here\n");
    }
    if (neighborIndex >= 0) {
        unsigned int neighborX, neighborY, neighborZ;
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

HOSTDEVICE void Grid::getNeighborCoords(unsigned int &neighborX, unsigned int &neighborY, unsigned int &neighborZ, const unsigned int x, const unsigned int y, const unsigned int z) const
{
	neighborX = x + 1 < nx ? x + 1 : 0;
	neighborY = y + 1 < ny ? y + 1 : 0;
	neighborZ = z + 1 < nz ? z + 1 : 0;
}


HOSTDEVICE bool Grid::isStopper(int index) const
{
    return isSolid(index) && previousCellHasFluid(index);
}

HOSTDEVICE bool Grid::previousCellHasFluid(int index) const
{
    unsigned int ux, uy, uz;
    this->transIndexToCoords(index, ux, uy, uz);

    int x(ux), y(uy), z(uz);
    int previousX = x - 1 >= 0 ? x - 1 : this->nx - 1;
    int previousY = y - 1 >= 0 ? y - 1 : this->ny - 1;
    int previousZ = z - 1 >= 0 ? z - 1 : this->nz - 1;

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
