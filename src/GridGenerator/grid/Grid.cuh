#ifndef GRID_H
#define GRID_H

#include "GridGenerator/global.h"


#include <stdio.h>
#include <sstream>
#include "cuda.h"
#include "cuda_runtime.h"

#include <GridGenerator/grid/distributions/Distribution.h>

struct Vertex;
struct Triangle;

extern CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

struct VF_PUBLIC Grid 
{
	char *field;
	int startX, startY, startZ;
	unsigned int nx, ny, nz;
	unsigned int size;
	Distribution d;

    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;
    unsigned int *matrixIndex;
    int reducedSize;

	HOSTDEVICE Grid();
	HOSTDEVICE Grid(char *field, int startX, int startY, int startZ, int nx, int ny, int nz, Distribution &d);

	HOSTDEVICE bool isFluid(int index) const;
	HOSTDEVICE bool isSolid(int index) const;
	HOSTDEVICE bool isQ(int index) const;
	HOSTDEVICE bool isRb(int index) const;
	HOSTDEVICE void setFieldEntryToFluid(unsigned int index);
	HOSTDEVICE void setFieldEntryToSolid(unsigned int index);
	HOSTDEVICE void setFieldEntry(const Vertex &v, char val);
	HOSTDEVICE char getFieldEntry(const Vertex &v) const;
	HOSTDEVICE int transCoordToIndex(const int &x, const int &y, const int &z) const;
	HOSTDEVICE int transCoordToIndex(const Vertex &v) const;
	HOSTDEVICE void transIndexToCoords(const int index, unsigned int &x, unsigned int &y, unsigned int &z) const;
	HOSTDEVICE void print() const;
	HOSTDEVICE void setDebugPoint(const Vertex &actualPoint, const int pointValue);
	HOSTDEVICE bool isOutOfRange(const Vertex &actualPoint) const;

	/*---------------------------------------------------------------------------------*/
	HOSTDEVICE void meshTriangle(const Triangle &actualTriangle);

	HOSTDEVICE void calculateQs(const Vertex &point, const Triangle &actualTriangle);

	char* toString(const char* name) const;

    /*---------------------------------------------------------------------------------*/
    HOSTDEVICE void setNeighborIndices(const int &index);
	HOSTDEVICE void getNeighborCoords(unsigned int &neighborX, unsigned int &neighborY, unsigned int &neighborZ, const unsigned int x, const unsigned int y, const unsigned int z) const;
    HOSTDEVICE void findNeighborIndex(int index);
    HOSTDEVICE int getNeighborIndex(const int &nodeIndex, int &nIndex, const int &x, const int &y, const int &z);

    HOSTDEVICE void setInvalidNode(const int &index, bool &invalidNodeFound);
    HOSTDEVICE bool isNeighborInvalid(const int &index);

    HOST void removeInvalidNodes();

    HOSTDEVICE bool isStopper(int index) const;
private:
    HOSTDEVICE bool previousCellHasFluid(int index) const;

};

#endif
