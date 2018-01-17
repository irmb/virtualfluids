#ifndef GRID_H
#define GRID_H

#include "GridGenerator/global.h"


#include <GridGenerator/grid/distributions/Distribution.h>

struct Geometry;
struct Vertex;
struct Triangle;
class GridStrategy;

extern CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

struct VF_PUBLIC Grid : enableSharedFromThis<Grid>
{
    real sX = 0.0, sY = 0.0, sZ = 0.0;
    real eX, eY, eZ;
    real delta = 1.0;

	char *field;
	unsigned int nx, ny, nz;
	unsigned int size;
	Distribution d;

    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;
    unsigned int *matrixIndex;
    int reducedSize;

    HOST Grid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution &d);
    HOST Grid();
    HOST Grid(char *field, int startX, int startY, int startZ, int nx, int ny, int nz, Distribution &d);
    static HOST SPtr<Grid> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, std::shared_ptr<GridStrategy> gridStrategy, Distribution &d);

    HOST void mesh(Geometry &geometry);
    HOST void freeMemory();

	HOSTDEVICE bool isFluid(int index) const;
	HOSTDEVICE bool isSolid(int index) const;
	HOSTDEVICE bool isQ(int index) const;
	HOSTDEVICE bool isRb(int index) const;
	HOSTDEVICE void setFieldEntryToFluid(unsigned int index);
	HOSTDEVICE void setFieldEntryToSolid(unsigned int index);
	HOSTDEVICE void setFieldEntry(const Vertex &v, char val);
	HOSTDEVICE char getFieldEntry(const Vertex &v) const;
	HOSTDEVICE int transCoordToIndex(const real &x, const real &y, const real &z) const;
	HOSTDEVICE int transCoordToIndex(const Vertex &v) const;
	HOSTDEVICE void transIndexToCoords(const int index, real &x, real &y, real &z) const;
	HOSTDEVICE void print() const;
	HOSTDEVICE void setDebugPoint(const Vertex &actualPoint, const int pointValue);
	HOSTDEVICE bool isOutOfRange(const Vertex &actualPoint) const;

	/*---------------------------------------------------------------------------------*/
	HOSTDEVICE void meshTriangle(const Triangle &actualTriangle);
    HOSTDEVICE void meshTriangleExact(const Triangle &triangle);

	HOSTDEVICE void calculateQs(const Vertex &point, const Triangle &actualTriangle);

	char* toString(const char* name) const;

    /*---------------------------------------------------------------------------------*/
    HOSTDEVICE void setNeighborIndices(const int &index);
	HOSTDEVICE void getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, const real x, const real y, const real z) const;
    HOSTDEVICE void findNeighborIndex(int index);
    HOSTDEVICE int getNeighborIndex(const int &nodeIndex, int &neighborIndex, const real &expectedX, const real &expectedY, const real &expectedZ);

    HOSTDEVICE void setInvalidNode(const int &index, bool &invalidNodeFound);
    HOSTDEVICE bool isNeighborInvalid(const int &index);

    HOST void removeInvalidNodes();

    HOSTDEVICE bool isStopper(int index) const;
private:
    HOSTDEVICE bool previousCellHasFluid(int index) const;

    SPtr<GridStrategy> gridStrategy;


};

#endif
