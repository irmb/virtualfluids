#ifndef GRID_H
#define GRID_H

#include "GridGenerator/global.h"


#include <GridGenerator/grid/distributions/Distribution.h>
#include "GridInterface.cuh"

struct Geometry;
struct Vertex;
struct Triangle;
class GridStrategy;


extern CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

struct VF_PUBLIC Grid : enableSharedFromThis<Grid>
{
    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX, endY, endZ;
    real delta = 1.0;

	char *field;
	uint nx, ny, nz;
	uint size;
	Distribution d;

    GridInterface gridInterface;

    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;

    int *matrixIndex;
    uint reducedSize;

    HOST Grid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution &d);
    HOST Grid();
    HOST Grid(char *field, int startX, int startY, int startZ, int nx, int ny, int nz, Distribution &d);
    
    static HOST SPtr<Grid> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, std::shared_ptr<GridStrategy> gridStrategy, Distribution &d);

    HOST void mesh(Geometry &geometry);
    HOST void freeMemory();

    HOST void removeOverlapNodes(SPtr<Grid> grid);
    HOSTDEVICE void setOverlapNodeToInvalid(uint index, const Grid& finerGrid);
    HOSTDEVICE bool isOverlapStopper(uint index) const;
    HOSTDEVICE bool isInside(uint index, const Grid& grid);

	HOSTDEVICE bool isFluid(uint index) const;
	HOSTDEVICE bool isSolid(uint index) const;
	HOSTDEVICE bool isQ(uint index) const;
    HOSTDEVICE bool isRb(uint index) const;
    HOSTDEVICE bool isInvalid(uint index) const;
	HOSTDEVICE void setFieldEntryToFluid(uint index);
	HOSTDEVICE void setFieldEntryToSolid(uint index);
    HOSTDEVICE void setFieldEntryToInvalid(uint index);
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
    HOSTDEVICE int getNeighborIndex(/*const int &nodeIndex, const int &neighborIndex, */const real &expectedX, const real &expectedY, const real &expectedZ);

    HOSTDEVICE void setInvalidNode(const int &index, bool &invalidNodeFound);
    HOSTDEVICE bool isNeighborInvalid(const int &index);

    HOST void removeInvalidNodes();

    HOSTDEVICE bool isStopper(int index) const;
private:
    HOSTDEVICE bool previousCellHasFluid(int index) const;

    HOSTDEVICE bool nodeInNextCellIsInvalid(int index) const;

    SPtr<GridStrategy> gridStrategy;


};

#endif
