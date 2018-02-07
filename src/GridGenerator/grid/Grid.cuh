#ifndef GRID_H
#define GRID_H

#include "GridGenerator/global.h"
#include "distributions/Distribution.h"

#define DIR_END_MAX 27

struct Geometry;
struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;

extern CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

class VF_PUBLIC Grid : public enableSharedFromThis<Grid>
{
private:
    HOST Grid(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d);
    HOST Grid();

public:
    HOST ~Grid();
    static HOST SPtr<Grid> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, std::shared_ptr<GridStrategy> gridStrategy, Distribution d);

    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX, endY, endZ;
    real delta = 1.0;

	char *field;
	uint nx, ny, nz;

private:
	uint size;
    uint reducedSize;
    bool periodicityX = true, periodicityY = true, periodicityZ = true;


public:
    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;

    int *matrixIndex;
    Distribution distribution;
    SPtr<GridStrategy> gridStrategy;


    HOSTDEVICE uint getSize() const;
    HOSTDEVICE uint getReducedSize() const;

    HOSTDEVICE void initalNodes(uint index);
    HOST void mesh(Geometry &geometry);
    HOST void freeMemory();
    HOST void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ);

    HOST void removeOverlapNodes(SPtr<Grid> grid);
    HOSTDEVICE void createGridInterface(uint index, const Grid& finerGrid);


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
    HOSTDEVICE void meshTriangle(const Triangle &triangle);
	HOSTDEVICE void calculateQs(const Vertex &point, const Triangle &actualTriangle);
    /*---------------------------------------------------------------------------------*/
    HOSTDEVICE void setNeighborIndices(const int &index);
	HOSTDEVICE void getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    HOSTDEVICE real getNeighhborCoord(bool periodicity, real actualCoord, real startCoord, real  endCoord) const;
    HOSTDEVICE void findNeighborIndex(int index);
    HOSTDEVICE void findForGridInterfaceNewIndexCF(uint index);
    HOSTDEVICE void findForGridInterfaceNewIndexFC(uint index);


    HOSTDEVICE void setInvalidNode(const int &index, bool &invalidNodeFound);
    HOSTDEVICE bool isNeighborInvalid(const int &index);

    HOST void removeInvalidNodes();

    //HOSTDEVICE bool isStopper(int index) const;
    HOSTDEVICE bool isEndOfGridStopper(uint index) const;



    HOSTDEVICE uint getNumberOfNodesCF() const;
    HOSTDEVICE uint getNumberOfNodesFC() const;
    HOST void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const;
    HOST uint* getCF_coarse() const;
    HOST uint* getCF_fine() const;
    HOST uint* getFC_coarse() const;
    HOST uint* getFC_fine() const;

    HOST std::string toString() const;

private:
    static void setGridInterface(uint* gridInterfaceList, const uint* oldGridInterfaceList, uint size);

    HOSTDEVICE void findGridInterface(uint index, const Grid& finerGrid);
    HOSTDEVICE void setOverlapNodeToInvalid(uint index, const Grid& finerGrid);
    HOSTDEVICE bool isInside(uint index, const Grid& grid) const;
    HOSTDEVICE bool isOverlapStopper(uint index) const;
    HOSTDEVICE bool nodeInNextCellIsInvalid(int index) const;
    HOSTDEVICE int getNeighborIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const;
    HOSTDEVICE void setStopperNeighborCoords(int index);


    GridInterface* gridInterface;

    friend class GridGpuStrategy;
    friend class GridCpuStrategy;

};

#endif
