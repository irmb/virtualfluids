#ifndef GRID_IMP_H
#define GRID_IMP_H

#include "GridGenerator/global.h"
#include "distributions/Distribution.h"

#include "core/LbmOrGks.h"

#include "Grid.h"
#include "Cell.h"
#include "Field.h"


class TriangularMesh;
struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;
class Object;
class BoundingBox;
class TriangularMeshDiscretizationStrategy;

extern CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

class VF_PUBLIC GridImp : public enableSharedFromThis<GridImp>, public Grid
{
private:
    HOST GridImp();
    HOST GridImp(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level);

public:
    virtual HOSTDEVICE ~GridImp();
    static HOST SPtr<GridImp> makeShared(Object* object, real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d, uint level);

private:
    HOST void initalNumberOfNodesAndSize();
    HOSTDEVICE Cell getOddCellFromIndex(uint index) const;
    HOSTDEVICE bool isValidInnerStopper(uint index) const;
    HOSTDEVICE bool isValidEndOfGridStopper(uint index) const;
    HOSTDEVICE bool isValidEndOfGridBoundaryStopper(uint index) const;
    HOSTDEVICE bool isOutSideOfGrid(Cell &cell) const;
    HOSTDEVICE bool contains(Cell &cell, char type) const;
    HOSTDEVICE void setNodeTo(Cell &cell, char type);

    HOSTDEVICE bool nodeInPreviousCellIs(int index, char type) const;
    HOSTDEVICE bool nodeInCellIs(Cell& cell, char type) const override;

    HOSTDEVICE uint getXIndex(real x) const;
    HOSTDEVICE uint getYIndex(real y) const;
    HOSTDEVICE uint getZIndex(real z) const;

    uint level;

    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX, endY, endZ;
    real delta = 1.0;

	uint nx, ny, nz;

	uint size;
    uint sparseSize;
    bool periodicityX = false, periodicityY = false, periodicityZ = false;

    Field field;
    Object* object;
    GridInterface* gridInterface;

    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;
    int *sparseIndices;
    SPtr<GridStrategy> gridStrategy;
    TriangularMeshDiscretizationStrategy* triangularMeshDiscretizationStrategy;


public:
    HOST void inital() override;
    HOSTDEVICE void removeOddBoundaryCellNode(uint index);

    HOST void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) override;
    void setPeriodicityX(bool periodicity) override;
    void setPeriodicityY(bool periodicity) override;
    void setPeriodicityZ(bool periodicity) override;

    HOSTDEVICE void setCellTo(uint index, char type);
    HOSTDEVICE void setNonStopperOutOfGridCellTo(uint index, char type);

    HOSTDEVICE uint transCoordToIndex(const real &x, const real &y, const real &z) const override;
    HOSTDEVICE void transIndexToCoords(uint index, real &x, real &y, real &z) const override;

    HOST virtual void findGridInterface(SPtr<Grid> grid, LbmOrGks lbmOrGks) override;
    HOST void freeMemory() override;

    HOST uint getLevel(real levelNull) const;
    HOST uint getLevel() const;
    HOST void setTriangularMeshDiscretizationStrategy(TriangularMeshDiscretizationStrategy* triangularMeshDiscretizationStrategy);
    
public:
    Distribution distribution;

    HOSTDEVICE void initalNodeToOutOfGrid(uint index);

    HOSTDEVICE void findInnerNode(uint index);

    bool isInside(const Cell& cell) const;

    HOSTDEVICE void findStopperNode(uint index);
	HOSTDEVICE void findEndOfGridStopperNode(uint index);
	HOSTDEVICE void findSolidStopperNode(uint index);

    HOSTDEVICE void findGridInterfaceCF(uint index, GridImp& finerGrid, LbmOrGks lbmOrGks);
    HOSTDEVICE void findGridInterfaceFC(uint index, GridImp& finerGrid);
    HOSTDEVICE void findOverlapStopper(uint index, GridImp& finerGrid);

    HOSTDEVICE void setNodeTo(uint index, char type);
    HOSTDEVICE bool isNode(uint index, char type) const;
    HOSTDEVICE bool nodeInNextCellIs(int index, char type) const;
    HOSTDEVICE bool hasAllNeighbors(uint index) const;
    HOSTDEVICE bool hasNeighborOfType(uint index, char type)const;

    HOSTDEVICE Field getField() const;
    HOSTDEVICE char getFieldEntry(uint index) const override;
    HOSTDEVICE void setFieldEntry(uint matrixIndex, char type) override;


    HOSTDEVICE real getDelta() const override;
    HOSTDEVICE uint getSize() const override;
    HOSTDEVICE uint getSparseSize() const override;
    HOSTDEVICE int getSparseIndex(uint matrixIndex) const override;
    HOST real* getDistribution() const override;
    HOST int* getDirection() const override;
    HOST int getStartDirection() const override;
    HOST int getEndDirection() const override;

    HOSTDEVICE Vertex getMinimumOnNode(Vertex exact) const;
    HOSTDEVICE Vertex getMaximumOnNode(Vertex exact) const;

    HOSTDEVICE real getStartX() const override;
    HOSTDEVICE real getStartY() const override;
    HOSTDEVICE real getStartZ() const override;
    HOSTDEVICE real getEndX() const override;
    HOSTDEVICE real getEndY() const override;
    HOSTDEVICE real getEndZ() const override;
    HOSTDEVICE uint getNumberOfNodesX() const override;
    HOSTDEVICE uint getNumberOfNodesY() const override;
    HOSTDEVICE uint getNumberOfNodesZ() const override;
    HOST void getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *geo) const override;

    HOSTDEVICE uint getNumberOfNodesCF() const override;
    HOSTDEVICE uint getNumberOfNodesFC() const override;
    HOST void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const override;
    HOST static void getGridInterface(uint* gridInterfaceList, const uint* oldGridInterfaceList, uint size);

    int* getNeighborsX() const override;
    int* getNeighborsY() const override;
    int* getNeighborsZ() const override;

    HOST uint* getCF_coarse() const override;
    HOST uint* getCF_fine() const override;
    HOST uint* getFC_coarse() const override;
    HOST uint* getFC_fine() const override;

    SPtr<GridStrategy> getGridStrategy() const override;


    HOSTDEVICE void print() const;


public:
    HOST virtual void findSparseIndices(SPtr<Grid> fineGrid);

    HOST void updateSparseIndices();
    HOSTDEVICE void setNeighborIndices(uint index);
    HOSTDEVICE real getFirstFluidNode(real coords[3], int direction, real startCoord) const;
    HOSTDEVICE real getLastFluidNode(real coords[3], int direction, real startCoord) const;
private:
    HOSTDEVICE void setStopperNeighborCoords(uint index);
    HOSTDEVICE void getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    HOSTDEVICE real getNeighborCoord(bool periodicity, real endCoord, real coords[3], int direction) const;
    

    HOSTDEVICE int getSparseIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const;

    HOSTDEVICE static real getMinimumOnNodes(const real& minExact, const real& decimalStart, const real& delta);
    HOSTDEVICE static real getMaximumOnNodes(const real& maxExact, const real& decimalStart, const real& delta);

public:
    HOSTDEVICE BoundingBox getBoundingBoxOnNodes(Triangle &triangle) const;

    HOST void mesh(Object* object) override;

    HOST void mesh(TriangularMesh &geometry) override;
    HOSTDEVICE void mesh(Triangle &triangle);


    HOST void findQs(Object* object) override;
    HOST void findQs(TriangularMesh &triangularMesh);
    HOSTDEVICE void findQs(Triangle &triangle);
private:
    HOSTDEVICE void setDebugPoint(uint index, int pointValue);
	HOSTDEVICE void calculateQs(const Vertex &point, const Triangle &actualTriangle) const;



private:
    //HOSTDEVICE bool isNeighborInside(const int &index) const;


private:
    friend class GridGpuStrategy;
    friend class GridCpuStrategy;
};

#endif
