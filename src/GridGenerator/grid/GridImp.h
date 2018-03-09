#ifndef GRID_IMP_H
#define GRID_IMP_H

#include "GridGenerator/global.h"
#include "distributions/Distribution.h"

#include "Grid.h"
#include "Cell.h"
#include "Field.h"


struct Geometry;
struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;
class Object;

extern CONSTANT int DIRECTIONS[DIR_END_MAX][DIMENSION];

class VF_PUBLIC GridImp : public enableSharedFromThis<GridImp>, public Grid
{
private:
    HOST GridImp();
    HOST GridImp(Object* object, real delta, SPtr<GridStrategy> gridStrategy, Distribution d);

public:
    virtual HOSTDEVICE ~GridImp();
    static HOST SPtr<GridImp> makeShared(Object* object, real delta, SPtr<GridStrategy> gridStrategy, Distribution d);

private:
    HOST void initalNumberOfNodesAndSize();
    HOST void initalBoundingBoxStartValues();
    HOSTDEVICE Cell getOddCellFromIndex(uint index) const;
    HOSTDEVICE bool isValidStartOfGridStopper(uint index) const;
    HOSTDEVICE bool isValidEndOfGridStopper(uint index) const;

    HOSTDEVICE bool nodeInNextCellIs(int index, char type) const;
    HOSTDEVICE bool nodeInPreviousCellIs(int index, char type) const;
    HOSTDEVICE bool nodeInCellIs(Cell& cell, char type) const override;

    HOSTDEVICE int getXIndex(real x) const;
    HOSTDEVICE int getYIndex(real y) const;
    HOSTDEVICE int getZIndex(real z) const;

    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX, endY, endZ;
    real delta = 1.0;

	uint nx, ny, nz;

	uint size;
    uint sparseSize;
    bool periodicityX = true, periodicityY = true, periodicityZ = true;

    Field field;
    Object* object;
    GridInterface* gridInterface;

    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;
    int *sparseIndices;
    SPtr<GridStrategy> gridStrategy;


public:
    HOST void inital() override;

    HOST void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) override;
    HOSTDEVICE void setCellTo(uint index, char type);

    HOSTDEVICE int transCoordToIndex(const real &x, const real &y, const real &z) const override;
    HOSTDEVICE void transIndexToCoords(int index, real &x, real &y, real &z) const override;

    HOST void findGridInterface(SPtr<Grid> grid) override;
    HOST void freeMemory() override;


public:
    Distribution distribution;

    HOSTDEVICE void findInnerNode(uint index);
    HOSTDEVICE void findStopperNode(uint index);

    HOSTDEVICE void findGridInterfaceCF(uint index, GridImp& finerGrid);
    HOSTDEVICE void findGridInterfaceFC(uint index, GridImp& finerGrid);
    HOSTDEVICE void findOverlapStopper(uint index, GridImp& finerGrid);


    HOSTDEVICE Field getField() const;
    HOSTDEVICE char getFieldEntry(uint index) const override;

    HOSTDEVICE real getDelta() const override;
    HOSTDEVICE uint getSize() const override;
    HOSTDEVICE uint getSparseSize() const override;
    HOSTDEVICE int getSparseIndex(uint matrixIndex) const override;
    HOST real* getDistribution() const override;
    HOST int* getDirection() const override;
    HOST int getStartDirection() const override;
    HOST int getEndDirection() const override;

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
    HOST void updateSparseIndices();
    HOSTDEVICE void setNeighborIndices(uint index);

    HOSTDEVICE void findForGridInterfaceSparseIndexCF(uint index);
    HOSTDEVICE void findForGridInterfaceSparseIndexFC(uint index);
private:
    HOSTDEVICE void setStopperNeighborCoords(uint index);
    HOSTDEVICE void getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    HOSTDEVICE real getNeighborCoord(bool periodicity, real endCoord, real coords[3], int direction) const;
    HOSTDEVICE real getFirstFluidNode(real coords[3], int direction, real startCoord) const;
    HOSTDEVICE int getSparseIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const;


public:
    HOST void mesh(Geometry &geometry) override;
    HOSTDEVICE void mesh(const Triangle &triangle);

private:
    HOSTDEVICE void setDebugPoint(uint index, int pointValue);
	HOSTDEVICE void calculateQs(const Vertex &point, const Triangle &actualTriangle) const;


public:
    HOSTDEVICE void setInvalidNode(const int &index, bool &invalidNodeFound);
private:
    HOSTDEVICE bool isNeighborInvalid(const int &index) const;


private:
    friend class GridGpuStrategy;
    friend class GridCpuStrategy;
};

#endif
