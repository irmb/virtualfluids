#ifndef GRID_IMP_H
#define GRID_IMP_H

#include "GridGenerator/global.h"
#include "distributions/Distribution.h"

#include "Grid.h"
#include "Cell.h"
#include "Field.h"

#define DIR_END_MAX 27

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
    HOST GridImp(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d);
    HOST GridImp();
    HOST GridImp(Object* object, real delta, SPtr<GridStrategy> gridStrategy, Distribution d);

    HOST void initalNumberOfNodesAndSize();
    HOST void initalBoundingBoxStartValues();

public:
    virtual HOSTDEVICE ~GridImp();
    static HOST SPtr<GridImp> makeShared(Object* object, real delta, SPtr<GridStrategy> gridStrategy, Distribution d);

    static HOST SPtr<GridImp> makeShared(real startX, real startY, real startZ, real endX, real endY, real endZ, real delta, SPtr<GridStrategy> gridStrategy, Distribution d);

    real startX = 0.0, startY = 0.0, startZ = 0.0;
    real endX, endY, endZ;
    real delta = 1.0;

	uint nx, ny, nz;

private:
	uint size;
    uint sparseSize;
    bool periodicityX = true, periodicityY = true, periodicityZ = true;

    Field field;
    Object* object;

    HOSTDEVICE Cell getOddCellFromIndex(uint index) const;
public:
    int *neighborIndexX, *neighborIndexY, *neighborIndexZ;

    int *sparseIndices;
    Distribution distribution;
    SPtr<GridStrategy> gridStrategy;

    HOSTDEVICE real getDelta() const;
    HOSTDEVICE uint getSize() const;
    HOSTDEVICE Field getField() const;
    HOSTDEVICE uint getSparseSize() const;

    HOSTDEVICE void setFieldEntry(const Vertex &v, char val);
    HOSTDEVICE char getFieldEntry(const Vertex &v) const;

    HOSTDEVICE void findInnerNode(uint index);
    HOSTDEVICE void findStopperNode(uint index);
    HOST void mesh(Geometry &geometry);

    HOST void freeMemory();
    HOST void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ);

    HOST void findGridInterface(SPtr<Grid> grid);

    HOSTDEVICE void findGridInterfaceCF(uint index, GridImp& finerGrid);
    HOSTDEVICE void findGridInterfaceFC(uint index, GridImp& finerGrid);
    HOSTDEVICE void findOverlapStopper(uint index, GridImp& finerGrid);

	HOSTDEVICE int transCoordToIndex(const real &x, const real &y, const real &z) const;
	HOSTDEVICE int transCoordToIndex(const Vertex &v) const;
	HOSTDEVICE void transIndexToCoords(int index, real &x, real &y, real &z) const;
	HOSTDEVICE void print() const;
	HOSTDEVICE void setDebugPoint(const Vertex &actualPoint, const int pointValue);
	HOSTDEVICE bool isOutOfRange(const Vertex &actualPoint) const;

	/*---------------------------------------------------------------------------------*/
    HOSTDEVICE void meshTriangle(const Triangle &triangle);
	HOSTDEVICE void calculateQs(const Vertex &point, const Triangle &actualTriangle);
    /*---------------------------------------------------------------------------------*/
    HOSTDEVICE void setNeighborIndices(const int &index);
	HOSTDEVICE void getNeighborCoords(real &neighborX, real &neighborY, real &neighborZ, real x, real y, real z) const;
    HOSTDEVICE void findNeighborIndex(int index);
    HOSTDEVICE void findForGridInterfaceNewIndexCF(uint index);
    HOSTDEVICE void findForGridInterfaceNewIndexFC(uint index);


    HOSTDEVICE void setInvalidNode(const int &index, bool &invalidNodeFound);
    HOSTDEVICE bool isNeighborInvalid(const int &index) const;

    HOST void removeInvalidNodes();

    //HOSTDEVICE bool isStopper(int index) const;
    HOSTDEVICE bool isValidStartOfGridStopper(uint index) const;
    HOSTDEVICE bool isValidEndOfGridStopper(uint index) const;

    HOSTDEVICE int getSparseIndex(uint matrixIndex) const;

    HOST real* getDistribution() const;
    HOST int* getDirection() const;
    HOST int getStartDirection() const;
    HOST int getEndDirection() const;

    HOST void setNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *geo) const;

    HOSTDEVICE uint getNumberOfNodesCF() const;
    HOSTDEVICE uint getNumberOfNodesFC() const;
    HOST void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const;
    HOST uint* getCF_coarse() const;
    HOST uint* getCF_fine() const;
    HOST uint* getFC_coarse() const;
    HOST uint* getFC_fine() const;

    HOST std::string toString() const;
    HOSTDEVICE void setCellTo(uint index, char type);

private:
    HOSTDEVICE uint getXIndex(real x) const;
    HOSTDEVICE uint getYIndex(real y) const;
    HOSTDEVICE uint getZIndex(real z) const;


    static void setGridInterface(uint* gridInterfaceList, const uint* oldGridInterfaceList, uint size);

    HOSTDEVICE bool isOverlapStopper(uint index) const;
    HOSTDEVICE bool nodeInNextCellIs(int index, char type) const;
    HOSTDEVICE bool nodeInPreviousCellIs(int index, char type) const;
    HOSTDEVICE int getNeighborIndex(const real &expectedX, const real &expectedY, const real &expectedZ) const;
    HOSTDEVICE void setStopperNeighborCoords(int index);

    HOSTDEVICE real getNeighborCoord(bool periodicity, real endCoord, real coords[3], int direction) const;
    HOSTDEVICE real getFirstFluidNode(real coords[3], int direction, real startCoord) const;

public:
    HOSTDEVICE real getStartX() const override;
    HOSTDEVICE real getStartY() const override;
    HOSTDEVICE real getStartZ() const override;
    HOSTDEVICE real getEndX() const override;
    HOSTDEVICE real getEndY() const override;
    HOSTDEVICE real getEndZ() const override;
    HOSTDEVICE uint getNumberOfNodesX() const override;
    HOSTDEVICE uint getNumberOfNodesY() const override;
    HOSTDEVICE uint getNumberOfNodesZ() const override;
    SPtr<GridStrategy> getGridStrategy() const override;
    void setStartX(real startX) override;
    void setStartY(real startY) override;
    void setStartZ(real startZ) override;
    void setEndX(real endX) override;
    void setEndY(real endY) override;
    void setEndZ(real endZ) override;

    int* getNeighborsX() const override;
    int* getNeighborsY() const override;
    int* getNeighborsZ() const override;
    void inital() override;
    HOSTDEVICE char getFieldEntry(uint index) const override;
private:
    GridInterface* gridInterface;

    friend class GridGpuStrategy;
    friend class GridCpuStrategy;

};

#endif
