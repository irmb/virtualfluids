#ifndef GRID_H
#define GRID_H

#include "GridGenerator/global.h"

#include "core/LbmOrGks.h"

#include "Cell.h"
#include "geometries/Vertex/Vertex.h"

class TriangularMesh;
struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;
class Object;

class VF_PUBLIC Grid
{
public:
    HOSTDEVICE virtual ~Grid() {}

    HOSTDEVICE virtual const Object* getObject() const = 0;

    HOSTDEVICE virtual real getDelta() const = 0;
    HOSTDEVICE virtual uint getSparseSize() const = 0;
    HOSTDEVICE virtual uint getSize() const = 0;

    HOSTDEVICE virtual real getStartX() const = 0;
    HOSTDEVICE virtual real getStartY() const = 0;
    HOSTDEVICE virtual real getStartZ() const = 0;

    HOSTDEVICE virtual real getEndX() const = 0;
    HOSTDEVICE virtual real getEndY() const = 0;
    HOSTDEVICE virtual real getEndZ() const = 0;

    HOSTDEVICE virtual Vertex getMinimumOnNode(Vertex exact) const = 0;
    HOSTDEVICE virtual Vertex getMaximumOnNode(Vertex exact) const = 0;

    HOSTDEVICE virtual uint getNumberOfNodesX() const = 0;
    HOSTDEVICE virtual uint getNumberOfNodesY() const = 0;
    HOSTDEVICE virtual uint getNumberOfNodesZ() const = 0;

    HOSTDEVICE virtual uint getNumberOfNodesCF() const = 0;
    HOSTDEVICE virtual uint getNumberOfNodesFC() const = 0;

    HOSTDEVICE virtual int getSparseIndex(uint matrixIndex) const = 0;
    HOSTDEVICE virtual char getFieldEntry(uint matrixIndex) const = 0;
    HOSTDEVICE virtual void setFieldEntry(uint matrixIndex, char type) = 0;

    HOST virtual void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const = 0;

    HOST virtual int *getNeighborsX() const = 0;
    HOST virtual int *getNeighborsY() const = 0;
    HOST virtual int *getNeighborsZ() const = 0;

    HOST virtual uint* getCF_coarse() const = 0;
    HOST virtual uint* getCF_fine()   const = 0;
    HOST virtual uint* getCF_offset() const = 0;

    HOST virtual uint* getFC_coarse() const = 0;
    HOST virtual uint* getFC_fine()   const = 0;
    HOST virtual uint* getFC_offset() const = 0;

    HOST virtual real* getDistribution() const = 0;
    HOST virtual int* getDirection() const = 0;
    HOST virtual int getStartDirection() const = 0;
    HOST virtual int getEndDirection() const = 0;

    HOST virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *neighborX, unsigned int *neighborY, unsigned int *neighborZ, unsigned int *geo) const = 0;

    HOST virtual SPtr<GridStrategy> getGridStrategy() const = 0;
    HOSTDEVICE virtual void transIndexToCoords(uint index, real &x, real &y, real &z) const = 0;
    HOSTDEVICE virtual uint transCoordToIndex(const real &x, const real &y, const real &z) const = 0;

    HOST virtual void inital() = 0;
    HOST virtual void inital(const SPtr<Grid> fineGrid) = 0;
    
    HOST virtual void setOddStart( bool xOddStart, bool yOddStart, bool zOddStart ) = 0;

    HOST virtual void findGridInterface(SPtr<Grid> grid, LbmOrGks lbmOrGks) = 0;

    HOST virtual void mesh(TriangularMesh& geometry) = 0;
    HOST virtual void mesh(Object* object) = 0;

    HOST virtual void findQs(Object* object) = 0;

    HOST virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) = 0;
    HOST virtual void setPeriodicityX(bool periodicity) = 0;
    HOST virtual void setPeriodicityY(bool periodicity) = 0;
    HOST virtual void setPeriodicityZ(bool periodicity) = 0;
    HOST virtual void freeMemory() = 0;


    HOSTDEVICE virtual bool nodeInCellIs(Cell& cell, char type) const = 0;

    HOST virtual void findSparseIndices(SPtr<Grid> fineGrid) = 0;

    HOSTDEVICE virtual real getFirstFluidNode(real coords[3], int direction, real startCoord) const = 0;
    HOSTDEVICE virtual real getLastFluidNode(real coords[3], int direction, real startCoord) const = 0;

	HOST virtual uint getNumberOfSolidBoundaryNodes() const = 0;
	HOST virtual void setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes) = 0;

	HOST virtual real getQValue(const uint index, const uint dir) const = 0;
};

#endif
