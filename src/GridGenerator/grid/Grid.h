#ifndef GRID_H
#define GRID_H

#include "GridGenerator/global.h"

struct Geometry;
struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;


class VF_PUBLIC Grid
{
public:
    HOSTDEVICE virtual ~Grid() {}

    HOSTDEVICE virtual real getDelta() const = 0;
    HOSTDEVICE virtual uint getReducedSize() const = 0;
    HOSTDEVICE virtual uint getSize() const = 0;

    HOSTDEVICE virtual real getStartX() const = 0;
    HOSTDEVICE virtual real getStartY() const = 0;
    HOSTDEVICE virtual real getStartZ() const = 0;

    HOSTDEVICE virtual real getEndX() const = 0;
    HOSTDEVICE virtual real getEndY() const = 0;
    HOSTDEVICE virtual real getEndZ() const = 0;

    virtual void setStartX(real startX) = 0;
    virtual void setStartY(real startY) = 0;
    virtual void setStartZ(real startZ) = 0;

    virtual void setEndX(real endX) = 0;
    virtual void setEndY(real endY) = 0;
    virtual void setEndZ(real endZ) = 0;

    HOSTDEVICE virtual uint getNumberOfNodesX() const = 0;
    HOSTDEVICE virtual uint getNumberOfNodesY() const = 0;
    HOSTDEVICE virtual uint getNumberOfNodesZ() const = 0;

    HOSTDEVICE virtual uint getNumberOfNodesCF() const = 0;
    HOSTDEVICE virtual uint getNumberOfNodesFC() const = 0;

    HOSTDEVICE virtual int getIndex(uint matrixIndex) const = 0;
    HOSTDEVICE virtual char getFieldEntry(uint matrixIndex) const = 0;
    virtual void setFieldEntry(uint index, char entry) = 0;

    HOST virtual void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const = 0;

    HOST virtual int *getNeighborsX() const = 0;
    HOST virtual int *getNeighborsY() const = 0;
    HOST virtual int *getNeighborsZ() const = 0;

    HOST virtual uint* getCF_coarse() const = 0;
    HOST virtual uint* getCF_fine() const = 0;
    HOST virtual uint* getFC_coarse() const = 0;
    HOST virtual uint* getFC_fine() const = 0;

    HOST virtual real* getDistribution() const = 0;
    HOST virtual int* getDirection() const = 0;
    HOST virtual int getStartDirection() const = 0;
    HOST virtual int getEndDirection() const = 0;

    HOST virtual void setNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *neighborX, unsigned int *neighborY, unsigned int *neighborZ, unsigned int *geo) const = 0;

    HOST virtual SPtr<GridStrategy> getGridStrategy() const = 0;
    HOSTDEVICE virtual void transIndexToCoords(int index, real &x, real &y, real &z) const = 0;
    HOSTDEVICE virtual int transCoordToIndex(const Vertex &v) const = 0;
    HOSTDEVICE virtual int transCoordToIndex(const real &x, const real &y, const real &z) const = 0;

    HOST virtual void allocateGridMemory() = 0;

    HOST virtual void findGridInterface(SPtr<Grid> grid) = 0;
    HOST virtual void mesh(Geometry& geometry) = 0;


    HOST virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) = 0;
    HOST virtual void freeMemory() = 0;

};

#endif
