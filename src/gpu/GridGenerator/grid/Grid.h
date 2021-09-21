#ifndef GRID_H
#define GRID_H

#include "Core/LbmOrGks.h"

#include "global.h"

#include "geometries/Vertex/Vertex.h"

#include "grid/Cell.h"

class TriangularMesh;
struct Vertex;
struct Triangle;
class GridStrategy;
class GridInterface;
class Object;
class BoundingBox;

class GRIDGENERATOR_EXPORT Grid
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

    CUDA_HOST virtual void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf) const = 0;
    CUDA_HOST virtual void getGridInterfaceIndicesFCBorderBulk(uint *iCellFccBorder, uint *&iCellFccBulk,
                                                                uint *iCellFcfBorder, uint *&iCellFcfBulk,
                                                                uint &intFCBorderKfc, uint &intFCBulkKfc,
                                                                int level) const = 0;

    CUDA_HOST virtual int *getNeighborsX() const = 0;
    CUDA_HOST virtual int *getNeighborsY() const = 0;
    CUDA_HOST virtual int *getNeighborsZ() const = 0;
    CUDA_HOST virtual int *getNeighborsNegative() const = 0;

    CUDA_HOST virtual uint* getCF_coarse() const = 0;
    CUDA_HOST virtual uint* getCF_fine()   const = 0;
    CUDA_HOST virtual uint* getCF_offset() const = 0;

    CUDA_HOST virtual uint* getFC_coarse() const = 0;
    CUDA_HOST virtual uint* getFC_fine()   const = 0;
    CUDA_HOST virtual uint* getFC_offset() const = 0;

    CUDA_HOST virtual real* getDistribution() const = 0;
    CUDA_HOST virtual int* getDirection() const = 0;
    CUDA_HOST virtual int getStartDirection() const = 0;
    CUDA_HOST virtual int getEndDirection() const = 0;

    CUDA_HOST virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, uint *geo) const = 0;

    CUDA_HOST virtual SPtr<GridStrategy> getGridStrategy() const = 0;
    HOSTDEVICE virtual void transIndexToCoords(uint index, real &x, real &y, real &z) const = 0;
    HOSTDEVICE virtual uint transCoordToIndex(const real &x, const real &y, const real &z) const = 0;

    CUDA_HOST virtual void inital(const SPtr<Grid> fineGrid, uint numberOfLayers) = 0;
    
    CUDA_HOST virtual void setOddStart( bool xOddStart, bool yOddStart, bool zOddStart ) = 0;

    CUDA_HOST virtual void findGridInterface(SPtr<Grid> grid, LbmOrGks lbmOrGks) = 0;

    HOSTDEVICE virtual void repairGridInterfaceOnMultiGPU(SPtr<Grid> fineGrid) = 0;

    CUDA_HOST virtual void limitToSubDomain(SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks) = 0;
    
    CUDA_HOST virtual void enableFindSolidBoundaryNodes() = 0;
    CUDA_HOST virtual void enableComputeQs() = 0;

    CUDA_HOST virtual void mesh(TriangularMesh& geometry) = 0;
    CUDA_HOST virtual void mesh(Object* object) = 0;

    CUDA_HOST virtual void closeNeedleCells() = 0;
    CUDA_HOST virtual void closeNeedleCellsThinWall() = 0;

    CUDA_HOST virtual void findQs(Object* object) = 0;

    CUDA_HOST virtual void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) = 0;
    CUDA_HOST virtual void setPeriodicityX(bool periodicity) = 0;
    CUDA_HOST virtual void setPeriodicityY(bool periodicity) = 0;
    CUDA_HOST virtual void setPeriodicityZ(bool periodicity) = 0;

    CUDA_HOST virtual bool getPeriodicityX() = 0;
    CUDA_HOST virtual bool getPeriodicityY() = 0;
    CUDA_HOST virtual bool getPeriodicityZ() = 0;

    CUDA_HOST virtual void setEnableFixRefinementIntoTheWall( bool enableFixRefinementIntoTheWall ) = 0;

    CUDA_HOST virtual void freeMemory() = 0;


    HOSTDEVICE virtual bool nodeInCellIs(Cell& cell, char type) const = 0;

    CUDA_HOST virtual void findSparseIndices(SPtr<Grid> fineGrid) = 0;

    HOSTDEVICE virtual real getFirstFluidNode(real coords[3], int direction, real startCoord) const = 0;
    HOSTDEVICE virtual real getLastFluidNode(real coords[3], int direction, real startCoord) const = 0;

	CUDA_HOST virtual uint getNumberOfSolidBoundaryNodes() const = 0;
	CUDA_HOST virtual void setNumberOfSolidBoundaryNodes(uint numberOfSolidBoundaryNodes) = 0;

	CUDA_HOST virtual real getQValue(const uint index, const uint dir) const = 0;
	CUDA_HOST virtual uint getQPatch(const uint index) const = 0;

    CUDA_HOST virtual void setInnerRegionFromFinerGrid( bool innerRegionFromFinerGrid ) = 0;

    CUDA_HOST virtual void setNumberOfLayers( uint numberOfLayers ) = 0;

    virtual void findCommunicationIndices(int direction, SPtr<BoundingBox> subDomainBox, LbmOrGks lbmOrGks) = 0;

    virtual uint getNumberOfSendNodes(int direction) = 0;
    virtual uint getNumberOfReceiveNodes(int direction)  = 0;

    virtual bool isSendNode(int index) const = 0;
    virtual bool isReceiveNode(int index) const = 0;

    virtual uint getSendIndex(int direction, uint index)  = 0;
    virtual uint getReceiveIndex(int direction, uint index)  = 0;

    virtual void repairCommunicationIndices(int direction) = 0;

    // needed for CUDA Streams 
    virtual void findFluidNodeIndices(bool onlyBulk) = 0;
    virtual uint getNumberOfFluidNodes() const = 0;
    virtual void getFluidNodeIndices(uint *fluidNodeIndices) const = 0;

    virtual void findFluidNodeIndicesBorder() = 0;

    virtual uint getNumberOfFluidNodesBorder() const = 0;
    virtual void getFluidNodeIndicesBorder(uint *fluidNodeIndicesBorder) const = 0;
};

#endif
