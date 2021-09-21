#ifndef GridBuilder_H
#define GridBuilder_H

#include <vector>
#include <string>
#include <memory>

#include "global.h"

#define GEOMQS 6
#define INLETQS 0
#define OUTLETQS 1
#define TOPQS 2
#define BOTTOMQS 3
#define FRONTQS 4
#define BACKQS 5

#define QFILES 7

struct Vertex;
class GridWrapper;
class Transformator;
class ArrowTransformator;
class PolyDataWriterWrapper;

class BoundingBox;
class Grid;

enum class SideType;

class BoundaryCondition;
class GeometryBoundaryCondition;

class GridBuilder
{
public:
    enum class GenerationDevice
    {
        CPU, GPU
    };

    virtual GRIDGENERATOR_EXPORT ~GridBuilder() {}
    virtual void getGridInformations(std::vector<int>& gridX, std::vector<int>& gridY, std::vector<int>& gridZ, std::vector<int>& distX, std::vector<int>& distY, std::vector<int>& distZ) = 0;
    virtual GRIDGENERATOR_EXPORT uint getNumberOfGridLevels() const = 0;


    virtual void writeArrows(std::string fileName) const = 0;

	virtual SPtr<Grid> getGrid(uint level) = 0;

    virtual unsigned int getNumberOfNodes(unsigned int level) const = 0;
    virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, 
                               uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, 
                               uint *geo, const int level) const = 0;
    virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const = 0;
    virtual uint getNumberOfNodesCF(int level) = 0;
    virtual uint getNumberOfNodesFC(int level) = 0;
    virtual void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const = 0;
    virtual void getGridInterfaceIndicesFCBorderBulk(uint *iCellFccBorder, uint *&iCellFccBulk, uint *iCellFcfBorder,
                                                      uint *&iCellFcfBulk, uint &intFCBorderKfc, uint &intFCBulkKfc,
                                                      int level) const           = 0;

    virtual void getOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) = 0;
    virtual void getOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) = 0;

    virtual uint getVelocitySize(int level) const = 0;
    virtual void getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const = 0;
    virtual uint getPressureSize(int level) const = 0;
    virtual void getPressureValues(real* rho, int* indices, int* neighborIndices, int level) const = 0;
    virtual void getVelocityQs(real* qs[27], int level) const = 0;
    virtual void getPressureQs(real* qs[27], int level) const = 0;

    virtual uint getGeometrySize(int level) const = 0;
    virtual void getGeometryIndices(int* indices, int level) const = 0;
    virtual void getGeometryQs(real* qs[27], int level) const = 0;
    virtual bool hasGeometryValues() const = 0;
    virtual void getGeometryValues(real* vx, real* vy, real* vz, int level) const = 0;

    virtual uint getCommunicationProcess(int direction) = 0;

    virtual SPtr<BoundaryCondition> getBoundaryCondition( SideType side, uint level ) const = 0;

    virtual SPtr<GeometryBoundaryCondition> getGeometryBoundaryCondition( uint level ) const = 0;

    virtual uint getNumberOfSendIndices( int direction, uint level ) = 0;
    virtual uint getNumberOfReceiveIndices( int direction, uint level ) = 0;
    virtual void getSendIndices( int* sendIndices, int direction, int level ) = 0;
    virtual void getReceiveIndices( int* sendIndices, int direction, int level ) = 0;

    virtual uint getNumberOfFluidNodes(unsigned int level) const = 0;
    virtual void getFluidNodeIndices(uint *fluidNodeIndices, const int level) const = 0;
    virtual uint getNumberOfFluidNodesBorder(unsigned int level) const = 0;
    virtual void getFluidNodeIndicesBorder(uint *fluidNodeIndices, const int level) const = 0;

};

#endif

