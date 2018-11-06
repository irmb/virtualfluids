#ifndef LEVEL_GRID_BUILDER_H
#define LEVEL_GRID_BUILDER_H

#include <vector>
#include <string>
#include <memory>
#include <array>

#include "global.h"

#include "grid/GridBuilder/GridBuilder.h"
#include "grid/Grid.h"
#include "grid/GridInterface.h"
#include "grid/NodeValues.h"

struct Vertex;
class  Grid;
class Transformator;
class ArrowTransformator;
class PolyDataWriterWrapper;
class BoundingBox;
enum class Device;

class Side;
class VelocityBoundaryCondition;
class PressureBoundaryCondition;
class GeometryBoundaryCondition;
enum class SideType;



class LevelGridBuilder : public GridBuilder
{
protected:
    VF_PUBLIC LevelGridBuilder(Device device, const std::string& d3qxx);

public:
    VF_PUBLIC static std::shared_ptr<LevelGridBuilder> makeShared(Device device, const std::string& d3qxx);

    VF_PUBLIC SPtr<Grid> getGrid(uint level) override;

    VF_PUBLIC void copyDataFromGpu();
    VF_PUBLIC virtual ~LevelGridBuilder();

    VF_PUBLIC void setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz);
    VF_PUBLIC void setPressureBoundaryCondition(SideType sideType, real rho);
    VF_PUBLIC void setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z);
    VF_PUBLIC void setNoSlipBoundaryCondition(SideType sideType);

    VF_PUBLIC void setCommunicationProcess(int direction, uint process);

    VF_PUBLIC uint getCommunicationProcess(int direction) override;

    VF_PUBLIC virtual std::shared_ptr<Grid> getGrid(int level, int box);


    VF_PUBLIC virtual unsigned int getNumberOfNodes(unsigned int level) const;


    VF_PUBLIC virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, 
                                         uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, 
                                         uint *geo, const int level) const override;
    VF_PUBLIC virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const;


    VF_PUBLIC uint getVelocitySize(int level) const;
    VF_PUBLIC virtual void getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const;
    VF_PUBLIC virtual void getVelocityQs(real* qs[27], int level) const;
    VF_PUBLIC uint getPressureSize(int level) const override;
    VF_PUBLIC void getPressureValues(real* rho, int* indices, int* neighborIndices, int level) const override;
    VF_PUBLIC virtual void getPressureQs(real* qs[27], int level) const;

    VF_PUBLIC virtual void getGeometryQs(real* qs[27], int level) const;
    VF_PUBLIC virtual uint getGeometrySize(int level) const;
    VF_PUBLIC virtual void getGeometryIndices(int* indices, int level) const;
    VF_PUBLIC virtual bool hasGeometryValues() const;
    VF_PUBLIC virtual void getGeometryValues(real* vx, real* vy, real* vz, int level) const;


    VF_PUBLIC void writeArrows(std::string fileName) const;

    VF_PUBLIC SPtr<BoundaryCondition> getBoundaryCondition( SideType side, uint level ) const override;
    VF_PUBLIC SPtr<GeometryBoundaryCondition> getGeometryBoundaryCondition(uint level) const override;

protected:
    

    struct BoundaryConditions
    {
		BoundaryConditions() : geometryBoundaryCondition(nullptr) {}

        std::vector<SPtr<VelocityBoundaryCondition> > velocityBoundaryConditions;
        std::vector<SPtr<PressureBoundaryCondition> > pressureBoundaryConditions;

        std::vector<SPtr<VelocityBoundaryCondition> > noSlipBoundaryConditions;

        SPtr<GeometryBoundaryCondition> geometryBoundaryCondition;
    };
    bool geometryHasValues = false;

    std::vector<std::shared_ptr<Grid> > grids;
    std::vector<SPtr<BoundaryConditions> > boundaryConditions;

    std::array<uint, 6> communicationProcesses;

    void checkLevel(int level);

protected:
    void setVelocityGeometryBoundaryCondition(real vx, real vy, real vz);

    void createBCVectors();
    void addShortQsToVector(int index);
    void addQsToVector(int index);
    void fillRBForNode(int index, int direction, int directionSign, int rb);

    Vertex getVertex(const int matrixIndex) const;

private:
    Device device;
    std::string d3qxx;

public:
    VF_PUBLIC void getGridInformations(std::vector<int>& gridX, std::vector<int>& gridY,
                                       std::vector<int>& gridZ, std::vector<int>& distX, std::vector<int>& distY,
                                       std::vector<int>& distZ) override;
    VF_PUBLIC uint getNumberOfGridLevels() const override;

    VF_PUBLIC uint getNumberOfNodesCF(int level) override;
    VF_PUBLIC uint getNumberOfNodesFC(int level) override;

    VF_PUBLIC void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const override;

    VF_PUBLIC void getOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) override;
    VF_PUBLIC void getOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) override;

    VF_PUBLIC uint getNumberOfSendIndices( int direction, uint level ) override;
    VF_PUBLIC uint getNumberOfReceiveIndices( int direction, uint level ) override;
    VF_PUBLIC void getSendIndices( int* sendIndices, int direction, int level ) override;
    VF_PUBLIC void getReceiveIndices( int* sendIndices, int direction, int level ) override;

};

#endif

