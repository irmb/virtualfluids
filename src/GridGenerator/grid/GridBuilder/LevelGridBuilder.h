#ifndef LEVEL_GRID_BUILDER_H
#define LEVEL_GRID_BUILDER_H

#include "GridGenerator/global.h"


#include <vector>
#include <string>
#include <memory>

#include "GridBuilder.h"
#include "grid/GridInterface.cuh"

struct Vertex;
class  Grid;
class Transformator;
class ArrowTransformator;
class PolyDataWriterWrapper;

template <typename T>
class BoundingBox;
//class GridStub;
enum class Device;

class LevelGridBuilder : public GridBuilder
{
protected:
    VF_PUBLIC LevelGridBuilder(Device device, const std::string& d3qxx);

public:
    //VF_PUBLIC void setGrids(std::vector<SPtr<GridStub> > grids);

    VF_PUBLIC static std::shared_ptr<LevelGridBuilder> makeShared(Device device, const std::string& d3qxx);

    VF_PUBLIC void addGrid(real minX, real minY, real minZ, real maxX, real maxY, real maxZ, bool periodictyX, bool periodictyY, bool periodictyZ);
    VF_PUBLIC SPtr<Grid> getGrid(uint level) override;

    VF_PUBLIC void copyDataFromGpu();
    VF_PUBLIC virtual ~LevelGridBuilder();
    VF_PUBLIC void verifyGridNeighbors();

    VF_PUBLIC virtual void addGrid(real minX, real minY, real minZ, real maxX, real maxY, real maxZ, real delta,
        Device device, const std::string& distribution, bool periodictyX, bool periodictyY, bool periodictyZ);
    VF_PUBLIC virtual void generateGrids();

    VF_PUBLIC virtual void meshGeometry(std::string input, int level);
    VF_PUBLIC virtual void deleteSolidNodes();

    VF_PUBLIC virtual void writeGridToVTK(std::string output, int level);
    VF_PUBLIC virtual void writeSimulationFiles(std::string output, BoundingBox<int> &nodesDelete,
                                                bool writeFilesBinary, int level);

    VF_PUBLIC virtual std::shared_ptr<Grid> getGrid(int level, int box);

    VF_PUBLIC virtual void createBoundaryConditions();

    VF_PUBLIC virtual unsigned int getNumberOfNodes(unsigned int level) const;;
    VF_PUBLIC virtual std::vector<std::vector<std::vector<real> > > getQsValues() const;

    VF_PUBLIC virtual int getBoundaryConditionSize(int rb) const;
    VF_PUBLIC virtual std::vector<std::string> getTypeOfBoundaryConditions() const;

    VF_PUBLIC virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *nx,
                                         unsigned int *ny, unsigned int *nz, unsigned int *geo, const int level) const;
    VF_PUBLIC virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const;

    VF_PUBLIC virtual void setQs(real** q27, int* k, int channelSide, unsigned int level) const;
    VF_PUBLIC virtual void setOutflowValues(real* RhoBC, int* kN, int channelSide, int level) const;
    VF_PUBLIC virtual void setVelocityValues(real* vx, real* vy, real* vz, int channelSide, int level) const;
    VF_PUBLIC virtual void setPressValues(real* RhoBC, int* kN, int channelSide, int level) const;

    VF_PUBLIC void writeArrows(std::string fileName, std::shared_ptr<ArrowTransformator> trans) const;

protected:

    std::vector<std::shared_ptr<Grid> > grids;


    std::vector<std::vector<std::vector<real> > > Qs;
    std::vector<std::string> channelBoundaryConditions;

    void checkLevel(int level);

protected:
    void removeOverlapNodes();

    void createBCVectors();
    void addShortQsToVector(int index);
    void addQsToVector(int index);
    void fillRBForNode(int index, int direction, int directionSign, int rb);

    int getMatrixIndex(const int i) const;
    Vertex getVertex(const int matrixIndex) const;
    void writeArrow(const int i, const int qi, const Vertex& startNode,
                    std::shared_ptr<const ArrowTransformator> trans/*, std::shared_ptr<PolyDataWriterWrapper> writer*/)
    const;

private:

    Device device;
    std::string d3qxx;

public:
    VF_PUBLIC void getGridInformations(std::vector<int>& gridX, std::vector<int>& gridY,
                                       std::vector<int>& gridZ, std::vector<int>& distX, std::vector<int>& distY,
                                       std::vector<int>& distZ) override;
    VF_PUBLIC uint getNumberOfGridLevels() override;

    VF_PUBLIC uint getNumberOfNodesCF(int level) override;
    VF_PUBLIC uint getNumberOfNodesFC(int level) override;

    VF_PUBLIC void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const;

    VF_PUBLIC uint* getCF_coarse(uint level) const override;
    VF_PUBLIC uint* getCF_fine(uint level) const override;
    VF_PUBLIC uint* getFC_coarse(uint level) const override;
    VF_PUBLIC uint* getFC_fine(uint level) const override;

    VF_PUBLIC void setOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) override;
    VF_PUBLIC void setOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) override;
    
};

#endif

