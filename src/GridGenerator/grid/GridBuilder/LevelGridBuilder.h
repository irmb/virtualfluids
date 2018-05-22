#ifndef LEVEL_GRID_BUILDER_H
#define LEVEL_GRID_BUILDER_H

#include "GridGenerator/global.h"

#include <vector>
#include <string>
#include <memory>
#include <map>

#include "GridBuilder.h"
#include "grid/GridInterface.h"

#include "grid/Grid.h"
#include "grid/NodeValues.h"

struct Vertex;
class  Grid;
class Transformator;
class ArrowTransformator;
class PolyDataWriterWrapper;
class BoundaryCondition;
class BoundingBox;
enum class Device;

class Side;

class BoundaryCondition
{
public:
    std::vector<uint> indices;

};

class PressureBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<PressureBoundaryCondition> make(real rho)
    {
        return SPtr<PressureBoundaryCondition>(new PressureBoundaryCondition(rho));
    }

    SPtr<Side> side;
    real rho;
private:
    PressureBoundaryCondition(real rho) : rho(rho)
    {

    }
};

class VelocityBoundaryCondition : public BoundaryCondition
{
public:
    static SPtr<VelocityBoundaryCondition> make(real vx, real vy, real vz)
    {
        return SPtr<VelocityBoundaryCondition>(new VelocityBoundaryCondition(vx, vy, vz));
    }

    SPtr<Side> side;
    real vx, vy, vz;
private:
    VelocityBoundaryCondition(real vx, real vy, real vz) : vx(vx), vy(vy), vz(vz)
    {
        
    } 
};

#define X_INDEX 0
#define Y_INDEX 1
#define Z_INDEX 2

#define POSITIVE_DIR 1
#define NEGATIVE_DIR -1

class Side
{
public:
    virtual void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition) = 0;
    virtual void setPeriodicy(SPtr<Grid> grid) = 0;
    
    virtual int getCoordinate() const = 0;
    virtual int getDirection() const = 0;

protected:
    static void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition, std::string coord, real constant, real startInner, real endInner, real startOuter, real endOuter) {
        for (int v1 = startInner; v1 < endInner; v1 += grid->getDelta())
        {
            for (int v2 = startOuter; v2 < endOuter; v2 += grid->getDelta())
            {
                uint index = getIndex(grid, coord, constant, v1, v2);
                if (grid->getFieldEntry(index) == FLUID)
                    boundaryCondition->indices.push_back(index);
            }
        }
    }

    static uint getIndex(SPtr<Grid> grid, std::string coord, real constant, real v1, real v2)
    {
        if (coord == "x")
            return grid->transCoordToIndex(constant, v1, v2);
        if (coord == "y")
            return grid->transCoordToIndex(v1, constant, v2);
        if (coord == "z")
            return grid->transCoordToIndex(v1, v2, constant);
    }
};

class MX : public Side
{
public:
    void setPeriodicy(SPtr<Grid> grid) override {
        grid->setPeriodicityX(false);
    }

    void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition) override {
        Side::addIndices(grid, boundaryCondition, "x", grid->getStartX() + grid->getDelta(),  grid->getStartY(), grid->getEndY(),grid->getStartZ(), grid->getEndZ());
    }

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }
};

class PX : public Side
{
public:
    void setPeriodicy(SPtr<Grid> grid) override {
        grid->setPeriodicityX(false);
    }

    void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition) override {
        Side::addIndices(grid, boundaryCondition, "x", grid->getEndX() - grid->getDelta(),  grid->getStartY(), grid->getEndY(), grid->getStartZ(), grid->getEndZ());
    }

    int getCoordinate() const override
    {
        return X_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }
};


class MY : public Side
{
public:
    void setPeriodicy(SPtr<Grid> grid) override {
        grid->setPeriodicityY(false);
    }

    void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition) override {
        Side::addIndices(grid, boundaryCondition, "y", grid->getStartY() + grid->getDelta(), grid->getStartX(), grid->getEndX(), grid->getStartZ(), grid->getEndZ());
    }

    int getCoordinate() const override
    {
        return Y_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }
};

class PY : public Side
{
public:
    void setPeriodicy(SPtr<Grid> grid) override {
        grid->setPeriodicityY(false);
    }

    void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition) override {
        Side::addIndices(grid, boundaryCondition, "y", grid->getEndY() - grid->getDelta(), grid->getStartX(), grid->getEndX(), grid->getStartZ(), grid->getEndZ());
    }

    int getCoordinate() const override
    {
        return Y_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }
};


class MZ : public Side
{
public:
    void setPeriodicy(SPtr<Grid> grid) override {
        grid->setPeriodicityZ(false);
    }

    void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition) override {
        Side::addIndices(grid, boundaryCondition, "z", grid->getStartZ() + grid->getDelta(), grid->getStartX(), grid->getEndX(), grid->getStartY(), grid->getEndY());
    }

    int getCoordinate() const override
    {
        return Z_INDEX;
    }

    int getDirection() const override
    {
        return NEGATIVE_DIR;
    }
};

class PZ : public Side
{
public:
    void setPeriodicy(SPtr<Grid> grid) override {
        grid->setPeriodicityZ(false);
    }

    void addIndices(SPtr<Grid> grid, SPtr<BoundaryCondition> boundaryCondition) override {
        Side::addIndices(grid, boundaryCondition, "z", grid->getEndZ() - grid->getDelta(), grid->getStartX(), grid->getEndX(), grid->getStartY(), grid->getEndY());
    }

    int getCoordinate() const override
    {
        return Z_INDEX;
    }

    int getDirection() const override
    {
        return POSITIVE_DIR;
    }
};


class LevelGridBuilder : public GridBuilder
{
protected:
    VF_PUBLIC LevelGridBuilder(Device device, const std::string& d3qxx);

public:
    VF_PUBLIC static std::shared_ptr<LevelGridBuilder> makeShared(Device device, const std::string& d3qxx);

    VF_PUBLIC SPtr<Grid> getGrid(uint level) override;

    VF_PUBLIC void copyDataFromGpu();
    VF_PUBLIC virtual ~LevelGridBuilder();

    VF_PUBLIC void setVelocityBoundaryCondition(SPtr<Side> side, real vx, real vy, real vz);
    VF_PUBLIC void setPressureBoundaryCondition(SPtr<Side> side, real rho);

    VF_PUBLIC virtual std::shared_ptr<Grid> getGrid(int level, int box);

    VF_PUBLIC virtual void createBoundaryConditions();

    VF_PUBLIC virtual unsigned int getNumberOfNodes(unsigned int level) const;;
    VF_PUBLIC virtual std::vector<std::vector<std::vector<real> > > getQsValues() const;

    VF_PUBLIC virtual int getBoundaryConditionSize(int rb) const;
    VF_PUBLIC virtual std::vector<std::string> getTypeOfBoundaryConditions() const;

    //VF_PUBLIC virtual void setInflowBoundaryCondition(BoundaryCondition boundaryCondition);
    //VF_PUBLIC virtual void setOutflowBoundaryCondition(BoundaryCondition boundaryCondition);
    //VF_PUBLIC virtual std::vector<BoundaryCondition> getTypeOfBoundaryCondition() const;


    VF_PUBLIC virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *nx,
                                         unsigned int *ny, unsigned int *nz, unsigned int *geo, const int level) const;
    VF_PUBLIC virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const;

    VF_PUBLIC virtual void setQs(real** q27, int* k, int channelSide, unsigned int level) const;
    VF_PUBLIC virtual void setOutflowValues(real* RhoBC, int* kN, int channelSide, int level) const;
    VF_PUBLIC virtual void setVelocityValues(real* vx, real* vy, real* vz, int channelSide, int level) const;

    VF_PUBLIC uint getVelocitySize(int level) const;
    VF_PUBLIC virtual void getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const;
    VF_PUBLIC virtual void getVelocityQs(real* qs[27], int level) const;
    VF_PUBLIC uint getPressureSize(int level) const override;
    VF_PUBLIC void getPressureValues(real* rho, int* indices, int level) const override;

    VF_PUBLIC virtual void setPressValues(real* RhoBC, int* kN, int channelSide, int level) const;

    VF_PUBLIC void writeArrows(std::string fileName, std::shared_ptr<ArrowTransformator> trans) const;

protected:

    std::vector<std::shared_ptr<Grid> > grids;
    std::vector<std::vector<std::vector<real> > > Qs;
    std::vector<std::string> channelBoundaryConditions;
    std::vector<SPtr<VelocityBoundaryCondition> > velocityBoundaryConditions;
    std::vector<SPtr<PressureBoundaryCondition> > pressureBoundaryConditions;

    //std::map<Side, BoundaryCondition> channelBoundaryConditionTypes;

    void checkLevel(int level);

protected:

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

    VF_PUBLIC void getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const override;

    VF_PUBLIC void setOffsetFC(real* xOffCf, real* yOffCf, real* zOffCf, int level) override;
    VF_PUBLIC void setOffsetCF(real* xOffFc, real* yOffFc, real* zOffFc, int level) override;

};

#endif

