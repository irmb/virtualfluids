#ifndef GridBuilder_H
#define GridBuilder_H

#include "GridGenerator/global.h"


#include <vector>
#include <string>
#include <memory>


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

template <typename T>
class BoundingBox;

struct Grid;

class GridBuilder
{
public:
    enum class GenerationDevice
    {
        CPU, GPU
    };

    virtual VF_PUBLIC ~GridBuilder() {};

    virtual void meshGeometry(std::string input, int level) = 0;
    virtual void deleteSolidNodes() = 0;

	virtual void flood(Vertex &startFlood, int level) = 0;

	virtual void writeGridToVTK(std::string output, int level) = 0;
    virtual void writeSimulationFiles(std::string output, BoundingBox<int> &nodesDelete, bool writeFilesBinary, int level) = 0;
    virtual void writeArrows(std::string fileName, std::shared_ptr<ArrowTransformator> trans) const = 0;

	virtual std::shared_ptr<Grid> getGrid(int level, int box) = 0;

    virtual void createBoundaryConditions() = 0;
    virtual std::vector<std::vector<std::vector<real> > > getQsValues() const = 0;
    virtual int getBoundaryConditionSize(int rb) const = 0;
    virtual std::vector<std::string> getTypeOfBoundaryConditions() const = 0;
    virtual unsigned int getNumberOfNodes(unsigned int level) const = 0;
    virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *nx, unsigned int *ny, unsigned int *nz, unsigned int *geo, const int level) const = 0;
    virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const = 0;
    virtual void setQs(real** q27, int* k, int channelSide, unsigned int level) const = 0;
    virtual void setOutflowValues(real* RhoBC, int* kN, int channelSide, int level) const = 0;
    virtual void setVelocityValues(real* vx, real* vy, real* vz, int channelSide, int level) const = 0;
    virtual void setPressValues(real* RhoBC, int* kN, int channelSide, int level) const = 0;

};

#endif

