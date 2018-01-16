#ifndef LEVEL_GRID_BUILDER_H
#define LEVEL_GRID_BUILDER_H

#include "GridGenerator/global.h"


#include <vector>
#include <string>
#include <memory>

#include "GridBuilder.h"

struct Vertex;
struct Grid;
class Transformator;
class ArrowTransformator;
class PolyDataWriterWrapper;

template <typename T>
class BoundingBox;


class LevelGridBuilder : public GridBuilder
{
public:

    VF_PUBLIC LevelGridBuilder(GenerationDevice device);
    VF_PUBLIC static std::shared_ptr<GridBuilder> make(std::string);

    VF_PUBLIC virtual ~LevelGridBuilder();

    VF_PUBLIC virtual void addGrid(uint minX, uint minY, uint minZ, uint maxX, uint maxY, uint maxZ, std::string distribution);

	VF_PUBLIC virtual void meshGeometry(std::string input, int level);
    VF_PUBLIC virtual void deleteSolidNodes();

	VF_PUBLIC virtual void flood(Vertex &startFlood, int level);

	VF_PUBLIC virtual void writeGridToVTK(std::string output, int level);
	VF_PUBLIC virtual void writeSimulationFiles(std::string output, BoundingBox<int> &nodesDelete, bool writeFilesBinary, int level);

	VF_PUBLIC virtual std::shared_ptr<Grid> getGrid(int level, int box);

    VF_PUBLIC virtual void createBoundaryConditions();

    VF_PUBLIC virtual unsigned int getNumberOfNodes(unsigned int level) const;;
    VF_PUBLIC virtual std::vector<std::vector<std::vector<real> > > getQsValues() const;

    VF_PUBLIC virtual int getBoundaryConditionSize(int rb) const;
    VF_PUBLIC virtual std::vector<std::string> getTypeOfBoundaryConditions() const;

    VF_PUBLIC virtual void getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *nx, unsigned int *ny, unsigned int *nz, unsigned int *geo, const int level) const;
    VF_PUBLIC virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const;

    VF_PUBLIC virtual void setQs(real** q27, int* k, int channelSide, unsigned int level) const;
    VF_PUBLIC virtual void setOutflowValues(real* RhoBC, int* kN, int channelSide, int level) const;
    VF_PUBLIC virtual void setVelocityValues(real* vx, real* vy, real* vz, int channelSide, int level) const;
    VF_PUBLIC virtual void setPressValues(real* RhoBC, int* kN, int channelSide, int level) const;

    VF_PUBLIC void writeArrows(std::string fileName, std::shared_ptr<ArrowTransformator> trans) const;

protected:
    GenerationDevice device;

    std::vector<std::shared_ptr<Grid> > grids;


    std::vector<std::vector<std::vector<real> > > Qs;
    std::vector<std::string> channelBoundaryConditions;

    void checkLevel(int level);

protected:
    void createBCVectors();
    void addShortQsToVector(int index);
    void addQsToVector(int index);
    void fillRBForNode(int x, int y, int z, int index, int direction, int directionSign, int rb);

    int getMatrixIndex(const int i) const;
    Vertex getVertex(const int matrixIndex) const;
    void writeArrow(const int i, const int qi, const Vertex& startNode, std::shared_ptr<const ArrowTransformator> trans/*, std::shared_ptr<PolyDataWriterWrapper> writer*/) const;

};

#endif

