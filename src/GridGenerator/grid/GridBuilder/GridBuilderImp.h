#ifndef GridBuilderImp_H
#define GridBuilderImp_H

#include "GridGenerator/global.h"


#include <vector>
#include <string>
#include <memory>

#include "GridBuilder.h"

struct Vertex;
class GridWrapper;
class Transformator;
class ArrowTransformator;
class PolyDataWriterWrapper;

template <typename T>
class BoundingBox;


class GridBuilderImp : public GridBuilder
{
public:
    VF_PUBLIC static std::shared_ptr<GridBuilder> make(std::string);

    VF_PUBLIC virtual ~GridBuilderImp();

    VF_PUBLIC virtual void addGrid(doubflo length, doubflo width, doubflo high, doubflo delta, std::string distribution, std::shared_ptr<Transformator> trans);

	VF_PUBLIC virtual void meshGeometry(std::string input, int level);
    VF_PUBLIC virtual void deleteSolidNodes();

	VF_PUBLIC virtual void flood(Vertex &startFlood, int level);

	VF_PUBLIC virtual void writeGridToVTK(std::string output, int level);
	VF_PUBLIC virtual void writeSimulationFiles(std::string output, BoundingBox<int> &nodesDelete, bool writeFilesBinary, int level);

	VF_PUBLIC virtual std::shared_ptr<GridWrapper> getKernel(int level, int box);

    VF_PUBLIC virtual void createBoundaryConditions();

    VF_PUBLIC virtual unsigned int getNumberOfNodes(unsigned int level) const;;
    VF_PUBLIC virtual std::vector<std::vector<std::vector<doubflo> > > getQsValues() const;

    VF_PUBLIC virtual int getBoundaryConditionSize(int rb) const;
    VF_PUBLIC virtual std::vector<std::string> getTypeOfBoundaryConditions() const;

    VF_PUBLIC virtual void getNodeValues(doubflo *xCoords, doubflo *yCoords, doubflo *zCoords, unsigned int *nx, unsigned int *ny, unsigned int *nz, unsigned int *geo, const int level) const;
    VF_PUBLIC virtual void getDimensions(int &nx, int &ny, int &nz, const int level) const;

    VF_PUBLIC virtual void setQs(doubflo** q27, int* k, int channelSide, unsigned int level) const;
    VF_PUBLIC virtual void setOutflowValues(doubflo* RhoBC, int* kN, int channelSide, int level) const;
    VF_PUBLIC virtual void setVelocityValues(doubflo* vx, doubflo* vy, doubflo* vz, int channelSide, int level) const;
    VF_PUBLIC virtual void setPressValues(doubflo* RhoBC, int* kN, int channelSide, int level) const;

    VF_PUBLIC void writeArrows(std::string fileName, std::shared_ptr<ArrowTransformator> trans) const;

protected:
    GridBuilderImp();

    std::vector<std::vector<std::shared_ptr<GridWrapper> > >gridKernels;
    std::vector<std::shared_ptr<Transformator>> transformators;
    std::vector<std::vector<BoundingBox<int>> > boxes;
    std::vector<int> rankTasks;
    std::vector<std::vector<int> > gridDimensions;


    std::vector<std::vector<std::vector<doubflo> > > Qs;
    std::vector<std::string> channelBoundaryConditions;

    void checkLevel(int level);

protected:
    virtual void createGridKernels(std::string distribution) = 0;

    void setNumberOfNodes(doubflo length, doubflo width, doubflo high, doubflo delta);
    void printMasterInformation(int nx, int ny, int nz);
    void setCudaDevice(int rank);
    void rebuildBoxes();
    void sendTasks();
    void receiveTasks();
    void writeBoxes(std::string name);

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

