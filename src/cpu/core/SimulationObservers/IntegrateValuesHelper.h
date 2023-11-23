#ifndef INTEGRATEVALUESHELPER_H
#define INTEGRATEVALUESHELPER_H

#include <PointerDefinitions.h>

#include "Block3D.h"
#include "CbArray2D.h"
#include <parallel/Communicator.h>
#include "D3Q27System.h"
#include "GbCuboid3D.h"
#include "Grid3D.h"

// struct CalcNodes
//{
//    SPtr<Block3D> block;
//    std::vector<UbTupleInt3> nodes;
//};
//
// struct Nodes
//{
//   SPtr<Block3D> block;
//   UbTupleInt3 nodes;
//};

class IntegrateValuesHelper
{
public:
    struct CalcNodes {
        SPtr<Block3D> block;
        std::vector<UbTupleInt3> nodes;
    };

    struct Node {
        SPtr<Block3D> block;
        UbTupleInt3 node;
    };

public:
    IntegrateValuesHelper(SPtr<Grid3D> grid, std::shared_ptr<vf::parallel::Communicator> comm, real minX1, real minX2, real minX3,
                          real maxX1, real maxX2, real maxX3);
    IntegrateValuesHelper(SPtr<Grid3D> grid, std::shared_ptr<vf::parallel::Communicator> comm, real minX1, real minX2, real minX3,
                          real maxX1, real maxX2, real maxX3, int level);
    virtual ~IntegrateValuesHelper();

    void calculateMQ();
    void calculateAV();
    void clearData();

    real getRho() { return sRho; }
    real getVx1() { return sVx1; }
    real getVx2() { return sVx2; }
    real getVx3() { return sVx3; }
    real getCellsVolume() { return sCellVolume; }
    //  real getVm() { return sVm; }
    // real getPress() {return sPress;}
    real getAvVx1() { return sAvVx1; }
    real getAvVx2() { return sAvVx2; }
    real getAvVx3() { return sAvVx3; }
    real getTSx1() { return sTSx1; }
    real getTSx2() { return sTSx2; }
    real getTSx3() { return sTSx3; }
    real getTSx1x3() { return sTSx1x3; }

    real getNumberOfFluidsNodes();
    real getNumberOfSolidNodes();
    GbCuboid3DPtr getBoundingBox();
    std::vector<CalcNodes> getCNodes();

protected:
private:
    void init(int level);

    bool root;
    SPtr<Grid3D> grid;
    real sVx1, sVx2, sVx3, sRho, sCellVolume; // sPress, sVm;
    real numberOfFluidsNodes, numberOfSolidNodes;
    real sAvVx1, sAvVx2, sAvVx3, sTSx1, sTSx2, sTSx3, sTSx1x3;
    std::vector<CalcNodes> cnodes;
    GbCuboid3DPtr boundingBox;
    std::shared_ptr<vf::parallel::Communicator> comm;
    CbArray2D<Node> cnodes2DMatrix;
    enum Values { AvVx = 0, AvVy = 1, AvVz = 2, AvVxx = 3, AvVyy = 4, AvVzz = 5, AvVxy = 6, AvVyz = 7, AvVxz = 8 };
};

#endif
