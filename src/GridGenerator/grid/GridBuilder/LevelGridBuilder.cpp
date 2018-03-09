#include "LevelGridBuilder.h"

#include <stdio.h>
#include <iostream>
#include <GridGenerator/grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h>
#include <GridGenerator/grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h>
#include <GridGenerator/grid/partition/Partition.h>

#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>


#include <GridGenerator/io/GridVTKWriter/GridVTKWriter.h>
#include <GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h>
#include <GridGenerator/io/VTKWriterWrapper/UnstructuredGridWrapper.h>
#include <GridGenerator/io/VTKWriterWrapper/PolyDataWriterWrapper.h>

#include <GridGenerator/grid/NodeValues.h>

#include <GridGenerator/geometries/Arrow/ArrowImp.h>
#include <GridGenerator/utilities/Transformator/ArrowTransformator.h>

#include <utilities/logger/Logger.h>


#include <GridGenerator/geometries/Geometry/Serialization/GeometryMemento.h>

#include <GridGenerator/grid/GridFactory.h>
#include "grid/GridInterface.h"
//#include "grid/GridMocks.h"
#include "grid/Grid.h"

#define GEOFLUID 19
#define GEOSOLID 16

LevelGridBuilder::LevelGridBuilder(Device device, const std::string& d3qxx) : device(device), d3qxx(d3qxx)
{
    this->Qs.resize(QFILES);
    this->channelBoundaryConditions.resize(6);
    channelBoundaryConditions[0] = "periodic";
    channelBoundaryConditions[1] = "periodic";
    channelBoundaryConditions[2] = "periodic";
    channelBoundaryConditions[3] = "periodic";
    channelBoundaryConditions[4] = "periodic";
    channelBoundaryConditions[5] = "periodic";
}
//
//void LevelGridBuilder::setGrids(std::vector<SPtr<GridStub> > grids)
//{
//    auto gridFactory = SPtr<GridFactory>(new GridFactory());
//    gridFactory->setGridStrategy(device);
//
//    for (int i = int(grids.size()) - 1; i >= 0; i--)
//    {
//        const auto grid = gridFactory->makeGrid(grids[i]->startX, grids[i]->startY, grids[i]->startZ, grids[i]->endX, grids[i]->endY, grids[i]->endZ, grids[i]->delta, d3qxx);
//        this->grids.insert(this->grids.begin(), grid);
//        if(i==0)
//            this->grids[0]->setPeriodicity(true, true, true);
//        else
//            grid->setPeriodicity(false, false, false);
//        this->findGridInterface();
//        this->grids[0]->print();
//    }
//}

std::shared_ptr<LevelGridBuilder> LevelGridBuilder::makeShared(Device device, const std::string& d3qxx)
{
    return SPtr<LevelGridBuilder>(new LevelGridBuilder(device, d3qxx));
}


void LevelGridBuilder::copyDataFromGpu()
{
    for (const auto grid : grids)
    {
        auto gridGpuStrategy = std::dynamic_pointer_cast<GridGpuStrategy>(grid->getGridStrategy());
        if(gridGpuStrategy)
            gridGpuStrategy->copyDataFromGPU(std::static_pointer_cast<GridImp>(grid));
    }
        
}

LevelGridBuilder::~LevelGridBuilder()
{
    for (const auto grid : grids)
        grid->freeMemory();
}

void LevelGridBuilder::verifyGridNeighbors()
{
    //for (const auto grid : grids)
    //    std::cout << grid->verifyNeighborIndices();
}

void LevelGridBuilder::addGrid(real minX, real minY, real minZ, real maxX, real maxY, real maxZ, bool periodictyX, bool periodictyY, bool periodictyZ)
{
    auto gridFactory = SPtr<GridFactory>(new GridFactory());
    gridFactory->setGridStrategy(device);

    const auto grid = gridFactory->makeGrid(new Cuboid(minX, minY, minZ, maxX, maxY, maxZ), -1.0, d3qxx);
    grid->setPeriodicity(periodictyX, periodictyY, periodictyZ);

    grids.insert(grids.begin(), grid);
}

SPtr<Grid> LevelGridBuilder::getGrid(uint level)
{
    return grids[level];
}

void LevelGridBuilder::generateGrids()
{

}


void LevelGridBuilder::addGrid(real minX, real minY, real minZ, real maxX, real maxY, real maxZ, real delta, Device device, const std::string& distribution, bool periodictyX, bool periodictyY, bool periodictyZ)
{
    auto gridFactory = SPtr<GridFactory>(new GridFactory());
    gridFactory->setGridStrategy(device);

    const auto grid = gridFactory->makeGrid(new Cuboid(minX, minY, minZ, maxX, maxY, maxZ), delta, distribution);
    grids.insert(grids.begin(), grid);

    grid->setPeriodicity(periodictyX, periodictyY, periodictyZ);

    this->removeOverlapNodes();
}


void LevelGridBuilder::getGridInformations(std::vector<int>& gridX, std::vector<int>& gridY,
    std::vector<int>& gridZ, std::vector<int>& distX, std::vector<int>& distY,
    std::vector<int>& distZ)
{
    for (const auto grid : grids)
    {
        gridX.push_back(int(grid->getNumberOfNodesX()));
        gridY.push_back(int(grid->getNumberOfNodesY()));
        gridZ.push_back(int(grid->getNumberOfNodesZ()));

        distX.push_back(int(grid->getStartX()));
        distY.push_back(int(grid->getStartY()));
        distZ.push_back(int(grid->getStartZ()));
    }
}

void LevelGridBuilder::removeOverlapNodes()
{
    const uint numberOfLevels = getNumberOfGridLevels();
    if(numberOfLevels > 1)
    {
        grids[0]->findGridInterface(grids[1]);
    }   
}


void LevelGridBuilder::meshGeometry(std::string input, int level)
{
    checkLevel(level);

    Geometry geometry(input);

    if (geometry.size > 0)
        this->grids[level]->mesh(geometry);

}

uint LevelGridBuilder::getNumberOfGridLevels()
{
    return uint(grids.size());
}

uint LevelGridBuilder::getNumberOfNodesCF(int level)
{
    return this->grids[level]->getNumberOfNodesCF();
}

uint LevelGridBuilder::getNumberOfNodesFC(int level)
{
    return this->grids[level]->getNumberOfNodesFC();
}

uint* LevelGridBuilder::getCF_coarse(uint level) const
{
    return this->grids[level]->getCF_coarse();
}

uint* LevelGridBuilder::getCF_fine(uint level) const
{
    return this->grids[level]->getCF_fine();
}

uint* LevelGridBuilder::getFC_coarse(uint level) const
{
    return this->grids[level]->getFC_coarse();
}

uint* LevelGridBuilder::getFC_fine(uint level) const
{
    return this->grids[level]->getFC_fine();
}

void LevelGridBuilder::getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const
{
    this->grids[level]->getGridInterfaceIndices(iCellCfc, iCellCff, iCellFcc, iCellFcf);
}

void LevelGridBuilder::setOffsetFC(real * xOffCf, real * yOffCf, real * zOffCf, int level)
{
    for (uint i = 0; i < getNumberOfNodesFC(level); i++)
    {
        xOffCf[i] = 0.0;
        yOffCf[i] = 0.0;
        zOffCf[i] = 0.0;
    }
}

void LevelGridBuilder::setOffsetCF(real * xOffFc, real * yOffFc, real * zOffFc, int level)
{
    for (uint i = 0; i < getNumberOfNodesCF(level); i++)
    {
        xOffFc[i] = 0.0;
        yOffFc[i] = 0.0;
        zOffFc[i] = 0.0;
    }
}


void LevelGridBuilder::deleteSolidNodes()
{
    //this->gridWrapper[0]->deleteSolidNodes();
    //this->gridWrapper[0]->copyDataFromGPU();
}



uint LevelGridBuilder::getNumberOfNodes(unsigned int level) const
{
    return grids[level]->getSparseSize();
}

std::vector<std::vector<std::vector<real> > > LevelGridBuilder::getQsValues() const
{
    return this->Qs;
}

int LevelGridBuilder::getBoundaryConditionSize(int rb) const
{
    return (int)Qs[rb].size();
}

std::vector<std::string> LevelGridBuilder::getTypeOfBoundaryConditions() const
{
    return this->channelBoundaryConditions;
}

void LevelGridBuilder::writeGridToVTK(std::string output, int level)
{
   checkLevel(level);
   GridVTKWriter::writeGridToVTKXML(grids[level], output);
   GridVTKWriter::writeSparseGridToVTK(grids[level], output);
}


void LevelGridBuilder::writeSimulationFiles(std::string output, BoundingBox<int> &nodesDelete, bool writeFilesBinary, int level)
{
    //checkLevel(level);
    //UnstructuredLevelGridBuilder builder;
    //builder.buildUnstructuredGrid(this->gridKernels[level]->grid, nodesDelete);

    //std::vector<Node> coords = builder.getCoordsVec();
    //std::vector<std::vector<std::vector<real> > > qs = builder.getQsValues();
    //SimulationFileWriter::writeSimulationFiles(output, coords, qs, writeFilesBinary, this->gridKernels[level]->grid, this->transformators[level]);
}

std::shared_ptr<Grid> LevelGridBuilder::getGrid(int level, int box)
{
    return this->grids[level];
}

void LevelGridBuilder::checkLevel(int level)
{
    if (level >= grids.size()) { std::cout << "wrong level input... return to caller\n"; return; }
}


void LevelGridBuilder::getDimensions(int &nx, int &ny, int &nz, const int level) const
{
    nx = grids[level]->getNumberOfNodesX();
    ny = grids[level]->getNumberOfNodesY();
    nz = grids[level]->getNumberOfNodesZ();
}

void LevelGridBuilder::getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *neighborX, unsigned int *neighborY, unsigned int *neighborZ, unsigned int *geo, const int level) const
{
    grids[level]->getNodeValues(xCoords, yCoords, zCoords, neighborX, neighborY, neighborZ, geo);
}

void LevelGridBuilder::setQs(real** q27, int* k, int channelSide, unsigned int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        k[index] = (int)Qs[channelSide][index][0];
        for (int column = 1; column < Qs[channelSide][index].size(); column++) {
            q27[column - 1][index] = Qs[channelSide][index][column];

        }
    }
}

void LevelGridBuilder::setOutflowValues(real* RhoBC, int* kN, int channelSide, int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        RhoBC[index] = 0.0;
        kN[index] = 0;
    }
}

void LevelGridBuilder::setVelocityValues(real* vx, real* vy, real* vz, int channelSide, int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        vx[index] = 0.0;
        vy[index] = 0.0;
        vz[index] = 0.0;
    }
}

void LevelGridBuilder::setPressValues(real* RhoBC, int* kN, int channelSide, int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        RhoBC[index] = 0.0;
        kN[index] = 0;
    }
}


void LevelGridBuilder::createBoundaryConditions()
{
    this->createBCVectors();
}

/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void LevelGridBuilder::createBCVectors()
{
    Grid* grid = grids[0].get();
    for (uint i = 0; i < grid->getSize(); i++)
    {
        real x, y, z;
        grid->transIndexToCoords(grid->getSparseIndex(i), x, y, z);

        if (grid->getFieldEntry(grid->getSparseIndex(i)) == Q) /*addShortQsToVector(i);*/ addQsToVector(i);
        if (x == 0 && y < grid->getNumberOfNodesY() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 0, -1, INLETQS);
        if (x == grid->getNumberOfNodesX() - 2 && y < grid->getNumberOfNodesY() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 0, 1, OUTLETQS);

        if (z == grid->getNumberOfNodesZ() - 2 && x < grid->getNumberOfNodesX() - 1 && y < grid->getNumberOfNodesY() - 1) fillRBForNode(i, 2, 1, TOPQS);
        if (z == 0 && x < grid->getNumberOfNodesX() - 1 && y < grid->getNumberOfNodesY() - 1) fillRBForNode(i, 2, -1, BOTTOMQS);

        if (y == 0 && x < grid->getNumberOfNodesX() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 1, -1, FRONTQS);
        if (y == grid->getNumberOfNodesY() - 2 && x < grid->getNumberOfNodesX() - 1 && z < grid->getNumberOfNodesZ() - 1) fillRBForNode(i, 1, 1, BACKQS);
    }
}

void LevelGridBuilder::addShortQsToVector(int index)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    Grid* grid = grids[0].get();

    for (int i = grid->getEndDirection(); i >= 0; i--)
    {
        int qIndex = i * grid->getSize() + grid->getSparseIndex(index);
        real q = grid->getDistribution()[qIndex];
        if (q > 0) {
            //printf("Q%d (old:%d, new:%d), : %2.8f \n", i, coordsVec[index].matrixIndex, index, grid.d.f[i * grid.size + coordsVec[index].matrixIndex]);
            qKey += (uint32_t)pow(2, 26 - i);
            qNode.push_back(q);
        }
    }
    if (qKey > 0) {
        real transportKey = *((real*)&qKey);
        qNode.push_back(transportKey);
        qNode.push_back((real)index);
        Qs[GEOMQS].push_back(qNode);
    }
    qNode.clear();
}

void LevelGridBuilder::addQsToVector(int index)
{
    std::vector<real> qNode;
    qNode.push_back((real)index);

    Grid* grid = grids[0].get();

    for (int i = grid->getEndDirection(); i >= 0; i--)
    {
        int qIndex = i * grid->getSize() + grid->getSparseIndex(index);
        real q = grid->getDistribution()[qIndex];
        if (q > 0)
            qNode.push_back(q);
        else
            qNode.push_back(-1);
    }
    Qs[GEOMQS].push_back(qNode);
    qNode.clear();
}

void LevelGridBuilder::fillRBForNode(int index, int direction, int directionSign, int rb)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    Grid* grid = grids[0].get();

    for (int i = grid->getEndDirection(); i >= 0; i--)
    {
        if (grid->getDirection()[i * DIMENSION + direction] != directionSign)
            continue;

        qKey += (uint32_t)pow(2, 26 - i);
        qNode.push_back(0.5f);
    }
    if (qKey > 0) {
        real transportKey = *((real*)&qKey);
        qNode.push_back(transportKey);
        qNode.push_back((real)index);
        Qs[rb].push_back(qNode);
    }
    qNode.clear();
}

void LevelGridBuilder::writeArrows(std::string fileName, std::shared_ptr<ArrowTransformator> trans) const
{
    Grid* grid = grids[0].get();

    //std::shared_ptr<PolyDataWriterWrapper> writer = std::shared_ptr<PolyDataWriterWrapper>(new PolyDataWriterWrapper());
    for (int index = 0; index < Qs[GEOMQS].size(); index++)
    {
        Vertex startNode = getVertex(getMatrixIndex(index));
        //for (int qi = grid.d.dir_start; qi <= grid.d.dir_end; qi++)
            //writeArrow(index, qi, startNode, trans, writer);
    }
    //writer->writePolyDataToFile(fileName);
}

void LevelGridBuilder::writeArrow(const int i, const int qi, const Vertex& startNode, std::shared_ptr<const ArrowTransformator> trans/*, std::shared_ptr<PolyDataWriterWrapper> writer*/) const
{
    Grid* grid = grids[0].get();

    real qval = Qs[GEOMQS][i][qi + 1];
    if (qval > 0.0f)
    {
        int qReverse = grid->getEndDirection() - qi;
        Vertex dir((real)grid->getDirection()[qReverse * DIMENSION + 0], (real)grid->getDirection()[qReverse* DIMENSION + 1], (real)grid->getDirection()[qReverse * DIMENSION + 2]);
        Vertex nodeOnGeometry(startNode + (dir * qval));
        std::shared_ptr<Arrow> arrow = ArrowImp::make(startNode, nodeOnGeometry);
        trans->transformGridToWorld(arrow);
        //writer->addVectorArrow(arrow);
    }
}



Vertex LevelGridBuilder::getVertex(int matrixIndex) const
{
    real x, y, z;
    this->grids[0]->transIndexToCoords(matrixIndex, x, y, z);
    return Vertex(x,y,z);
}

int LevelGridBuilder::getMatrixIndex(int i) const
{
    int index = (int)Qs[GEOMQS][i][0];
    return this->grids[0]->getSparseIndex(index);
}

