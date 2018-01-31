#include "LevelGridBuilder.h"

#include <stdio.h>
#include <iostream>
#include <GridGenerator/grid/Grid.cuh>
#include <GridGenerator/grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h>
#include <GridGenerator/grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h>
#include <GridGenerator/grid/partition/Partition.h>

#include <GridGenerator/geometries/Triangle/Triangle.cuh>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.cuh>
#include <GridGenerator/geometries/Geometry/Geometry.cuh>

#include <GridGenerator/utilities/Transformator/Transformator.h>

#include <GridGenerator/io/GridVTKWriter/GridVTKWriter.h>
#include <GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h>
#include <GridGenerator/io/VTKWriterWrapper/UnstructuredGridWrapper.h>
#include <GridGenerator/io/VTKWriterWrapper/PolyDataWriterWrapper.h>

#include <GridGenerator/grid/NodeValues.h>

#include <GridGenerator/geometries/Arrow/ArrowImp.h>
#include <GridGenerator/utilities/Transformator/ArrowTransformator.h>

#include <utilities/logger/Logger.h>

#include <utilities/StringUtil/StringUtil.h>

#include <GridGenerator/geometries/Geometry/Serialization/GeometryMemento.h>

#include <GridGenerator/grid/GridFactory.h>
#include "grid/GridInterface.cuh"


#define GEOFLUID 19
#define GEOSOLID 16

LevelGridBuilder::LevelGridBuilder(const std::string& device, const std::string& d3qxx) : device(device), d3qxx(d3qxx)
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

std::shared_ptr<LevelGridBuilder> LevelGridBuilder::makeShared(const std::string& device, const std::string& d3qxx)
{
    return SPtr<LevelGridBuilder>(new LevelGridBuilder(device, d3qxx));
}


void LevelGridBuilder::copyDataFromGpu()
{
    for (const auto grid : grids)
    {
        auto gridGpuStrategy = std::dynamic_pointer_cast<GridGpuStrategy>(grid->gridStrategy);
        if(gridGpuStrategy)
            gridGpuStrategy->copyDataFromGPU(grid);
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
    const auto grid = GridFactory::makeGrid(minX, minY, minZ, maxX, maxY, maxZ, -1.0, device, d3qxx);
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


void LevelGridBuilder::addGrid(real minX, real minY, real minZ, real maxX, real maxY, real maxZ, real delta, const std::string& device, const std::string& distribution, bool periodictyX, bool periodictyY, bool periodictyZ)
{
    const auto grid = GridFactory::makeGrid(minX, minY, minZ, maxX, maxY, maxZ, delta, device, distribution);
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
        gridX.push_back(int(grid->nx));
        gridY.push_back(int(grid->ny));
        gridZ.push_back(int(grid->nz));

        distX.push_back(int(grid->startX));
        distY.push_back(int(grid->startY));
        distZ.push_back(int(grid->startZ));
    }
}

void LevelGridBuilder::removeOverlapNodes()
{
    const uint numberOfLevels = getNumberOfGridLevels();
    if(numberOfLevels > 1)
    {
        grids[0]->removeOverlapNodes(grids[1]);
        grids[0]->print();
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
    return grids[level]->getReducedSize();
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
   GridVTKWriter::writeGridToVTKXML(*grids[level].get(), output);
   GridVTKWriter::writeSparseGridToVTK(*grids[level].get(), output);
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
    nx = grids[level]->nx;
    ny = grids[level]->ny;
    nz = grids[level]->nz;
}

void LevelGridBuilder::getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *neighborX, unsigned int *neighborY, unsigned int *neighborZ, unsigned int *geo, const int level) const
{
    xCoords[0] = 0;
    yCoords[0] = 0;
    zCoords[0] = 0;
    neighborX[0] = 0;
    neighborY[0] = 0;
    neighborZ[0] = 0;
    geo[0] = GEOSOLID;

    Grid grid = *grids[level].get();

    int nodeNumber = 0;
    for (uint i = 0; i < grid.getSize(); i++)
    {
        if (grid.matrixIndex[i] == -1)
            continue;

        real x, y, z;
        grid.transIndexToCoords(i, x, y, z);

        xCoords[nodeNumber + 1] = x;
        yCoords[nodeNumber + 1] = y;
        zCoords[nodeNumber + 1] = z;
        neighborX[nodeNumber + 1] = uint(grid.neighborIndexX[i] + 1);
        neighborY[nodeNumber + 1] = uint(grid.neighborIndexY[i] + 1);
        neighborZ[nodeNumber + 1] = uint(grid.neighborIndexZ[i] + 1);
        geo[nodeNumber + 1] = uint(grid.isSolid(i) ? GEOSOLID : GEOFLUID);
        nodeNumber++;
    }
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
    Grid grid = *grids[0].get();
    for (uint i = 0; i < grid.getSize(); i++)
    {
        real x, y, z;
        grid.transIndexToCoords(grid.matrixIndex[i], x, y, z);

        if (grid.field[grid.matrixIndex[i]] == Q) /*addShortQsToVector(i);*/ addQsToVector(i);
        if (x == 0 && y < grid.ny - 1 && z < grid.nz - 1) fillRBForNode(i, 0, -1, INLETQS);
        if (x == grid.nx - 2 && y < grid.ny - 1 && z < grid.nz - 1) fillRBForNode(i, 0, 1, OUTLETQS);

        if (z == grid.nz - 2 && x < grid.nx - 1 && y < grid.ny - 1) fillRBForNode(i, 2, 1, TOPQS);
        if (z == 0 && x < grid.nx - 1 && y < grid.ny - 1) fillRBForNode(i, 2, -1, BOTTOMQS);

        if (y == 0 && x < grid.nx - 1 && z < grid.nz - 1) fillRBForNode(i, 1, -1, FRONTQS);
        if (y == grid.ny - 2 && x < grid.nx - 1 && z < grid.nz - 1) fillRBForNode(i, 1, 1, BACKQS);
    }
}

void LevelGridBuilder::addShortQsToVector(int index)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    Grid grid = *grids[0].get();

    for (int i = grid.distribution.dir_end; i >= 0; i--)
    {
        int qIndex = i * grid.getSize() + grid.matrixIndex[index];
        real q = grid.distribution.f[qIndex];
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

    Grid grid = *grids[0].get();

    for (int i = grid.distribution.dir_end; i >= 0; i--)
    {
        int qIndex = i * grid.getSize() + grid.matrixIndex[index];
        real q = grid.distribution.f[qIndex];
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

    Grid grid = *grids[0].get();

    for (int i = grid.distribution.dir_end; i >= 0; i--)
    {
        if (grid.distribution.dirs[i * DIMENSION + direction] != directionSign)
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
    Grid grid = *grids[0].get();

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
    Grid grid = *grids[0].get();

    real qval = Qs[GEOMQS][i][qi + 1];
    if (qval > 0.0f)
    {
        int qReverse = grid.distribution.dir_end - qi;
        Vertex dir((real)grid.distribution.dirs[qReverse * DIMENSION + 0], (real)grid.distribution.dirs[qReverse* DIMENSION + 1], (real)grid.distribution.dirs[qReverse * DIMENSION + 2]);
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
    return this->grids[0]->matrixIndex[index];
}

