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


#define GEOFLUID 19
#define GEOSOLID 16

LevelGridBuilder::LevelGridBuilder()
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


LevelGridBuilder::~LevelGridBuilder()
{
    for (auto grid : grids)
        grid->freeMemory();
}

void LevelGridBuilder::addGrid(real minX, real minY, real minZ, real maxX, real maxY, real maxZ, real delta, const std::string& device, const std::string& distribution)
{
    const auto grid = GridFactory::makeGrid(minX, minY, minZ, maxX, maxY, maxZ, delta, device, distribution);
    grid->print();
    grids.push_back(grid);
    this->removeOverlapNodes();
}


void LevelGridBuilder::removeOverlapNodes()
{
    for (int level = 0; level < grids.size() - 1; level++)
        grids[level]->removeOverlapNodes(grids[level + 1]);
}


void LevelGridBuilder::meshGeometry(std::string input, int level)
{
    checkLevel(level);

    Geometry geometry(input);

    if (geometry.size > 0)
        this->grids[level]->mesh(geometry);

}

void LevelGridBuilder::deleteSolidNodes()
{
    //this->gridWrapper[0]->deleteSolidNodes();
    //this->gridWrapper[0]->copyDataFromGPU();
}

void LevelGridBuilder::flood(Vertex &startFlood, int level)
{
    checkLevel(level);
    //this->gridWrapper[level]->floodFill(startFlood);
}

void LevelGridBuilder::createBoundaryConditions()
{
    this->createBCVectors();
}


unsigned int LevelGridBuilder::getNumberOfNodes(unsigned int level) const
{
    return (unsigned int) grids[level]->reducedSize;
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

    for (uint i = 0; i < grid.size; i++)
    {
        if (grid.matrixIndex[i] == -1)
            continue;

        real x, y, z;
        grid.transIndexToCoords(i, x, y, z);

        xCoords[i + 1] = x;
        yCoords[i + 1] = y;
        zCoords[i + 1] = z;
        neighborX[i + 1] = (unsigned int)(grid.neighborIndexX[i] + 1);
        neighborY[i + 1] = (unsigned int)(grid.neighborIndexY[i] + 1);
        neighborZ[i + 1] = (unsigned int)(grid.neighborIndexZ[i] + 1);
        geo[i + 1] = (unsigned int)grid.isSolid(i) ? GEOSOLID : GEOFLUID;
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



/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void LevelGridBuilder::createBCVectors()
{
    Grid grid = *grids[0].get();
    for (uint i = 0; i < grid.reducedSize; i++)
    {
        real x, y, z;
        grid.transIndexToCoords(grid.matrixIndex[i], x, y, z);

        if (grid.field[grid.matrixIndex[i]] == Q) /*addShortQsToVector(i);*/ addQsToVector(i);
        if (x == 0 && y < grid.ny - 1 && z < grid.nz - 1) fillRBForNode(x, y, z, i, 0, -1, INLETQS);
        if (x == grid.nx - 2 && y < grid.ny - 1 && z < grid.nz - 1) fillRBForNode(x, y, z, i, 0, 1, OUTLETQS);

        if (z == grid.nz - 2 && x < grid.nx - 1 && y < grid.ny - 1) fillRBForNode(x, y, z, i, 2, 1, TOPQS);
        if (z == 0 && x < grid.nx - 1 && y < grid.ny - 1) fillRBForNode(x, y, z, i, 2, -1, BOTTOMQS);

        if (y == 0 && x < grid.nx - 1 && z < grid.nz - 1) fillRBForNode(x, y, z, i, 1, -1, FRONTQS);
        if (y == grid.ny - 2 && x < grid.nx - 1 && z < grid.nz - 1) fillRBForNode(x, y, z, i, 1, 1, BACKQS);
    }
}

void LevelGridBuilder::addShortQsToVector(int index)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    Grid grid = *grids[0].get();

    for (int i = grid.d.dir_end; i >= 0; i--)
    {
        int qIndex = i * grid.size + grid.matrixIndex[index];
        real q = grid.d.f[qIndex];
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

    for (int i = grid.d.dir_end; i >= 0; i--)
    {
        int qIndex = i * grid.size + grid.matrixIndex[index];
        real q = grid.d.f[qIndex];
        if (q > 0)
            qNode.push_back(q);
        else
            qNode.push_back(-1);
    }
    Qs[GEOMQS].push_back(qNode);
    qNode.clear();
}

void LevelGridBuilder::fillRBForNode(int x, int y, int z, int index, int direction, int directionSign, int rb)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    Grid grid = *grids[0].get();

    for (int i = grid.d.dir_end; i >= 0; i--)
    {
        if (grid.d.dirs[i * DIMENSION + direction] != directionSign)
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
        int qReverse = grid.d.dir_end - qi;
        Vertex dir((real)grid.d.dirs[qReverse * DIMENSION + 0], (real)grid.d.dirs[qReverse* DIMENSION + 1], (real)grid.d.dirs[qReverse * DIMENSION + 2]);
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

