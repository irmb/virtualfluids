#include "GridBuilderImp.h"
#include "mpi.h"

#include <stdio.h>
#include <iostream>
#include <GridGenerator/grid/Grid.cuh>
#include <GridGenerator/grid/GridWrapper/GridWrapperCPU/GridWrapperCPU.h>
#include <GridGenerator/grid/GridWrapper/GridWrapperGPU/GridWrapperGPU.h>
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

#include <GridGenerator/grid/GridBuilder/GridCpuBuilder/GridCpuBuilder.h>
#include <GridGenerator/grid/GridBuilder/GridGpuBuilder/GridGpuBuilder.h>

#include <utilities/StringUtil/StringUtil.h>

#include <GridGenerator/geometries/Geometry/Serialization/GeometryMemento.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

void serialize(const GeometryMemento &memento, const std::string &filename)
{
    std::ofstream ofs(filename, std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << memento;
}
void deserialize(GeometryMemento &memento, const std::string &filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    ia >> memento;
}

#define GEOFLUID 19
#define GEOSOLID 16


GridBuilderImp::GridBuilderImp(GenerationDevice device) : device(device)
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


GridBuilderImp::~GridBuilderImp()
{

}

void GridBuilderImp::meshGeometry(std::string input, int level)
{
    checkLevel(level);

    for (int i = 0; i < rankTasks.size(); i += 2)
    {
        int levelTask = rankTasks[i];

        int index = rankTasks[i + 1];

        if (levelTask == level)
        {
            Geometry* geometry;
            if (StringUtil::endsWith(input, ".stl") || StringUtil::endsWith(input, ".STL"))
            {
                geometry = new Geometry(input, boxes[level][index], transformators[level].get());
                printf("Serialize State\n");
                GeometryMemento memento = geometry->getState();
                serialize(memento, input + "state");
            }
            else {
                *logging::out << logging::Logger::INTERMEDIATE << "start reading memento ...\n";

                clock_t begin = clock();
                geometry = new Geometry();
                GeometryMemento memento;
                deserialize(memento, input);
                geometry->setState(memento);
                clock_t end = clock();

                real time = real(end - begin) / CLOCKS_PER_SEC;
                *logging::out << logging::Logger::INTERMEDIATE << "time reading memento: " << SSTR(time) << "s\n";
            }
            
            if (geometry->size > 0)
                this->gridKernels[level][index]->meshGrid(*geometry);
            //this->gridKernels[level][index]->copyDataFromGPU();
            delete geometry;
        }
    }
}

void GridBuilderImp::deleteSolidNodes()
{
    this->gridKernels[0][0]->deleteSolidNodes();
    this->gridKernels[0][0]->copyDataFromGPU();
}

void GridBuilderImp::flood(Vertex &startFlood, int level)
{
    checkLevel(level);
    this->gridKernels[level][0]->floodFill(startFlood);
}

void GridBuilderImp::createBoundaryConditions()
{
    this->createBCVectors();
}


unsigned int GridBuilderImp::getNumberOfNodes(unsigned int level) const
{
    return (unsigned int) this->gridKernels[level][0]->grid.reducedSize;
}

std::vector<std::vector<std::vector<real> > > GridBuilderImp::getQsValues() const
{
    return this->Qs;
}

int GridBuilderImp::getBoundaryConditionSize(int rb) const
{
    return (int)Qs[rb].size();
}

std::vector<std::string> GridBuilderImp::getTypeOfBoundaryConditions() const
{
    return this->channelBoundaryConditions;
}

void GridBuilderImp::writeGridToVTK(std::string output, int level)
{
    checkLevel(level);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < rankTasks.size(); i += 2) 
    {
        int levelTask = rankTasks[i];
        int index = rankTasks[i + 1];

        if (levelTask == level) {
            std::stringstream ss;
            ss << output << rank;
            GridVTKWriter::writeGridToVTK(this->gridKernels[level][index]->grid, ss.str(), this->transformators[level], true);
        }
    }
}


void GridBuilderImp::writeSimulationFiles(std::string output, BoundingBox<int> &nodesDelete, bool writeFilesBinary, int level)
{
    //checkLevel(level);
    //UnstructuredGridBuilderImp builder;
    //builder.buildUnstructuredGrid(this->gridKernels[level]->grid, nodesDelete);

    //std::vector<Node> coords = builder.getCoordsVec();
    //std::vector<std::vector<std::vector<real> > > qs = builder.getQsValues();
    //SimulationFileWriter::writeSimulationFiles(output, coords, qs, writeFilesBinary, this->gridKernels[level]->grid, this->transformators[level]);
}

std::shared_ptr<GridWrapper> GridBuilderImp::getGridWrapper(int level, int box)
{
    return this->gridKernels[level][box];
}

void GridBuilderImp::checkLevel(int level)
{
    if (level >= boxes.size()) { std::cout << "wrong level input... return to caller\n"; return; }
}


void GridBuilderImp::addGrid(real length, real width, real high, real delta, std::string distribution, std::shared_ptr<Transformator> trans)
{
    this->transformators.push_back(trans);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    this->setCudaDevice(rank);
    this->setNumberOfNodes(length, width, high, delta);
    this->rebuildBoxes();

    if (rank == 0)
        this->sendTasks();
    else
        this->receiveTasks();

    this->createGridKernels(distribution);
}

void GridBuilderImp::createGridKernels(std::string distribution)
{
    for (int i = 0; i < rankTasks.size(); i += 2)
    {
        int level = rankTasks[i];
        int index = rankTasks[i + 1];

        switch (this->device)
        {
        case GenerationDevice::CPU:
            this->gridKernels[level][index] = std::shared_ptr<GridWrapperCPU>(new GridWrapperCPU(this->boxes[level][index], distribution));
            break;
        case GenerationDevice::GPU:
            this->gridKernels[level][index] = std::shared_ptr<GridWrapperGPU>(new GridWrapperGPU(this->boxes[level][index], distribution));
            break;
        }
    }
}


void GridBuilderImp::setNumberOfNodes(real length, real width, real high, real delta)
{
    int nx = (int)(length / delta);
    int ny = (int)(width / delta);
    int nz = (int)(high / delta);

    std::vector<int> dim = { nx,ny,nz };
    this->gridDimensions.push_back(dim);

    this->printMasterInformation(nx, ny, nz);
}

void GridBuilderImp::printMasterInformation(int nx, int ny, int nz)
{
    *logging::out << "global field dimension : " << SSTR(nx) << ", " << SSTR(ny) << ", " << SSTR(nz) << "\n";
    *logging::out << "global field size : " << SSTR(nx * ny * nz) << "\n";
    *logging::out << "------------------------------------------- \n";
}

void GridBuilderImp::setCudaDevice(int rank)
{
    int countDevices;
    cudaGetDeviceCount(&countDevices);
#ifdef __unix__
    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);

    std::string hostnameStr(hostname);
    int device = rank % 4;
    cudaSetDevice(device);
    logging::Logger::getInstance()->logTerminal("Hostname: " + hostnameStr + ", device: " + SSTR(device) + "/" + SSTR(countDevices) + "\n");
#else
    cudaSetDevice(1);
#endif
}

void GridBuilderImp::rebuildBoxes()
{
    int numProcess;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
    int actualLevelSize = (int)this->boxes.size() + 1;
    this->boxes.resize(actualLevelSize);
    this->gridKernels.resize(actualLevelSize);

    //for each level rebuild boxes 
    for (int i = 0; i < actualLevelSize; i++) {
        this->boxes[i] = Partition::getProcessBoxes(numProcess / actualLevelSize, gridDimensions[i][0], gridDimensions[i][1], gridDimensions[i][2]);
        this->gridKernels[i].resize(this->boxes[i].size());
    }
}

void GridBuilderImp::sendTasks()
{
    int numProcess;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
    std::vector<std::vector<int> >tasks = Partition::partitionBoxes(this->boxes, numProcess, transformators);

    rankTasks.resize(tasks[0].size());
    for (int i = 0; i < rankTasks.size(); i++) {
        rankTasks[i] = tasks[0][i];
    }

    for (int i = 1; i < numProcess; i++) {
        int numOfBoxes = (int)tasks[i].size();
        MPI_Send(&numOfBoxes, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&tasks[i][0], numOfBoxes, MPI_INT, i, 1, MPI_COMM_WORLD);
    }
}

void GridBuilderImp::receiveTasks()
{
    int numOfBoxes;
    MPI_Status status;
    MPI_Recv(&numOfBoxes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    rankTasks.resize(numOfBoxes);

    int *boxNumbers = new int[numOfBoxes];
    MPI_Recv(boxNumbers, numOfBoxes, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

    for (int i = 0; i < numOfBoxes; i++) {
        rankTasks[i] = boxNumbers[i];
    }
    delete[] boxNumbers;
}

void GridBuilderImp::writeBoxes(std::string name)
{
    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //UnstructuredGridWrapper writer;
    //for (int i = 0; i < rankTasks.size(); i += 2) {
    //    int level = rankTasks[i];
    //    int index = rankTasks[i + 1];
    //    BoundingBox<real> box = boxes[level][index];
    //    transformators[level]->transformGridToWorld(box);
    //    writer.addBoundingBox(box, rank);
    //}
    //writer.writeUnstructuredGridToFile(name);
}


void GridBuilderImp::getDimensions(int &nx, int &ny, int &nz, const int level) const
{
    nx = this->gridKernels[level][0]->grid.nx;
    ny = this->gridKernels[level][0]->grid.ny;
    nz = this->gridKernels[level][0]->grid.nz;
}

void GridBuilderImp::getNodeValues(real *xCoords, real *yCoords, real *zCoords, unsigned int *neighborX, unsigned int *neighborY, unsigned int *neighborZ, unsigned int *geo, const int level) const
{
    xCoords[0] = 0;
    yCoords[0] = 0;
    zCoords[0] = 0;
    neighborX[0] = 0;
    neighborY[0] = 0;
    neighborZ[0] = 0;
    geo[0] = GEOSOLID;

    Grid grid = this->gridKernels[level][0]->grid;

    for (int i = 0; i < grid.reducedSize; i++)
    {
        unsigned int x, y, z;
        grid.transIndexToCoords(grid.matrixIndex[i], x, y, z);

        xCoords[i + 1] = (real)x;
        yCoords[i + 1] = (real)y;
        zCoords[i + 1] = (real)z;
        neighborX[i + 1] = (unsigned int)(grid.neighborIndexX[grid.matrixIndex[i]] + 1);
        neighborY[i + 1] = (unsigned int)(grid.neighborIndexY[grid.matrixIndex[i]] + 1);
        neighborZ[i + 1] = (unsigned int)(grid.neighborIndexZ[grid.matrixIndex[i]] + 1);
        geo[i + 1] = (unsigned int)grid.isSolid(grid.matrixIndex[i]) ? GEOSOLID : GEOFLUID;
    }
}

void GridBuilderImp::setQs(real** q27, int* k, int channelSide, unsigned int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        k[index] = (int)Qs[channelSide][index][0];
        for (int column = 1; column < Qs[channelSide][index].size(); column++) {
            q27[column - 1][index] = Qs[channelSide][index][column];

        }
    }
}

void GridBuilderImp::setOutflowValues(real* RhoBC, int* kN, int channelSide, int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        RhoBC[index] = 0.0;
        kN[index] = 0;
    }
}

void GridBuilderImp::setVelocityValues(real* vx, real* vy, real* vz, int channelSide, int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        vx[index] = 0.0;
        vy[index] = 0.0;
        vz[index] = 0.0;
    }
}

void GridBuilderImp::setPressValues(real* RhoBC, int* kN, int channelSide, int level) const
{
    for (int index = 0; index < Qs[channelSide].size(); index++) {
        RhoBC[index] = 0.0;
        kN[index] = 0;
    }
}


/*#################################################################################*/
/*---------------------------------private methods---------------------------------*/
/*---------------------------------------------------------------------------------*/
void GridBuilderImp::createBCVectors()
{
    Grid grid = this->gridKernels[0][0]->grid;
    for (int i = 0; i < grid.reducedSize; i++)
    {
        unsigned int x, y, z;
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

void GridBuilderImp::addShortQsToVector(int index)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    Grid grid = this->gridKernels[0][0]->grid;
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

void GridBuilderImp::addQsToVector(int index)
{
    std::vector<real> qNode;
    qNode.push_back((real)index);

    Grid grid = this->gridKernels[0][0]->grid;
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

void GridBuilderImp::fillRBForNode(int x, int y, int z, int index, int direction, int directionSign, int rb)
{
    uint32_t qKey = 0;
    std::vector<real> qNode;

    Grid grid = this->gridKernels[0][0]->grid;
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

void GridBuilderImp::writeArrows(std::string fileName, std::shared_ptr<ArrowTransformator> trans) const
{
    Grid grid = this->gridKernels[0][0]->grid;
    //std::shared_ptr<PolyDataWriterWrapper> writer = std::shared_ptr<PolyDataWriterWrapper>(new PolyDataWriterWrapper());
    for (int index = 0; index < Qs[GEOMQS].size(); index++)
    {
        Vertex startNode = getVertex(getMatrixIndex(index));
        //for (int qi = grid.d.dir_start; qi <= grid.d.dir_end; qi++)
            //writeArrow(index, qi, startNode, trans, writer);
    }
    //writer->writePolyDataToFile(fileName);
}

void GridBuilderImp::writeArrow(const int i, const int qi, const Vertex& startNode, std::shared_ptr<const ArrowTransformator> trans/*, std::shared_ptr<PolyDataWriterWrapper> writer*/) const
{
    Grid grid = this->gridKernels[0][0]->grid;
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

Vertex GridBuilderImp::getVertex(int matrixIndex) const
{
    unsigned int x, y, z;
    this->gridKernels[0][0]->grid.transIndexToCoords(matrixIndex, x, y, z);
    return Vertex((real)x, (real)y, (real)z);
}

int GridBuilderImp::getMatrixIndex(int i) const
{
    int index = (int)Qs[GEOMQS][i][0];
    return this->gridKernels[0][0]->grid.matrixIndex[index];
}

