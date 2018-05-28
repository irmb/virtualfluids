#include "LevelGridBuilder.h"

#include <stdio.h>
#include <iostream>
#include <GridGenerator/grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h>
#include <GridGenerator/grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h>
#include <GridGenerator/grid/partition/Partition.h>

#include <GridGenerator/geometries/Triangle/Triangle.h>
#include <GridGenerator/geometries/BoundingBox/BoundingBox.h>
#include <GridGenerator/geometries/TriangularMesh/TriangularMesh.h>


#include <GridGenerator/io/GridVTKWriter/GridVTKWriter.h>
#include <GridGenerator/io/SimulationFileWriter/SimulationFileWriter.h>
#include <GridGenerator/io/VTKWriterWrapper/UnstructuredGridWrapper.h>
#include <GridGenerator/io/VTKWriterWrapper/PolyDataWriterWrapper.h>

#include <GridGenerator/grid/NodeValues.h>

#include <GridGenerator/geometries/Arrow/ArrowImp.h>
#include <GridGenerator/utilities/transformator/ArrowTransformator.h>

#include <utilities/logger/Logger.h>

#include <GridGenerator/grid/GridFactory.h>
#include "grid/GridInterface.h"
#include "grid/Grid.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/BoundaryConditions/Side.h"


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


    sideIsSet = { { SideType::MX, false },{ SideType::PX, false },{ SideType::MY, false },{ SideType::PY, false },{ SideType::MZ, false },{ SideType::PZ, false } };

}

std::shared_ptr<LevelGridBuilder> LevelGridBuilder::makeShared(Device device, const std::string& d3qxx)
{
    return SPtr<LevelGridBuilder>(new LevelGridBuilder(device, d3qxx));
}

void LevelGridBuilder::setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz)
{
    sideIsSet[sideType] = true;

    SPtr<VelocityBoundaryCondition> velocityBoundaryCondition = VelocityBoundaryCondition::make(vx, vy, vz);

    auto side = SideFactory::make(sideType);

    side->setPeriodicy(grids[0]);
    velocityBoundaryConditions.push_back(velocityBoundaryCondition);
    velocityBoundaryCondition->side = side;
}

void LevelGridBuilder::setPressureBoundaryCondition(SideType sideType, real rho)
{
    sideIsSet[sideType] = true;

    SPtr<PressureBoundaryCondition> pressureBoundaryCondition = PressureBoundaryCondition::make(rho);

    auto side = SideFactory::make(sideType);

    side->setPeriodicy(grids[0]);
    pressureBoundaryConditions.push_back(pressureBoundaryCondition);
    pressureBoundaryCondition->side = side;
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

SPtr<Grid> LevelGridBuilder::getGrid(uint level)
{
    return grids[level];
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


std::shared_ptr<Grid> LevelGridBuilder::getGrid(int level, int box)
{
    return this->grids[level];
}

void LevelGridBuilder::checkLevel(int level)
{
    if (level >= grids.size())
    { 
        std::cout << "wrong level input... return to caller\n";
        return; 
    }
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


uint LevelGridBuilder::getVelocitySize(int level) const
{
    uint size = 0;
    for (auto boundaryCondition : this->velocityBoundaryConditions)
    {
        size += uint(boundaryCondition->indices.size());
    }
    return size;
}

void LevelGridBuilder::getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : this->velocityBoundaryConditions)
    {
        for(int i = 0; i < boundaryCondition->indices.size(); i++)
        {
            indices[allIndicesCounter] = boundaryCondition->indices[i];  

            vx[allIndicesCounter] = boundaryCondition->vx;
            vy[allIndicesCounter] = boundaryCondition->vy;
            vz[allIndicesCounter] = boundaryCondition->vz;
            allIndicesCounter++;
        }
    }
}

uint LevelGridBuilder::getPressureSize(int level) const
{
    uint size = 0;
    for (auto boundaryCondition : this->pressureBoundaryConditions)
    {
        size += uint(boundaryCondition->indices.size());
    }
    return size;
}

void LevelGridBuilder::getPressureValues(real* rho, int* indices, int* neighborIndices, int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : this->pressureBoundaryConditions)
    {
        for (int i = 0; i < boundaryCondition->indices.size(); i++)
        {
            indices[allIndicesCounter] = boundaryCondition->indices[i];

            neighborIndices[allIndicesCounter] = boundaryCondition->neighborIndices[i];;

            rho[allIndicesCounter] = boundaryCondition->rho;
            allIndicesCounter++;
        }
    }
}


void LevelGridBuilder::getVelocityQs(real* qs[27], int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : this->velocityBoundaryConditions)
    {
        for (int i = 0; i < boundaryCondition->indices.size(); i++)
        {

            for (int dir = 0; dir < grids[level]->getEndDirection(); dir++)
            {
                if (grids[level]->getDirection()[dir * DIMENSION + boundaryCondition->side->getCoordinate()] != boundaryCondition->side->getDirection())
                    qs[dir][allIndicesCounter] = -1.0;
                else
                    qs[dir][allIndicesCounter] = 0.5;
            }
            allIndicesCounter++;
        }
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

