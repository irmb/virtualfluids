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

#include <basics/utilities/UbTuple.h>


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

std::shared_ptr<LevelGridBuilder> LevelGridBuilder::makeShared(Device device, const std::string& d3qxx)
{
    return SPtr<LevelGridBuilder>(new LevelGridBuilder(device, d3qxx));
}

void LevelGridBuilder::setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz)
{
    if (sideType == SideType::GEOMETRY)
        setVelocityGeometryBoundaryCondition(vx, vy, vz);
    else
    {

        SPtr<VelocityBoundaryCondition> velocityBoundaryCondition = VelocityBoundaryCondition::make(vx, vy, vz);

        auto side = SideFactory::make(sideType);

        velocityBoundaryConditions.push_back(velocityBoundaryCondition);
        velocityBoundaryCondition->side = side;
    }
}

void LevelGridBuilder::setVelocityGeometryBoundaryCondition(real vx, real vy, real vz)
{
    geometryBoundaryCondition->hasValues = true;
    geometryBoundaryCondition->vx = vx;
    geometryBoundaryCondition->vy = vy;
    geometryBoundaryCondition->vz = vz;
}

void LevelGridBuilder::setPressureBoundaryCondition(SideType sideType, real rho)
{
    SPtr<PressureBoundaryCondition> pressureBoundaryCondition = PressureBoundaryCondition::make(rho);

    auto side = SideFactory::make(sideType);

    pressureBoundaryConditions.push_back(pressureBoundaryCondition);
    pressureBoundaryCondition->side = side;
}

void LevelGridBuilder::setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z)
{
    grids[0]->setPeriodicity(periodic_X, periodic_Y, periodic_Z);
}

void LevelGridBuilder::setNoSlipBoundaryCondition(SideType sideType)
{
    SPtr<VelocityBoundaryCondition> noSlipBoundaryCondition = VelocityBoundaryCondition::make(0.0, 0.0, 0.0);

    auto side = SideFactory::make(sideType);

    noSlipBoundaryConditions.push_back(noSlipBoundaryCondition);
    noSlipBoundaryCondition->side = side;
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
            indices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->indices[i]) +1;  

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
            indices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->indices[i]) + 1;

            neighborIndices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->neighborIndices[i]) + 1;

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

void LevelGridBuilder::getPressureQs(real* qs[27], int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : this->pressureBoundaryConditions)
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


uint LevelGridBuilder::getGeometrySize(int level) const
{
    if (geometryBoundaryCondition)
        return geometryBoundaryCondition->indices.size();
    
    return 0;
}

void LevelGridBuilder::getGeometryIndices(int* indices, int level) const
{
    for (uint i = 0; i < geometryBoundaryCondition->indices.size(); i++)
    {
        indices[i] = grids[level]->getSparseIndex(geometryBoundaryCondition->indices[i]) + 1;
    }
}

bool LevelGridBuilder::hasGeometryValues() const
{
    return geometryBoundaryCondition->hasValues;
}


void LevelGridBuilder::getGeometryValues(real* vx, real* vy, real* vz, int level) const
{
    for (uint i = 0; i < geometryBoundaryCondition->indices.size(); i++)
    {
        vx[i] = geometryBoundaryCondition->vx;
        vy[i] = geometryBoundaryCondition->vy;
        vz[i] = geometryBoundaryCondition->vz;
    }
}


void LevelGridBuilder::getGeometryQs(real* qs[27], int level) const
{
    for (int i = 0; i < geometryBoundaryCondition->indices.size(); i++)
    {
        for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
        {
            qs[dir][i] = geometryBoundaryCondition->qs[i][dir];
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






#include <fstream>
using namespace std;
void writeLines(std::string filename, std::vector<UbTupleFloat3> nodes, std::vector<UbTupleInt2> lines)
{
    string vtkfilename = filename + ".bin.vtu";

    ofstream out(vtkfilename.c_str(), ios::out | ios::binary);

    int nofNodes = (int)nodes.size();
    int nofCells = (int)lines.size();

    int bytesPerByteVal = 4; //==sizeof(int)
    int bytesPoints = 3 /*x1/x2/x3        */ * nofNodes * sizeof(float);
    int bytesCellConnectivty = 2 /*nodes per line */ * nofCells * sizeof(int);
    int bytesCellOffsets = 1 /*offset per line */ * nofCells * sizeof(int);
    int bytesCellTypes = 1 /*type of line */ * nofCells * sizeof(unsigned char);

    int offset = 0;
    //VTK FILE
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >" << "\n";
    out << "   <UnstructuredGrid>" << "\n";
    out << "      <Piece NumberOfPoints=\"" << nofNodes << "\" NumberOfCells=\"" << nofCells << "\">\n";

    //POINTS SECTION
    out << "         <Points>\n";
    out << "            <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"  />\n";
    out << "         </Points>\n";
    offset += (bytesPerByteVal + bytesPoints);

    //CELLS SECTION
    out << "         <Cells>\n";
    out << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += (bytesPerByteVal + bytesCellConnectivty);
    out << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << offset << "\" />\n";
    offset += (bytesPerByteVal + bytesCellOffsets);
    out << "            <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\" />\n ";
    offset += (bytesPerByteVal + bytesCellTypes);
    out << "         </Cells>\n";

    out << "      </Piece>\n";
    out << "   </UnstructuredGrid>\n";

    // AppendedData SECTION
    out << "   <AppendedData encoding=\"raw\">\n";
    out << "_";

    //POINTS SECTION
    out.write((char*)&bytesPoints, bytesPerByteVal);
    for (int n = 0; n < nofNodes; n++)
    {
        out.write((char*)&val<1>(nodes[n]), sizeof(float));
        out.write((char*)&val<2>(nodes[n]), sizeof(float));
        out.write((char*)&val<3>(nodes[n]), sizeof(float));
    }

    //CELLS SECTION
    //cellConnectivity
    out.write((char*)&bytesCellConnectivty, bytesPerByteVal);
    for (int c = 0; c < nofCells; c++)
    {
        out.write((char*)&val<1>(lines[c]), sizeof(int));
        out.write((char*)&val<2>(lines[c]), sizeof(int));

    }

    //cellOffsets
    out.write((char*)&bytesCellOffsets, bytesPerByteVal);
    int itmp;
    for (int c = 1; c <= nofCells; c++)
    {
        itmp = 2 * c;
        out.write((char*)&itmp, sizeof(int));
    }

    //cellTypes
    out.write((char*)&bytesCellTypes, bytesPerByteVal);
    unsigned char vtkCellType = 3;
    for (int c = 0; c < nofCells; c++)
    {
        out.write((char*)&vtkCellType, sizeof(unsigned char));
    }
    out << "\n</AppendedData>\n";
    out << "</VTKFile>";
    out << endl;
    out.close();
}

void LevelGridBuilder::writeArrows(std::string fileName) const
{
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleInt2> cells;

    int actualNodeNumber = 0;
    for (int index = 0; index < geometryBoundaryCondition->indices.size(); index++)
    {
        Vertex startNode = getVertex(geometryBoundaryCondition->indices[index]);
        for (int qi = 0; qi <= 26; qi++)
        {
            real qval = geometryBoundaryCondition->qs[index][qi];
            if (qval > 0.0f)
            {
                Vertex dir((real)grids[0]->getDirection()[qi * DIMENSION + 0], (real)grids[0]->getDirection()[qi* DIMENSION + 1], (real)grids[0]->getDirection()[qi * DIMENSION + 2]);
                Vertex nodeOnGeometry(startNode + (dir * qval));

                nodes.push_back(makeUbTuple(float(startNode.x), float(startNode.y), float(startNode.z)));
                nodes.push_back(makeUbTuple(float(nodeOnGeometry.x), float(nodeOnGeometry.y), float(nodeOnGeometry.z)));
                actualNodeNumber += 2;
                cells.push_back(makeUbTuple(actualNodeNumber - 2, actualNodeNumber - 1));
            }
        }
    }

    writeLines(fileName, nodes, cells);
}


Vertex LevelGridBuilder::getVertex(int matrixIndex) const
{
    real x, y, z;
    this->grids[0]->transIndexToCoords(matrixIndex, x, y, z);
    return Vertex(x,y,z);
}

