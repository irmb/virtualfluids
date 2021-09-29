#include "LevelGridBuilder.h"

#include <stdio.h>
#include <iostream>
#include <algorithm>

#include "geometries/Arrow/ArrowImp.h"
#include "geometries/Triangle/Triangle.h"
#include "geometries/BoundingBox/BoundingBox.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/GridStrategy/GridGpuStrategy/GridGpuStrategy.h"
#include "grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "grid/partition/Partition.h"
#include "grid/NodeValues.h"
#include "grid/GridFactory.h"
#include "grid/GridInterface.h"
#include "grid/Grid.h"

#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "io/VTKWriterWrapper/UnstructuredGridWrapper.h"
#include "io/VTKWriterWrapper/PolyDataWriterWrapper.h"
#include "io/QLineWriter.h"

#include "utilities/transformator/ArrowTransformator.h"
#include "utilities/communication.h"

#define GEOFLUID 19
#define GEOSOLID 16

using namespace vf::gpu;

LevelGridBuilder::LevelGridBuilder(Device device, const std::string& d3qxx) : device(device), d3qxx(d3qxx)
{
    this->communicationProcesses[CommunicationDirections::MX] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::PX] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::MY] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::PY] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::MZ] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::PZ] = INVALID_INDEX;
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
        for (uint level = 0; level < getNumberOfGridLevels(); level++)
        {
            SPtr<VelocityBoundaryCondition> velocityBoundaryCondition = VelocityBoundaryCondition::make(vx, vy, vz);

            auto side = SideFactory::make(sideType);

            velocityBoundaryCondition->side = side;
            velocityBoundaryCondition->side->addIndices(grids, level, velocityBoundaryCondition);

            velocityBoundaryCondition->fillVelocityLists();

            boundaryConditions[level]->velocityBoundaryConditions.push_back(velocityBoundaryCondition);

            *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Velocity BC on level " << level << " with " << (int)velocityBoundaryCondition->indices.size() <<"\n";
        }
    }
}

void LevelGridBuilder::setVelocityGeometryBoundaryCondition(real vx, real vy, real vz)
{
    geometryHasValues = true;

    for (uint level = 0; level < getNumberOfGridLevels(); level++)
    {
		if (boundaryConditions[level]->geometryBoundaryCondition != nullptr)
		{
			boundaryConditions[level]->geometryBoundaryCondition->vx = vx;
			boundaryConditions[level]->geometryBoundaryCondition->vy = vy;
			boundaryConditions[level]->geometryBoundaryCondition->vz = vz;
			boundaryConditions[level]->geometryBoundaryCondition->side->addIndices(grids, level, boundaryConditions[level]->geometryBoundaryCondition);

            boundaryConditions[level]->geometryBoundaryCondition->fillVelocityLists();

            *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Geometry Velocity BC on level " << level << " with " << (int)boundaryConditions[level]->geometryBoundaryCondition->indices.size() <<"\n";
		}
    }
}

void LevelGridBuilder::setPressureBoundaryCondition(SideType sideType, real rho)
{
    for (uint level = 0; level < getNumberOfGridLevels(); level++)
    {
        SPtr<PressureBoundaryCondition> pressureBoundaryCondition = PressureBoundaryCondition::make(rho);

        auto side = SideFactory::make(sideType);
        pressureBoundaryCondition->side = side;
        pressureBoundaryCondition->side->addIndices(grids, level, pressureBoundaryCondition);

        boundaryConditions[level]->pressureBoundaryConditions.push_back(pressureBoundaryCondition);

        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Pressure BC on level " << level << " with " << (int)pressureBoundaryCondition->indices.size() <<"\n";
    }
}

void LevelGridBuilder::setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z)
{
    for( uint level = 0; level < this->grids.size(); level++ )
        grids[level]->setPeriodicity(periodic_X, periodic_Y, periodic_Z);
}

void LevelGridBuilder::setNoSlipBoundaryCondition(SideType sideType)
{
    for (uint level = 0; level < getNumberOfGridLevels(); level++)
    {
        SPtr<VelocityBoundaryCondition> noSlipBoundaryCondition = VelocityBoundaryCondition::make(0.0, 0.0, 0.0);

        auto side = SideFactory::make(sideType);

        noSlipBoundaryCondition->side = side;
        noSlipBoundaryCondition->side->addIndices(grids, level, noSlipBoundaryCondition);

        boundaryConditions[level]->noSlipBoundaryConditions.push_back(noSlipBoundaryCondition);
    }
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall)
{
    for( uint level = 0; level < this->grids.size(); level++ )
        grids[level]->setEnableFixRefinementIntoTheWall( enableFixRefinementIntoTheWall );
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::setCommunicationProcess(int direction, uint process)
{
    this->communicationProcesses[direction] = process;
}

GRIDGENERATOR_EXPORT uint LevelGridBuilder::getCommunicationProcess(int direction)
{
    return this->communicationProcesses[direction];
}



void LevelGridBuilder::copyDataFromGpu()
{
    for (const auto& grid : grids)
    {
        auto gridGpuStrategy = std::dynamic_pointer_cast<GridGpuStrategy>(grid->getGridStrategy());
        if(gridGpuStrategy)
            gridGpuStrategy->copyDataFromGPU(std::static_pointer_cast<GridImp>(grid));
    }
        
}

LevelGridBuilder::~LevelGridBuilder()
{
    for (const auto& grid : grids)
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
    for (const auto& grid : grids)
    {
        gridX.push_back(int(grid->getNumberOfNodesX()));
        gridY.push_back(int(grid->getNumberOfNodesY()));
        gridZ.push_back(int(grid->getNumberOfNodesZ()));

        distX.push_back(int(grid->getStartX()));
        distY.push_back(int(grid->getStartY()));
        distZ.push_back(int(grid->getStartZ()));
    }
}


uint LevelGridBuilder::getNumberOfGridLevels() const
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

void LevelGridBuilder::getGridInterfaceIndicesFCBorderBulk(uint *iCellFccBorder, uint *&iCellFccBulk, uint *iCellFcfBorder, uint *&iCellFcfBulk,
                                                            uint &intFCBorderKfc, uint &intFCBulkKfc, int level) const
{
    this->grids[level]->getGridInterfaceIndicesFCBorderBulk(iCellFccBorder, iCellFccBulk,
                                                             iCellFcfBorder, iCellFcfBulk,
                                                             intFCBorderKfc, intFCBulkKfc, level);
}

void LevelGridBuilder::getOffsetFC(real * xOffFC, real * yOffFC, real * zOffFC, int level)
{
    for (uint i = 0; i < getNumberOfNodesFC(level); i++)
    {
        uint offset = this->grids[level]->getFC_offset()[i];

        xOffFC[i] = - this->grids[level]->getDirection()[ 3*offset + 0 ];
        yOffFC[i] = - this->grids[level]->getDirection()[ 3*offset + 1 ];
        zOffFC[i] = - this->grids[level]->getDirection()[ 3*offset + 2 ];
    }
}

void LevelGridBuilder::getOffsetCF(real * xOffCF, real * yOffCF, real * zOffCF, int level)
{
    for (uint i = 0; i < getNumberOfNodesCF(level); i++)
    {
        uint offset = this->grids[level]->getCF_offset()[i];

        xOffCF[i] = - this->grids[level]->getDirection()[ 3*offset + 0 ];
        yOffCF[i] = - this->grids[level]->getDirection()[ 3*offset + 1 ];
        zOffCF[i] = - this->grids[level]->getDirection()[ 3*offset + 2 ];
    }
}

GRIDGENERATOR_EXPORT uint LevelGridBuilder::getNumberOfSendIndices(int direction, uint level)
{
    return this->grids[level]->getNumberOfSendNodes(direction);
}

GRIDGENERATOR_EXPORT uint LevelGridBuilder::getNumberOfReceiveIndices(int direction, uint level)
{
    return this->grids[level]->getNumberOfReceiveNodes(direction);
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::getSendIndices(int * sendIndices, int direction, int level)
{
    SPtr<Grid> grid = this->grids[level];
    for( uint i = 0; i < getNumberOfSendIndices(direction, level); i++ )
    {
        sendIndices[i] = grid->getSparseIndex( grid->getSendIndex(direction, i) ) + 1;
    }
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::getReceiveIndices(int * receiveIndices, int direction, int level)
{
    SPtr<Grid> grid = this->grids[level];
    for( uint i = 0; i < getNumberOfReceiveIndices(direction, level); i++ )
    {
        receiveIndices[i] = grid->getSparseIndex( grid->getReceiveIndex(direction, i) ) + 1;
    }
}

GRIDGENERATOR_EXPORT std::vector<uint>
LevelGridBuilder::getAndReorderSendIndices(int *sendIndices, uint &numberOfSendNeighborsAfterFtoC, uint *iCellFCC,
                                           uint sizeOfICellFCCBorder, uint *iCellCFC, uint sizeOfICellCFC,
                                           uint *neighborX, uint *neighborY, uint *neighborZ, int direction, int level)
{
    std::vector<uint> sendIndicesForCommAfterFtoCPositions;
    getSendIndices(sendIndices, direction, level);
    reorderSendIndicesForCommAfterFtoC(sendIndices, numberOfSendNeighborsAfterFtoC, iCellFCC, sizeOfICellCFC, iCellCFC,
                                       sizeOfICellCFC, neighborX, neighborY, neighborZ, direction, level,
                                       sendIndicesForCommAfterFtoCPositions);
    return sendIndicesForCommAfterFtoCPositions;
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::getAndReorderReceiveIndices(int *recvIndices,
                                                                        uint &numberOfRecvNeighborsAfterFtoC,
                                                                        uint *iCellFCCBorder, uint sizeOfICellFCCBorder,
                                                                        int direction, int level,
                                                                        bool receiveIndicesNeedToBeReordered)
{
    getReceiveIndices(recvIndices, direction, level);
    if (receiveIndicesNeedToBeReordered)
        reorderRecvIndexForCommAfterFtoC(recvIndices, numberOfRecvNeighborsAfterFtoC, iCellFCCBorder,
                                         sizeOfICellFCCBorder, direction, level);
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::reorderSendIndicesForCommAfterFtoC(
    int *sendIndices, uint &numberOfSendNeighborsAfterFtoC, uint *iCellFCC, uint sizeOfICellFCC, uint *iCellCFC,
    uint sizeOfICellCFC, uint *neighborX, uint *neighborY, uint *neighborZ, int direction, int level,
    std::vector<uint> &sendIndicesForCommAfterFtoCPositions)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "reorder send indices for communication after fine to coarse: level: " << level
                  << " direction: " << direction;
    if (sizeOfICellFCC == 0 || sizeOfICellCFC == 0)
        *logging::out << logging::Logger::LOGGER_ERROR
                      << "reorderSendIndicesForCommAfterFtoC(): iCellFCC needs to be inititalized before calling "
                         "this function "
                      << "\n";
    uint numberOfSendIndices = getNumberOfSendIndices(direction, level);
    if (numberOfSendIndices == 0) {
        numberOfSendNeighborsAfterFtoC = 0;
        return;
    }

    int sparseIndexSend;
    bool isInICells;
    std::vector<int> sendIndicesAfterFtoC;
    std::vector<int> sendIndicesOther;
    int neighborToAddX, neighborToAddY, neighborToAddZ;

    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        neighborToAddX = neighborToAddY = neighborToAddZ = -1;
        sparseIndexSend = sendIndices[posInSendIndices];

        // check if sparse index is in ICellFCC
        isInICells = false;
        for (uint j = 0; j < sizeOfICellFCC; j++) {
            if (sparseIndexSend < 0)
                break;
            if (iCellFCC[j] == (uint) sparseIndexSend) {
                isInICells = true;
                break;
            }
        }
        // check if sparse index is in ICellCFC 
        for (uint j = 0; j < sizeOfICellCFC; j++) {
            if (sparseIndexSend < 0)
                break;
            if (iCellCFC[j] == (uint)sparseIndexSend) {
                isInICells = true;                   
                // also find neighbors
                if (direction != 0 && direction != 1)
                    neighborToAddX = neighborX[sparseIndexSend];
                if (direction != 2 && direction != 3)
                    neighborToAddY = neighborY[sparseIndexSend];
                if (direction != 4 && direction != 5)
                    neighborToAddZ = neighborZ[sparseIndexSend];
                break;
            }
        }

        // add index to corresponding vectors but omit indices which are already in sendIndicesAfterFtoC
        if (isInICells) {
            if (std::find(sendIndicesAfterFtoC.begin(), sendIndicesAfterFtoC.end(), sparseIndexSend) ==
                sendIndicesAfterFtoC.end()) {
                sendIndicesAfterFtoC.push_back(sparseIndexSend);
                sendIndicesForCommAfterFtoCPositions.push_back(posInSendIndices);
            }
        }   

        // also add neighbors
        if (neighborToAddX != -1)
            findIfSparseIndexIsInSendIndicesAndAddToVectors(neighborToAddX, sendIndices, numberOfSendIndices,
                                                            sendIndicesAfterFtoC, sendIndicesForCommAfterFtoCPositions);
        if (neighborToAddY != -1)
            findIfSparseIndexIsInSendIndicesAndAddToVectors(neighborToAddY, sendIndices, numberOfSendIndices,
                                                            sendIndicesAfterFtoC, sendIndicesForCommAfterFtoCPositions);
        if (neighborToAddZ != -1)
            findIfSparseIndexIsInSendIndicesAndAddToVectors(neighborToAddZ, sendIndices, numberOfSendIndices,
                                                            sendIndicesAfterFtoC, sendIndicesForCommAfterFtoCPositions);
    }

    numberOfSendNeighborsAfterFtoC = (uint) sendIndicesAfterFtoC.size();
    
    // add sparseIndices not in sendIndicesAfterFtoC to sendIndicesOther
    for (uint posInSendIndices = 0; posInSendIndices < numberOfSendIndices; posInSendIndices++) {
        sparseIndexSend = sendIndices[posInSendIndices];
        if (std::find(sendIndicesAfterFtoC.begin(), sendIndicesAfterFtoC.end(), sparseIndexSend) ==
            sendIndicesAfterFtoC.end())
            sendIndicesOther.push_back(sparseIndexSend);
    }


    // copy new vectors back to sendIndices array
    for (uint i = 0; i < numberOfSendNeighborsAfterFtoC; i++)
        sendIndices[i] = sendIndicesAfterFtoC[i];
    for (uint i = 0; i < sendIndicesOther.size(); i++)
        sendIndices[i + numberOfSendNeighborsAfterFtoC] = sendIndicesOther[i];

    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "... numberOfSendNeighborsAfterFtoC: " << numberOfSendNeighborsAfterFtoC << "\n";

    if (numberOfSendNeighborsAfterFtoC + sendIndicesOther.size() != numberOfSendIndices) {
        *logging::out << logging::Logger::LOGGER_ERROR
                      << "reorderSendIndicesForCommAfterFtoC(): incorrect number of nodes"
                      << "\n";
        std::cout << "numberOfSendNeighborsAfterFtoC = " << numberOfSendNeighborsAfterFtoC
                  << ", sendIndicesOther.size() = " << sendIndicesOther.size()
                  << ", numberOfSendIndices = " << numberOfSendIndices << std::endl;
    }
}

void LevelGridBuilder::findIfSparseIndexIsInSendIndicesAndAddToVectors(
    int sparseIndex, int *sendIndices, uint numberOfSendIndices,
    std::vector<int>& sendIndicesAfterFtoC,
    std::vector<uint>& sendIndicesForCommAfterFtoCPositions) const
{
    int sparseIndexSendForNeighborSearch;
    for (uint j = 0; j < numberOfSendIndices; j++) {
        sparseIndexSendForNeighborSearch = sendIndices[j];
        if (sparseIndex == sparseIndexSendForNeighborSearch) {
            if (std::find(sendIndicesAfterFtoC.begin(), sendIndicesAfterFtoC.end(), sparseIndexSendForNeighborSearch) ==
                sendIndicesAfterFtoC.end()) {
                sendIndicesAfterFtoC.push_back(sparseIndexSendForNeighborSearch);
                sendIndicesForCommAfterFtoCPositions.push_back(j);
            } 
            break;
        }
    }
}


GRIDGENERATOR_EXPORT void LevelGridBuilder::reorderRecvIndexForCommAfterFtoC(int *recvIndices, uint &numberOfRecvNeighborsAfterFtoC,
                                                                             uint *iCellFCCBorder, uint sizeOfICellFCCBorder, 
                                                                             int direction, int level)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "reorder receive indices for communication after fine to coarse: level: " << level
                  << " direction: " << direction;
    if (sizeOfICellFCCBorder == 0)
        *logging::out << logging::Logger::LOGGER_ERROR
                      << "reorderRecvIndexForCommAfterFtoC(): iCellFCC needs to be inititalized before calling "
                         "this function "
                      << "\n";

    uint numberOfSendIndices = getNumberOfReceiveIndices(direction, level);
    if (numberOfSendIndices == 0) {
        numberOfRecvNeighborsAfterFtoC = 0;
        return;
    }

    int sparseIndexRecv;
    bool isInICellFCCBorder;
    std::vector<int> recvIndicesAfterFtoC;
    std::vector<int> recvIndicesOther;

    for (uint i = 0; i < numberOfSendIndices; i++) {
        sparseIndexRecv = recvIndices[i];

        // check if sparse index is in ICellFCC border
        isInICellFCCBorder = false;
        for (uint j = 0; j < sizeOfICellFCCBorder; j++) {
            if (sparseIndexRecv < 0)
                continue;
            if (iCellFCCBorder[j] == (uint) sparseIndexRecv) {
                isInICellFCCBorder = true;
                break;
            }
        }

        // add index to corresponding vector
        if (isInICellFCCBorder)
            recvIndicesAfterFtoC.push_back(sparseIndexRecv);
        else
            recvIndicesOther.push_back(sparseIndexRecv);
    }

    numberOfRecvNeighborsAfterFtoC = (uint) recvIndicesAfterFtoC.size();

    // copy new vectors back to receiveIndices array
    for (uint i = 0; i < numberOfRecvNeighborsAfterFtoC; i++)
        recvIndices[i] = recvIndicesAfterFtoC[i];
    for (uint i = 0; i < recvIndicesOther.size(); i++)
        recvIndices[i + numberOfRecvNeighborsAfterFtoC] = recvIndicesOther[i];

    *logging::out << logging::Logger::INFO_INTERMEDIATE
                  << "... numberOfRecvNeighborsAfterFtoC: " << numberOfRecvNeighborsAfterFtoC << "\n";

    bool numberOfNodesIsCorrect = numberOfRecvNeighborsAfterFtoC + recvIndicesOther.size() == numberOfSendIndices;
    std::cout << "correct number of nodes?: " << numberOfNodesIsCorrect << std::endl;
}

uint LevelGridBuilder::getNumberOfNodes(unsigned int level) const
{
    return grids[level]->getSparseSize();
}

std::shared_ptr<Grid> LevelGridBuilder::getGrid(int level, int box)
{
    return this->grids[level];
}

void LevelGridBuilder::checkLevel(int level)
{
    if (level >= (int)grids.size())
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

void LevelGridBuilder::getNodeValues(real *xCoords, real *yCoords, real *zCoords, 
                                     uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, 
                                     uint *geo, const int level) const
{
    grids[level]->getNodeValues(xCoords, yCoords, zCoords, neighborX, neighborY, neighborZ, neighborNegative, geo);
}

//TODO: add getSlipSize...

GRIDGENERATOR_EXPORT void LevelGridBuilder::getFluidNodeIndices(uint *fluidNodeIndices, const int level) const 
{ 
    grids[level]->getFluidNodeIndices(fluidNodeIndices);
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::getFluidNodeIndicesBorder(uint *fluidNodeIndices, const int level) const
{
    grids[level]->getFluidNodeIndicesBorder(fluidNodeIndices);
}

uint LevelGridBuilder::getNumberOfFluidNodes(unsigned int level) const 
{
    return grids[level]->getNumberOfFluidNodes(); 
}

GRIDGENERATOR_EXPORT uint LevelGridBuilder::getNumberOfFluidNodesBorder(unsigned int level) const
{
    return grids[level]->getNumberOfFluidNodesBorder();
}

uint LevelGridBuilder::getVelocitySize(int level) const
{
    uint size = 0;
    for (auto boundaryCondition : boundaryConditions[level]->velocityBoundaryConditions)
    {
        size += uint(boundaryCondition->indices.size());
    }
    return size;
}

//TODO: add getSlipIndices...


void LevelGridBuilder::getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->velocityBoundaryConditions)
    {
        for(std::size_t i = 0; i < boundaryCondition->indices.size(); i++)
        {
            indices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->indices[i]) +1;  
            vx[allIndicesCounter] = (real)boundaryCondition->getVx((uint)i);
            vy[allIndicesCounter] = (real)boundaryCondition->getVy((uint)i);
            vz[allIndicesCounter] = (real)boundaryCondition->getVz((uint)i);
            allIndicesCounter++;
        }
    }
}

uint LevelGridBuilder::getPressureSize(int level) const
{
    uint size = 0;
    for (auto boundaryCondition : boundaryConditions[level]->pressureBoundaryConditions)
    {
        size += uint(boundaryCondition->indices.size());
    }
    return size;
}

void LevelGridBuilder::getPressureValues(real* rho, int* indices, int* neighborIndices, int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->pressureBoundaryConditions)
    {
        for (std::size_t i = 0; i < boundaryCondition->indices.size(); i++)
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
    for (auto boundaryCondition : boundaryConditions[level]->velocityBoundaryConditions)
    {
        for ( uint index = 0; index < boundaryCondition->indices.size(); index++ )
        {
            for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
            {
                qs[dir][allIndicesCounter] = boundaryCondition->qs[index][dir];

                //real x,y,z;
                //grids[level]->transIndexToCoords( index, x, y, z );

                //x += grids[level]->getDirection()[dir * DIMENSION + 0] * grids[level]->getDelta();
                //y += grids[level]->getDirection()[dir * DIMENSION + 1] * grids[level]->getDelta();
                //z += grids[level]->getDirection()[dir * DIMENSION + 2] * grids[level]->getDelta();

                //uint neighborIndex = grids[level]->transCoordToIndex( x, y, z );

                //if( grids[level]->getFieldEntry(neighborIndex) == STOPPER_OUT_OF_GRID_BOUNDARY ||
                //    grids[level]->getFieldEntry(neighborIndex) == STOPPER_OUT_OF_GRID || 
                //    grids[level]->getFieldEntry(neighborIndex) == STOPPER_SOLID)
                //    qs[dir][allIndicesCounter] = 0.5;
                //else
                //    qs[dir][allIndicesCounter] = -1.0;
            }
            allIndicesCounter++;
        }
    }
}

void LevelGridBuilder::getPressureQs(real* qs[27], int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->pressureBoundaryConditions)
    {
        for ( uint index = 0; index < boundaryCondition->indices.size(); index++ )
        {
            for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
            {
                qs[dir][allIndicesCounter] = boundaryCondition->qs[index][dir];

                //real x,y,z;
                //grids[level]->transIndexToCoords( index, x, y, z );

                //x += grids[level]->getDirection()[dir * DIMENSION + 0] * grids[level]->getDelta();
                //y += grids[level]->getDirection()[dir * DIMENSION + 1] * grids[level]->getDelta();
                //z += grids[level]->getDirection()[dir * DIMENSION + 2] * grids[level]->getDelta();

                //uint neighborIndex = grids[level]->transCoordToIndex( x, y, z );

                //if( grids[level]->getFieldEntry(neighborIndex) == STOPPER_OUT_OF_GRID_BOUNDARY ||
                //    grids[level]->getFieldEntry(neighborIndex) == STOPPER_OUT_OF_GRID || 
                //    grids[level]->getFieldEntry(neighborIndex) == STOPPER_SOLID)
                //    qs[dir][allIndicesCounter] = 0.5;
                //else
                //    qs[dir][allIndicesCounter] = -1.0;
            }
            allIndicesCounter++;
        }
    }
}


uint LevelGridBuilder::getGeometrySize(int level) const
{
    if (boundaryConditions[level]->geometryBoundaryCondition)
        return  (uint)boundaryConditions[level]->geometryBoundaryCondition->indices.size();
    
    return 0;
}

void LevelGridBuilder::getGeometryIndices(int* indices, int level) const
{
    for (uint i = 0; i <  boundaryConditions[level]->geometryBoundaryCondition->indices.size(); i++)
    {
        indices[i] = grids[level]->getSparseIndex(boundaryConditions[level]->geometryBoundaryCondition->indices[i]) + 1;
    }
}

bool LevelGridBuilder::hasGeometryValues() const
{
    return geometryHasValues;
}


void LevelGridBuilder::getGeometryValues(real* vx, real* vy, real* vz, int level) const
{
    for (uint i = 0; i < boundaryConditions[level]->geometryBoundaryCondition->indices.size(); i++)
    {
		vx[i] = boundaryConditions[level]->geometryBoundaryCondition->getVx(i);
		vy[i] = boundaryConditions[level]->geometryBoundaryCondition->getVy(i);
		vz[i] = boundaryConditions[level]->geometryBoundaryCondition->getVz(i);
    }
}


void LevelGridBuilder::getGeometryQs(real* qs[27], int level) const
{
    for (std::size_t i = 0; i < boundaryConditions[level]->geometryBoundaryCondition->indices.size(); i++)
    {
        for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
        {
            qs[dir][i] = boundaryConditions[level]->geometryBoundaryCondition->qs[i][dir];
        }
    }
}





void LevelGridBuilder::writeArrows(std::string fileName) const 
{
    QLineWriter::writeArrows(fileName, boundaryConditions[getNumberOfGridLevels() - 1]->geometryBoundaryCondition, grids[getNumberOfGridLevels() - 1]);
}

GRIDGENERATOR_EXPORT SPtr<BoundaryCondition> LevelGridBuilder::getBoundaryCondition(SideType side, uint level) const
{
    for( auto bc : this->boundaryConditions[level]->pressureBoundaryConditions )
        if( bc->isSide(side) )
            return bc;
    
    for( auto bc : this->boundaryConditions[level]->velocityBoundaryConditions )
        if( bc->isSide(side) )
            return bc;

    auto bc = this->boundaryConditions[level]->geometryBoundaryCondition;
    
    if( bc && bc->isSide(side) )
        return bc;

    return nullptr;
}

GRIDGENERATOR_EXPORT SPtr<GeometryBoundaryCondition> LevelGridBuilder::getGeometryBoundaryCondition(uint level) const
{
    return this->boundaryConditions[level]->geometryBoundaryCondition;
}
