#include "LevelGridBuilder.h"

#include <stdio.h>
#include <iostream>

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

            vx[allIndicesCounter] = boundaryCondition->getVx(i);
            vy[allIndicesCounter] = boundaryCondition->getVy(i);
            vz[allIndicesCounter] = boundaryCondition->getVz(i);
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
        return  boundaryConditions[level]->geometryBoundaryCondition->indices.size();
    
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
