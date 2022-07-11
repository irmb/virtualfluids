//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LevelGridBuilder.cpp
//! \ingroup grid
//! \author Soeren Peters, Stephan Lenz, Martin Sch�nherr
//=======================================================================================
#include "LevelGridBuilder.h"

#include <stdio.h>
#include <iostream>
#include <algorithm>

#include "geometries/Arrow/ArrowImp.h"
#include "geometries/BoundingBox/BoundingBox.h"
#include "geometries/Triangle/Triangle.h"
#include "geometries/TriangularMesh/TriangularMesh.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/Grid.h"
#include "grid/GridFactory.h"
#include "grid/GridInterface.h"
#include "grid/NodeValues.h"

#include "io/GridVTKWriter/GridVTKWriter.h"
#include "io/QLineWriter.h"
#include "io/SimulationFileWriter/SimulationFileWriter.h"
#include "VirtualFluids_GPU/PreCollisionInteractor/VelocitySetter.h"

#include "utilities/communication.h"
#include "utilities/transformator/ArrowTransformator.h"

#define GEOFLUID 19
#define GEOSOLID 16

using namespace vf::gpu;

LevelGridBuilder::LevelGridBuilder()
{
    this->communicationProcesses[CommunicationDirections::MX] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::PX] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::MY] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::PY] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::MZ] = INVALID_INDEX;
    this->communicationProcesses[CommunicationDirections::PZ] = INVALID_INDEX;
}

std::shared_ptr<LevelGridBuilder> LevelGridBuilder::makeShared()
{
    return SPtr<LevelGridBuilder>(new LevelGridBuilder());
}

void LevelGridBuilder::setSlipBoundaryCondition(SideType sideType, real normalX, real normalY, real normalZ)
{
    if(sideType == SideType::GEOMETRY){
        setSlipGeometryBoundaryCondition(normalX, normalY, normalZ);
    }else{
        SPtr<SlipBoundaryCondition> slipBoundaryCondition = SlipBoundaryCondition::make(normalX, normalY, normalZ);

        auto side = SideFactory::make(sideType);

        slipBoundaryCondition->side = side;
        slipBoundaryCondition->side->addIndices(grids, 0, slipBoundaryCondition);

        slipBoundaryCondition->fillSlipNormalLists();
        boundaryConditions[0]->slipBoundaryConditions.push_back(slipBoundaryCondition);

        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Slip BC on level " << 0 << " with " << (int)slipBoundaryCondition->indices.size() << "\n";
    }
}

void LevelGridBuilder::setSlipGeometryBoundaryCondition(real normalX, real normalY, real normalZ)
{
    geometryHasValues = true;

    for (uint level = 0; level < getNumberOfGridLevels(); level++)
    {
		if (boundaryConditions[level]->geometryBoundaryCondition != nullptr)
		{
			boundaryConditions[level]->geometryBoundaryCondition->normalX = normalX;
			boundaryConditions[level]->geometryBoundaryCondition->normalY = normalY;
			boundaryConditions[level]->geometryBoundaryCondition->normalZ = normalZ;
			boundaryConditions[level]->geometryBoundaryCondition->side->addIndices(grids, level, boundaryConditions[level]->geometryBoundaryCondition);

            boundaryConditions[level]->geometryBoundaryCondition->fillSlipNormalLists();

            *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Geometry Slip BC on level " << level << " with " << (int)boundaryConditions[level]->geometryBoundaryCondition->indices.size() <<"\n";
		}
    }
}

void LevelGridBuilder::setStressBoundaryCondition(  SideType sideType, 
                                                    real nomalX, real normalY, real normalZ, 
                                                    uint samplingOffset, real z0)
{
    SPtr<StressBoundaryCondition> stressBoundaryCondition = StressBoundaryCondition::make(nomalX, normalY, normalZ, samplingOffset, z0);

    auto side = SideFactory::make(sideType);

    stressBoundaryCondition->side = side;
    stressBoundaryCondition->side->addIndices(grids, 0, stressBoundaryCondition);

    stressBoundaryCondition->fillStressNormalLists();
    stressBoundaryCondition->fillSamplingOffsetLists();
    stressBoundaryCondition->fillZ0Lists();
    // stressBoundaryCondition->fillSamplingIndices(grids, 0, samplingOffset); //redundant with Side::setStressSamplingIndices but potentially a better approach for cases with complex geometries

    boundaryConditions[0]->stressBoundaryConditions.push_back(stressBoundaryCondition);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Stress BC on level " << 0 << " with " << (int)stressBoundaryCondition->indices.size() << "\n";
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
    if (sideType == SideType::GEOMETRY)
        setNoSlipGeometryBoundaryCondition();
    else {
        for (uint level = 0; level < getNumberOfGridLevels(); level++) {
            SPtr<VelocityBoundaryCondition> noSlipBoundaryCondition = VelocityBoundaryCondition::make(0.0, 0.0, 0.0);

            auto side = SideFactory::make(sideType);

            noSlipBoundaryCondition->side = side;
            noSlipBoundaryCondition->side->addIndices(grids, level, noSlipBoundaryCondition);

            noSlipBoundaryCondition->fillVelocityLists();

            // now effectively just a wrapper for velocityBC with zero velocity. No distinction in Gridgenerator.
            boundaryConditions[level]->velocityBoundaryConditions.push_back(noSlipBoundaryCondition); 
        }
    }
}

void LevelGridBuilder::setNoSlipGeometryBoundaryCondition()
{
    geometryHasValues = true;

    for (uint level = 0; level < getNumberOfGridLevels(); level++)
    {
		if (boundaryConditions[level]->geometryBoundaryCondition != nullptr)
		{
			boundaryConditions[level]->geometryBoundaryCondition->side->addIndices(grids, level, boundaryConditions[level]->geometryBoundaryCondition);

            *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Geometry No-Slip BC on level " << level << " with " << (int)boundaryConditions[level]->geometryBoundaryCondition->indices.size() <<"\n";
		}
    }
}

void LevelGridBuilder::setPrecursorBoundaryCondition(SideType sideType, real vx, real vy, real vz, SPtr<VelocityFileCollection> fileCollection, int nTRead)
{
    for (uint level = 0; level < getNumberOfGridLevels(); level++)
    {
        auto reader = fileCollection->createReaderForCollection();
        SPtr<PrecursorBoundaryCondition> precursorBoundaryCondition = PrecursorBoundaryCondition::make(vx, vy, vz, reader, nTRead);

        auto side = SideFactory::make(sideType);

        precursorBoundaryCondition->side = side;
        precursorBoundaryCondition->side->addIndices(grids, level, precursorBoundaryCondition);

        boundaryConditions[level]->precursorBoundaryConditions.push_back(precursorBoundaryCondition);

        *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Precursor BC on level " << level << " with " << (int)precursorBoundaryCondition->indices.size() << "\n";
    }
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::setEnableFixRefinementIntoTheWall(bool enableFixRefinementIntoTheWall)
{
    for (uint level = 0; level < this->grids.size(); level++)
        grids[level]->setEnableFixRefinementIntoTheWall(enableFixRefinementIntoTheWall);
}

GRIDGENERATOR_EXPORT void LevelGridBuilder::setCommunicationProcess(int direction, uint process)
{
    this->communicationProcesses[direction] = process;
}

GRIDGENERATOR_EXPORT uint LevelGridBuilder::getCommunicationProcess(int direction)
{
    return this->communicationProcesses[direction];
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
    for (const auto &grid : grids)
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

uint LevelGridBuilder::getSlipSize(int level) const
{
    uint size = 0;
    for (auto boundaryCondition : boundaryConditions[level]->slipBoundaryConditions)
    {
        size += uint(boundaryCondition->indices.size());
    }
    return size;
}

void LevelGridBuilder::getSlipValues(real* normalX, real* normalY, real* normalZ, int* indices, int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->slipBoundaryConditions)
    {
        for (uint index = 0; index < boundaryCondition->indices.size(); index++)
        {
            indices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->indices[index]) + 1;
            
            normalX[allIndicesCounter] = boundaryCondition->getNormalx(index);
            normalY[allIndicesCounter] = boundaryCondition->getNormaly(index);
            normalZ[allIndicesCounter] = boundaryCondition->getNormalz(index);
            allIndicesCounter++;
        }
    }
}

void LevelGridBuilder::getSlipQs(real* qs[27], int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->slipBoundaryConditions)
    {
        for (uint index = 0; index < boundaryCondition->indices.size(); index++)
        {
            for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
            {
                qs[dir][allIndicesCounter] = boundaryCondition->qs[index][dir];
            }
            allIndicesCounter++;
        }
    }
}

uint LevelGridBuilder::getStressSize(int level) const
{
    uint size = 0;
    for (auto boundaryCondition : boundaryConditions[level]->stressBoundaryConditions)
    {
        size += uint(boundaryCondition->indices.size());
    }
    return size;
}

void LevelGridBuilder::getStressValues( real* normalX, real* normalY, real* normalZ, 
                                        real* vx,      real* vy,      real* vz, 
                                        real* vx1,     real* vy1,     real* vz1, 
                                        int* indices, int* samplingIndices, int* samplingOffset, real* z0, int level) const
{

    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->stressBoundaryConditions)
    {
        for (uint index = 0; index < boundaryCondition->indices.size(); index++)
        {
            indices[allIndicesCounter]          = grids[level]->getSparseIndex(boundaryCondition->indices[index]) + 1;
            samplingIndices[allIndicesCounter]  = grids[level]->getSparseIndex(boundaryCondition->velocitySamplingIndices[index]) + 1;

            normalX[allIndicesCounter] = boundaryCondition->getNormalx(index);
            normalY[allIndicesCounter] = boundaryCondition->getNormaly(index);
            normalZ[allIndicesCounter] = boundaryCondition->getNormalz(index);

            samplingOffset[allIndicesCounter] = boundaryCondition->getSamplingOffset(index);
            z0[allIndicesCounter] = boundaryCondition->getZ0(index);
            allIndicesCounter++;
        }
    }
}

void LevelGridBuilder::getStressQs(real* qs[27], int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->stressBoundaryConditions)
    {
        for (uint index = 0; index < boundaryCondition->indices.size(); index++)
        {
            for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
            {
                qs[dir][allIndicesCounter] = boundaryCondition->qs[index][dir];
            }
            allIndicesCounter++;
        }
    }
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

void LevelGridBuilder::getVelocityValues(real* vx, real* vy, real* vz, int* indices, int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->velocityBoundaryConditions)
    {
        for(std::size_t i = 0; i < boundaryCondition->indices.size(); i++)
        {
            indices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->indices[i]) +1;  

            vx[allIndicesCounter] = boundaryCondition->getVx((uint)i);
            vy[allIndicesCounter] = boundaryCondition->getVy((uint)i);
            vz[allIndicesCounter] = boundaryCondition->getVz((uint)i);
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
            }
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
            }
            allIndicesCounter++;
        }
    }
}

uint LevelGridBuilder::getPrecursorSize(int level) const
{
    uint size = 0;
    for (auto boundaryCondition : boundaryConditions[level]->precursorBoundaryConditions)
    {
        size += uint(boundaryCondition->indices.size());
    }
    return size;
}

void LevelGridBuilder::getPrecursorValues(  uint* neighborNT, uint* neighborNB, uint* neighborST, uint* neighborSB, 
                                            real* weightsNT, real* weightsNB, real* weightsST, real* weightsSB, 
                                            int* indices, std::vector<SPtr<VelocityReader>>& reader, 
                                            int& nVelocityPoints, uint& nTRead, 
                                            real& velocityX, real& velocityY, real& velocityZ, int level) const
{
    int allIndicesCounter = 0;
    int allVelocitiesCounter = 0;
    uint tmpNTRead = 0;

    for (auto boundaryCondition : boundaryConditions[level]->precursorBoundaryConditions)
    {
        if( tmpNTRead == 0 )
            tmpNTRead = boundaryCondition->nTRead;
        if( tmpNTRead != boundaryCondition->nTRead )
            throw std::runtime_error("All precursor boundary conditions must have the same NTRead value");

        auto BCreader = boundaryCondition->getReader();
        BCreader->setWritingOffset(allIndicesCounter);
        reader.push_back(BCreader);

        std::vector<real> y, z;
        real xTmp, yTmp, zTmp;
        for(uint i = 0; i<boundaryCondition->indices.size(); i++)
        {
            indices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->indices[i]) + 1;
            grids[level]->transIndexToCoords(boundaryCondition->indices[i], xTmp, yTmp, zTmp);
            y.push_back(yTmp);
            z.push_back(zTmp);
            allIndicesCounter++;
        }
        BCreader->fillArrays(y, z);
        BCreader->getNeighbors(neighborNT, neighborNB, neighborST, neighborSB);
        BCreader->getWeights(weightsNT, weightsNB, weightsST, weightsSB);
        allVelocitiesCounter += BCreader->getNPointsRead();
        velocityX = boundaryCondition->getVelocityX();
        velocityY = boundaryCondition->getVelocityY();
        velocityZ = boundaryCondition->getVelocityZ();
    }
    nVelocityPoints = allVelocitiesCounter;

    if (tmpNTRead == 0)
        throw std::runtime_error("NTRead of precursor needs to be larger than 0.");

    nTRead = tmpNTRead;
}

void LevelGridBuilder::getPrecursorQs(real* qs[27], int level) const
{
    int allIndicesCounter = 0;
    for (auto boundaryCondition : boundaryConditions[level]->precursorBoundaryConditions)
    {
        for ( uint index = 0; index < boundaryCondition->indices.size(); index++ )
        {
            for (int dir = 0; dir <= grids[level]->getEndDirection(); dir++)
            {
                qs[dir][allIndicesCounter] = boundaryCondition->qs[index][dir];
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

GRIDGENERATOR_EXPORT SPtr<gg::BoundaryCondition> LevelGridBuilder::getBoundaryCondition(SideType side, uint level) const
{
    for (auto bc : this->boundaryConditions[level]->slipBoundaryConditions)
        if (bc->isSide(side))
            return bc;

    for (auto bc : this->boundaryConditions[level]->velocityBoundaryConditions)
        if( bc->isSide(side) )
            return bc;

    for (auto bc : this->boundaryConditions[level]->pressureBoundaryConditions)
        if (bc->isSide(side))
            return bc;

    auto bc = this->boundaryConditions[level]->geometryBoundaryCondition;

    if (bc && bc->isSide(side))
        return bc;

    return nullptr;
}

GRIDGENERATOR_EXPORT SPtr<GeometryBoundaryCondition> LevelGridBuilder::getGeometryBoundaryCondition(uint level) const
{
    return this->boundaryConditions[level]->geometryBoundaryCondition;
}
