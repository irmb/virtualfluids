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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \author Soeren Peters, Stephan Lenz
//=======================================================================================
#include "LevelGridBuilder.h"

#include <stdio.h>
#include <iostream>

#include "geometries/BoundingBox/BoundingBox.h"

#include "grid/BoundaryConditions/BoundaryCondition.h"
#include "grid/BoundaryConditions/Side.h"
#include "grid/GridStrategy/GridCpuStrategy/GridCpuStrategy.h"
#include "grid/NodeValues.h"
#include "grid/GridFactory.h"
#include "grid/Grid.h"

#define GEOFLUID 19
#define GEOSOLID 16

LevelGridBuilder::LevelGridBuilder(Device device, const std::string& d3qxx) : device(device), d3qxx(d3qxx)
{
}

std::shared_ptr<LevelGridBuilder> LevelGridBuilder::makeShared(Device device, const std::string& d3qxx)
{
    return SPtr<LevelGridBuilder>(new LevelGridBuilder(device, d3qxx));
}

void LevelGridBuilder::setVelocityBoundaryCondition(SideType sideType, real vx, real vy, real vz)
{
    SPtr<VelocityBoundaryCondition> velocityBoundaryCondition = VelocityBoundaryCondition::make(vx, vy, vz);

    auto side = SideFactory::make(sideType);

    velocityBoundaryCondition->side = side;
    velocityBoundaryCondition->side->addIndices(grids, 0, velocityBoundaryCondition);

    velocityBoundaryCondition->fillVelocityLists();

    boundaryConditions[0]->velocityBoundaryConditions.push_back(velocityBoundaryCondition);

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Set Velocity BC on level " << 0 << " with " << (int)velocityBoundaryCondition->indices.size() <<"\n";
}

void LevelGridBuilder::setPeriodicBoundaryCondition(bool periodic_X, bool periodic_Y, bool periodic_Z)
{
    for( uint level = 0; level < this->grids.size(); level++ )
        grids[level]->setPeriodicity(periodic_X, periodic_Y, periodic_Z);
}

void LevelGridBuilder::setNoSlipBoundaryCondition(SideType sideType)
{
    for (int level = 0; level < getNumberOfGridLevels(); level++)
    {
        SPtr<VelocityBoundaryCondition> noSlipBoundaryCondition = VelocityBoundaryCondition::make(0.0, 0.0, 0.0);

        auto side = SideFactory::make(sideType);

        noSlipBoundaryCondition->side = side;
        noSlipBoundaryCondition->side->addIndices(grids, level, noSlipBoundaryCondition);

        boundaryConditions[level]->noSlipBoundaryConditions.push_back(noSlipBoundaryCondition);
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
    return 0;
}

uint LevelGridBuilder::getNumberOfNodesFC(int level)
{
    return 0;
}

void LevelGridBuilder::getGridInterfaceIndices(uint* iCellCfc, uint* iCellCff, uint* iCellFcc, uint* iCellFcf, int level) const
{
}

void LevelGridBuilder::getOffsetFC(real * xOffFC, real * yOffFC, real * zOffFC, int level)
{
}

void LevelGridBuilder::getOffsetCF(real * xOffCF, real * yOffCF, real * zOffCF, int level)
{
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

void LevelGridBuilder::getNodeValues(real *xCoords, real *yCoords, real *zCoords, 
                                     uint *neighborX, uint *neighborY, uint *neighborZ, uint *neighborNegative, 
                                     uint *geo, const int level) const
{
    grids[level]->getNodeValues(xCoords, yCoords, zCoords, neighborX, neighborY, neighborZ, neighborNegative, geo);
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
        for(int i = 0; i < boundaryCondition->indices.size(); i++)
        {
            indices[allIndicesCounter] = grids[level]->getSparseIndex(boundaryCondition->indices[i]) +1;  

            vx[allIndicesCounter] = boundaryCondition->getVx(i);
            vy[allIndicesCounter] = boundaryCondition->getVy(i);
            vz[allIndicesCounter] = boundaryCondition->getVz(i);
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

GRIDGENERATOR_EXPORT SPtr<BoundaryCondition> LevelGridBuilder::getBoundaryCondition(SideType side, uint level) const
{
    for( auto bc : this->boundaryConditions[level]->velocityBoundaryConditions )
        if( bc->isSide(side) )
            return bc;

    return nullptr;
}

