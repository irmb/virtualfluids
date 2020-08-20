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
//! \file DataBase.cpp
//! \ingroup DataBase
//! \author Stephan Lenz
//=======================================================================================
#include "DataBase.h"

#include <iostream>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "DataBaseAllocator.h"
#include "DataBaseStruct.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

DataBase::DataBase( std::string type ) 
        : myAllocator    ( DataBaseAllocator::create( type ) ),
          numberOfNodes      (0),
          numberOfCells      (0),
          numberOfFaces      (0),
          numberOfLevels     (0),
          cellToCell     (nullptr),
          faceToCell     (nullptr),
          parentCell     (nullptr),
          faceCenter     (nullptr),
          cellCenter     (nullptr),
          cellProperties (nullptr),
          faceOrientation(nullptr),
          data           (nullptr),
          dataUpdate     (nullptr),
          massFlux       (nullptr),
          diffusivity      (nullptr)
{
}

DataBase::~DataBase()
{
    this->myAllocator->freeMemory( *this );
}

void DataBase::setMesh(GksMeshAdapter & adapter)
{
    this->numberOfNodes      = adapter.nodes.size();

    this->numberOfCells      = adapter.cells.size();

    this->numberOfFaces      = adapter.faces.size();

    this->numberOfLevels     = adapter.numberOfLevels;

    this->perLevelCount.resize( this->numberOfLevels );

    for( uint level = 0; level < this->numberOfLevels; level++ )
    {
        perLevelCount[ level ].numberOfCells = adapter.numberOfCellsPerLevel[ level ];
        perLevelCount[ level ].startOfCells  = adapter.startOfCellsPerLevel [ level ];

        perLevelCount[ level ].numberOfBulkCells = adapter.numberOfBulkCellsPerLevel[ level ];

        perLevelCount[ level ].numberOfFacesX = adapter.numberOfFacesPerLevelXYZ[ 3 * level     ];
        perLevelCount[ level ].startOfFacesX  = adapter.startOfFacesPerLevelXYZ [ 3 * level     ];

        perLevelCount[ level ].numberOfFacesY = adapter.numberOfFacesPerLevelXYZ[ 3 * level + 1 ];
        perLevelCount[ level ].startOfFacesY  = adapter.startOfFacesPerLevelXYZ [ 3 * level + 1 ];

        perLevelCount[ level ].numberOfFacesZ = adapter.numberOfFacesPerLevelXYZ[ 3 * level + 2 ];
        perLevelCount[ level ].startOfFacesZ  = adapter.startOfFacesPerLevelXYZ [ 3 * level + 2 ];

        perLevelCount[ level ].numberOfFaces = perLevelCount[ level ].numberOfFacesX
                                             + perLevelCount[ level ].numberOfFacesY
                                             + perLevelCount[ level ].numberOfFacesZ;

        perLevelCount[ level ].numberOfInnerFaces = adapter.numberOfInnerFacesPerLevel[ level ];
    }

    this->myAllocator->allocateMemory( shared_from_this() );

    this->myAllocator->copyMesh( shared_from_this(), adapter );
}

void DataBase::copyDataHostToDevice()
{
    this->myAllocator->copyDataHostToDevice( shared_from_this() );
}

void DataBase::copyDataDeviceToHost()
{
    this->myAllocator->copyDataDeviceToHost( shared_from_this(), this->dataHost.data() );
}

void DataBase::copyDataDeviceToHost( real* dataHost )
{
    this->myAllocator->copyDataDeviceToHost( shared_from_this(), dataHost );
}

int DataBase::getCrashCellIndex()
{
    return this->myAllocator->getCrashCellIndex(shared_from_this());
}

DataBaseStruct DataBase::toStruct()
{
    DataBaseStruct dataBase;

    dataBase.numberOfCells            = this->numberOfCells;
    dataBase.numberOfFaces            = this->numberOfFaces;

    dataBase.cellToCell               = this->cellToCell;
    dataBase.faceToCell               = this->faceToCell;

    dataBase.parentCell               = this->parentCell;

    dataBase.faceCenter               = this->faceCenter;
    dataBase.cellCenter               = this->cellCenter;

    dataBase.cellProperties           = this->cellProperties;

    dataBase.faceOrientation          = this->faceOrientation;

    dataBase.data                     = this->data;
    dataBase.dataUpdate               = this->dataUpdate;

    dataBase.massFlux                 = this->massFlux;

    dataBase.diffusivity              = this->diffusivity;

    dataBase.crashCellIndex           = this->crashCellIndex;

    return dataBase;
}

uint DataBase::getCellLevel(uint cellIdx)
{
    uint level = 0;

    while( cellIdx >= this->perLevelCount[level].startOfCells
                   + this->perLevelCount[level].numberOfCells ) level++;

    return level;
}

uint DataBase::getFaceLevel(uint faceIdx)
{
    uint level = 0;

    while( faceIdx >= this->perLevelCount[level].startOfFacesX
                   + this->perLevelCount[level].numberOfFaces ) level++;

    return level;
}

Vec3 DataBase::getCellCenter(uint cellIdx)
{
    Vec3 cellCenter;

    for( uint node = 0; node < 8; node++ )
    {
        cellCenter = cellCenter + this->nodeCoordinates[ this->cellToNode[ cellIdx ][ node ] ];
    }

    cellCenter.x /= c8o1;
    cellCenter.y /= c8o1;
    cellCenter.z /= c8o1;

    return cellCenter;
}

bool DataBase::isGhostCell(uint cellIdx)
{
    uint level = this->getCellLevel( cellIdx );

    return ( cellIdx >= this->perLevelCount[ level ].startOfCells + this->perLevelCount[ level ].numberOfBulkCells )
           ||
           ( isCellProperties( this->cellPropertiesHost[cellIdx], CELL_PROPERTIES_FINE_GHOST ) );

}

std::string DataBase::getDeviceType()
{
    return this->myAllocator->getDeviceType();
}
