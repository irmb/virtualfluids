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
//! \file DataBaseAllocatorCPU.cpp
//! \ingroup DataBase
//! \author Stephan Lenz
//=======================================================================================
#include "DataBaseAllocatorCPU.h"

#include <cstring>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"

#include "CellProperties/CellProperties.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

#include "Definitions/MemoryAccessPattern.h"

void DataBaseAllocatorCPU::freeMemory( DataBase& dataBase)
{
    dataBase.cellToNode.clear();
    dataBase.faceToNode.clear();

    dataBase.cellPropertiesHost.clear();

    delete [] dataBase.cellToCell;

    delete [] dataBase.faceToCell;

    delete [] dataBase.parentCell;

    delete [] dataBase.faceCenter;
    delete [] dataBase.cellCenter;

    delete [] dataBase.cellProperties;

    delete [] dataBase.faceOrientation;

    delete [] dataBase.data;
    delete [] dataBase.dataUpdate;

    delete [] dataBase.massFlux;

    delete [] dataBase.diffusivity;

    delete [] dataBase.crashCellIndex;

    dataBase.dataHost.clear();
}

void DataBaseAllocatorCPU::allocateMemory(SPtr<DataBase> dataBase)
{
    dataBase->cellToNode.resize( dataBase->numberOfCells );
    dataBase->faceToNode.resize( dataBase->numberOfFaces );

    dataBase->cellPropertiesHost.resize( dataBase->numberOfCells );

    dataBase->cellToCell = new uint [ LENGTH_CELL_TO_CELL * dataBase->numberOfCells ];

    dataBase->faceToCell = new uint [ LENGTH_FACE_TO_CELL * dataBase->numberOfFaces ];

    dataBase->parentCell = new uint [ dataBase->numberOfCells ];

    dataBase->faceCenter = new real [ LENGTH_VECTOR * dataBase->numberOfFaces ];
    dataBase->cellCenter = new real [ LENGTH_VECTOR * dataBase->numberOfCells ];

    dataBase->cellProperties = new CellProperties [ dataBase->numberOfCells ];

    dataBase->faceOrientation = new char [ dataBase->numberOfFaces ];

    dataBase->data       = new real            [ LENGTH_CELL_DATA * dataBase->numberOfCells ];
    dataBase->dataUpdate = new realAccumulator [ LENGTH_CELL_DATA * dataBase->numberOfCells ];

    dataBase->massFlux   = new real [ LENGTH_VECTOR    * dataBase->numberOfCells ];

    dataBase->diffusivity  = new real [ dataBase->numberOfCells ];

    dataBase->crashCellIndex = new int;

    dataBase->dataHost.resize( LENGTH_CELL_DATA * dataBase->numberOfCells );

    dataBase->diffusivityHost.resize( dataBase->numberOfCells );
}

void DataBaseAllocatorCPU::copyMesh(SPtr<DataBase> dataBase, GksMeshAdapter & adapter)
{
    dataBase->nodeCoordinates = adapter.nodes;

    //////////////////////////////////////////////////////////////////////////

    for( uint cellIdx = 0; cellIdx < dataBase->numberOfCells; cellIdx++ )
    {
        dataBase->cellToNode[ cellIdx ][ 0 ] = adapter.cells[ cellIdx ].cellToNode[ 7 ];
        dataBase->cellToNode[ cellIdx ][ 1 ] = adapter.cells[ cellIdx ].cellToNode[ 3 ];
        dataBase->cellToNode[ cellIdx ][ 2 ] = adapter.cells[ cellIdx ].cellToNode[ 1 ];
        dataBase->cellToNode[ cellIdx ][ 3 ] = adapter.cells[ cellIdx ].cellToNode[ 5 ];
        dataBase->cellToNode[ cellIdx ][ 4 ] = adapter.cells[ cellIdx ].cellToNode[ 6 ];
        dataBase->cellToNode[ cellIdx ][ 5 ] = adapter.cells[ cellIdx ].cellToNode[ 2 ];
        dataBase->cellToNode[ cellIdx ][ 6 ] = adapter.cells[ cellIdx ].cellToNode[ 0 ];
        dataBase->cellToNode[ cellIdx ][ 7 ] = adapter.cells[ cellIdx ].cellToNode[ 4 ];
        
        for( uint neighbordx = 0; neighbordx < LENGTH_CELL_TO_CELL; neighbordx++ )
            dataBase->cellToCell[ CELL_TO_CELL( cellIdx, neighbordx, dataBase->numberOfCells ) ] 
                = adapter.cells[ cellIdx ].cellToCell[ neighbordx ];

        dataBase->parentCell[ cellIdx ] = adapter.cells[ cellIdx ].parent;

        dataBase->cellCenter[ VEC_X( cellIdx, dataBase->numberOfCells ) ] = adapter.cells[ cellIdx ].cellCenter.x;
        dataBase->cellCenter[ VEC_Y( cellIdx, dataBase->numberOfCells ) ] = adapter.cells[ cellIdx ].cellCenter.y;
        dataBase->cellCenter[ VEC_Z( cellIdx, dataBase->numberOfCells ) ] = adapter.cells[ cellIdx ].cellCenter.z;

        dataBase->cellPropertiesHost[ cellIdx ] = CELL_PROPERTIES_DEFAULT;

        if( adapter.cells[ cellIdx ].isWall )
            setCellProperties( dataBase->cellPropertiesHost[ cellIdx ], CELL_PROPERTIES_WALL ); 

        if( adapter.cells[ cellIdx ].isFluxBC )
            setCellProperties( dataBase->cellPropertiesHost[ cellIdx ], CELL_PROPERTIES_IS_FLUX_BC );

        if( adapter.cells[ cellIdx ].isInsulated )
            setCellProperties( dataBase->cellPropertiesHost[ cellIdx ], CELL_PROPERTIES_IS_INSULATED ); 

        if( adapter.cells[ cellIdx ].isGhostCell )
            setCellProperties( dataBase->cellPropertiesHost[ cellIdx ], CELL_PROPERTIES_GHOST ); 

        if( adapter.cells[ cellIdx ].isFineGhostCell() )
            setCellProperties( dataBase->cellPropertiesHost[ cellIdx ], CELL_PROPERTIES_FINE_GHOST ); 
    }

    for( uint faceIdx = 0; faceIdx < dataBase->numberOfFaces; faceIdx++ )
    {
        for( uint nodeIdx = 0; nodeIdx < 4; nodeIdx++ )
            dataBase->faceToNode[ faceIdx ][ nodeIdx ]
                = adapter.faces[ faceIdx ].faceToNode[ nodeIdx ];

        dataBase->faceToCell[ NEG_CELL( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].negCell;
        dataBase->faceToCell[ POS_CELL( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].posCell;

        dataBase->faceCenter[ VEC_X( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].faceCenter.x;
        dataBase->faceCenter[ VEC_Y( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].faceCenter.y;
        dataBase->faceCenter[ VEC_Z( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].faceCenter.z;

        dataBase->faceOrientation[ faceIdx ] = adapter.faces[ faceIdx ].orientation;
    }

    //////////////////////////////////////////////////////////////////////////

    memcpy ( dataBase->cellProperties, dataBase->cellPropertiesHost.data(), sizeof(CellProperties) * dataBase->numberOfCells );

    //////////////////////////////////////////////////////////////////////////

    *dataBase->crashCellIndex = -1;
}

void DataBaseAllocatorCPU::copyDataHostToDevice(SPtr<DataBase> dataBase)
{
    memcpy( dataBase->data, dataBase->dataHost.data(), sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells );
}

void DataBaseAllocatorCPU::copyDataDeviceToHost(SPtr<DataBase> dataBase, real* hostData)
{
    memcpy( hostData, dataBase->data, sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells );
    
    memcpy( dataBase->diffusivityHost.data(), dataBase->diffusivity, sizeof(real) * dataBase->numberOfCells );
}

int DataBaseAllocatorCPU::getCrashCellIndex(SPtr<DataBase> dataBase)
{
    return *dataBase->crashCellIndex;
}

void DataBaseAllocatorCPU::freeMemory(BoundaryCondition& boundaryCondition)
{
    delete [] boundaryCondition.ghostCells ;
    delete [] boundaryCondition.domainCells;
    delete [] boundaryCondition.secondCells;
}

void DataBaseAllocatorCPU::allocateMemory(SPtr<BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells)
{
    boundaryCondition->ghostCells  = new uint[ ghostCells.size()  ];
    boundaryCondition->domainCells = new uint[ domainCells.size() ];
    boundaryCondition->secondCells = new uint[ secondCells.size() ];

    memcpy ( boundaryCondition->ghostCells , ghostCells.data() , sizeof(uint) * ghostCells.size()  );
    memcpy ( boundaryCondition->domainCells, domainCells.data(), sizeof(uint) * domainCells.size() );
    memcpy ( boundaryCondition->secondCells, secondCells.data(), sizeof(uint) * secondCells.size() );
}

std::string DataBaseAllocatorCPU::getDeviceType()
{
    return std::string("CPU");
}
