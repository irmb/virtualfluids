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
//! \file DataBaseAllocatorGPU.cpp
//! \ingroup DataBase
//! \author Stephan Lenz
//=======================================================================================
#include "DataBaseAllocatorGPU.h"

#include <cstring>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/DataTypes.h"
#include "PointerDefinitions.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"

#include "CellProperties/CellProperties.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

#include "Definitions/MemoryAccessPattern.h"

#include "CudaUtility/CudaUtility.h"

void DataBaseAllocatorGPU::freeMemory( DataBase& dataBase )
{
    dataBase.cellToNode.clear();
    dataBase.faceToNode.clear();

    dataBase.cellPropertiesHost.clear();

    checkCudaErrors( cudaFree ( dataBase.cellToCell ) );

    checkCudaErrors( cudaFree ( dataBase.faceToCell ) );

    checkCudaErrors( cudaFree ( dataBase.parentCell ) );

    checkCudaErrors( cudaFree ( dataBase.faceCenter ) );
    checkCudaErrors( cudaFree ( dataBase.cellCenter ) );

    checkCudaErrors( cudaFree ( dataBase.cellProperties ) );

    checkCudaErrors( cudaFree ( dataBase.faceOrientation ) );

    checkCudaErrors( cudaFree ( dataBase.data ) );
    checkCudaErrors( cudaFree ( dataBase.dataUpdate ) );

    checkCudaErrors( cudaFree ( dataBase.massFlux ) );

    checkCudaErrors( cudaFree ( dataBase.diffusivity ) );

    checkCudaErrors( cudaFree ( dataBase.crashCellIndex ) );

    dataBase.dataHost.clear();
}

void DataBaseAllocatorGPU::allocateMemory(SPtr<DataBase> dataBase)
{
    dataBase->cellToNode.resize( dataBase->numberOfCells );
    dataBase->faceToNode.resize( dataBase->numberOfFaces );

    dataBase->cellPropertiesHost.resize( dataBase->numberOfCells );

    checkCudaErrors( cudaMalloc ( &dataBase->cellToCell, sizeof(uint) * LENGTH_CELL_TO_CELL * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->faceToCell, sizeof(uint) * LENGTH_FACE_TO_CELL * dataBase->numberOfFaces ) );

    checkCudaErrors( cudaMalloc ( &dataBase->parentCell, sizeof(uint) * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->faceCenter, sizeof(real) * LENGTH_VECTOR * dataBase->numberOfFaces ) );
    checkCudaErrors( cudaMalloc ( &dataBase->cellCenter, sizeof(real) * LENGTH_VECTOR * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->cellProperties, sizeof(CellProperties) * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->faceOrientation, sizeof(char) * dataBase->numberOfFaces ) );

    checkCudaErrors( cudaMalloc ( &dataBase->data,       sizeof(real) *            LENGTH_CELL_DATA * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &dataBase->dataUpdate, sizeof(realAccumulator) * LENGTH_CELL_DATA * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->massFlux ,  sizeof(real) * LENGTH_VECTOR    * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->diffusivity,  sizeof(real) * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->crashCellIndex,  sizeof(int) ) );

    dataBase->dataHost.resize( LENGTH_CELL_DATA * dataBase->numberOfCells );

    dataBase->diffusivityHost.resize( dataBase->numberOfCells );
}

void DataBaseAllocatorGPU::copyMesh(SPtr<DataBase> dataBase, GksMeshAdapter & adapter)
{
    dataBase->nodeCoordinates = adapter.nodes;

    //////////////////////////////////////////////////////////////////////////

    std::vector<uint> cellToCellBuffer   ( LENGTH_CELL_TO_CELL * dataBase->numberOfCells );

    std::vector<uint> faceToCellBuffer   ( LENGTH_FACE_TO_CELL * dataBase->numberOfFaces );

    std::vector<uint> parentCellBuffer   ( dataBase->numberOfCells );

    std::vector<real> faceCenterBuffer   ( LENGTH_VECTOR * dataBase->numberOfFaces );
    std::vector<real> cellCenterBuffer   ( LENGTH_VECTOR * dataBase->numberOfCells );

    std::vector<char> faceOrientationBuffer( dataBase->numberOfFaces );

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
            cellToCellBuffer[ CELL_TO_CELL( cellIdx, neighbordx, dataBase->numberOfCells ) ] 
                = adapter.cells[ cellIdx ].cellToCell[ neighbordx ];

        parentCellBuffer[ cellIdx ] = adapter.cells[ cellIdx ].parent;

        cellCenterBuffer[ VEC_X( cellIdx, dataBase->numberOfCells ) ] = adapter.cells[ cellIdx ].cellCenter.x;
        cellCenterBuffer[ VEC_Y( cellIdx, dataBase->numberOfCells ) ] = adapter.cells[ cellIdx ].cellCenter.y;
        cellCenterBuffer[ VEC_Z( cellIdx, dataBase->numberOfCells ) ] = adapter.cells[ cellIdx ].cellCenter.z;

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

        faceToCellBuffer[ NEG_CELL( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].negCell;
        faceToCellBuffer[ POS_CELL( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].posCell;

        faceCenterBuffer[ VEC_X( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].faceCenter.x;
        faceCenterBuffer[ VEC_Y( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].faceCenter.y;
        faceCenterBuffer[ VEC_Z( faceIdx, dataBase->numberOfFaces ) ] = adapter.faces[ faceIdx ].faceCenter.z;

        faceOrientationBuffer[ faceIdx ] = adapter.faces[ faceIdx ].orientation;
    }

    //////////////////////////////////////////////////////////////////////////

    checkCudaErrors( cudaMemcpy ( dataBase->cellToCell,     cellToCellBuffer.data(),     sizeof(uint) * LENGTH_CELL_TO_CELL * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    
    checkCudaErrors( cudaMemcpy ( dataBase->faceToCell,     faceToCellBuffer.data(),     sizeof(uint) * LENGTH_FACE_TO_CELL * dataBase->numberOfFaces, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->parentCell,     parentCellBuffer.data(),     sizeof(uint) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->faceCenter,     faceCenterBuffer.data(),     sizeof(real) * LENGTH_VECTOR * dataBase->numberOfFaces, cudaMemcpyHostToDevice ) );
    checkCudaErrors( cudaMemcpy ( dataBase->cellCenter,     cellCenterBuffer.data(),     sizeof(real) * LENGTH_VECTOR * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->cellProperties, dataBase->cellPropertiesHost.data(), sizeof(CellProperties) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->faceOrientation, faceOrientationBuffer.data(), sizeof(char) * dataBase->numberOfFaces, cudaMemcpyHostToDevice ) );

    //////////////////////////////////////////////////////////////////////////

    checkCudaErrors( cudaMemset( dataBase->crashCellIndex, -1, sizeof(int) ) );

    //////////////////////////////////////////////////////////////////////////
}

void DataBaseAllocatorGPU::copyDataHostToDevice(SPtr<DataBase> dataBase)
{
    checkCudaErrors( cudaMemcpy( dataBase->data, dataBase->dataHost.data(), sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
}

void DataBaseAllocatorGPU::copyDataDeviceToHost(SPtr<DataBase> dataBase,  real* hostData )
{
    checkCudaErrors( cudaMemcpy( hostData, dataBase->data, sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );

    checkCudaErrors( cudaMemcpy( dataBase->diffusivityHost.data(), dataBase->diffusivity, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
}

int DataBaseAllocatorGPU::getCrashCellIndex(SPtr<DataBase> dataBase)
{
    int crashCellIndex;

    checkCudaErrors( cudaMemcpy( &crashCellIndex, dataBase->crashCellIndex, sizeof(int), cudaMemcpyDeviceToHost ) );

    return crashCellIndex;
}

void DataBaseAllocatorGPU::freeMemory(BoundaryCondition& boundaryCondition)
{
    checkCudaErrors( cudaFree ( boundaryCondition.ghostCells  ) );
    checkCudaErrors( cudaFree ( boundaryCondition.domainCells ) );
    checkCudaErrors( cudaFree ( boundaryCondition.secondCells ) );
}

void DataBaseAllocatorGPU::allocateMemory(SPtr<BoundaryCondition> boundaryCondition, std::vector<uint> ghostCells, std::vector<uint> domainCells, std::vector<uint> secondCells)
{
    checkCudaErrors( cudaMalloc ( &boundaryCondition->ghostCells , sizeof(uint) * ghostCells.size()  ) );
    checkCudaErrors( cudaMalloc ( &boundaryCondition->domainCells, sizeof(uint) * domainCells.size() ) );
    checkCudaErrors( cudaMalloc ( &boundaryCondition->secondCells, sizeof(uint) * secondCells.size() ) );

    checkCudaErrors( cudaMemcpy ( boundaryCondition->ghostCells , ghostCells.data() , sizeof(uint) * ghostCells.size() , cudaMemcpyHostToDevice ) );
    checkCudaErrors( cudaMemcpy ( boundaryCondition->domainCells, domainCells.data(), sizeof(uint) * domainCells.size(), cudaMemcpyHostToDevice ) );
    checkCudaErrors( cudaMemcpy ( boundaryCondition->secondCells, secondCells.data(), sizeof(uint) * secondCells.size(), cudaMemcpyHostToDevice ) );
}

std::string DataBaseAllocatorGPU::getDeviceType()
{
    return std::string("GPU");
}
