#include "DataBaseAllocatorGPU.h"

#include <cstring>
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/DataTypes.h"
#include "Core/PointerDefinitions.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"

#include "CellProperties/CellProperties.cuh"

#include "BoundaryConditions/BoundaryCondition.h"

#include "Communication/Communicator.h"

#include "Definitions/MemoryAccessPattern.h"

void DataBaseAllocatorGPU::freeMemory( DataBase& dataBase )
{
    dataBase.cellToNode.clear();
    dataBase.faceToNode.clear();

    checkCudaErrors( cudaFree ( dataBase.cellToCell ) );

    checkCudaErrors( cudaFree ( dataBase.faceToCell ) );

    checkCudaErrors( cudaFree ( dataBase.parentCell ) );

    checkCudaErrors( cudaFree ( dataBase.faceCenter ) );
    checkCudaErrors( cudaFree ( dataBase.cellCenter ) );

    checkCudaErrors( cudaFree ( dataBase.cellProperties ) );

    checkCudaErrors( cudaFree ( dataBase.faceOrientation ) );

    checkCudaErrors( cudaFree ( dataBase.fineToCoarse ) );
    checkCudaErrors( cudaFree ( dataBase.coarseToFine ) );

    checkCudaErrors( cudaFree ( dataBase.data ) );
    checkCudaErrors( cudaFree ( dataBase.dataUpdate ) );

    checkCudaErrors( cudaFree ( dataBase.massFlux ) );

    dataBase.dataHost.clear();
}

void DataBaseAllocatorGPU::allocateMemory(SPtr<DataBase> dataBase)
{
    dataBase->cellToNode.resize( dataBase->numberOfCells );
    dataBase->faceToNode.resize( dataBase->numberOfFaces );

    checkCudaErrors( cudaMalloc ( &dataBase->cellToCell, sizeof(uint) * LENGTH_CELL_TO_CELL * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->faceToCell, sizeof(uint) * LENGTH_FACE_TO_CELL * dataBase->numberOfFaces ) );

    checkCudaErrors( cudaMalloc ( &dataBase->parentCell, sizeof(uint) * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->faceCenter, sizeof(real) * LENGTH_VECTOR * dataBase->numberOfFaces ) );
    checkCudaErrors( cudaMalloc ( &dataBase->cellCenter, sizeof(real) * LENGTH_VECTOR * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->cellProperties, sizeof(CellProperties) * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->faceOrientation, sizeof(char) * dataBase->numberOfFaces ) );

    checkCudaErrors( cudaMalloc ( &dataBase->fineToCoarse, sizeof(uint) * LENGTH_FINE_TO_COARSE * dataBase->numberOfCoarseGhostCells ) );
    checkCudaErrors( cudaMalloc ( &dataBase->coarseToFine, sizeof(uint) * LENGTH_COARSE_TO_FINE * dataBase->numberOfFineGhostCells   ) );

    checkCudaErrors( cudaMalloc ( &dataBase->data,       sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &dataBase->dataUpdate, sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells ) );

    checkCudaErrors( cudaMalloc ( &dataBase->massFlux ,  sizeof(real) * LENGTH_VECTOR    * dataBase->numberOfCells ) );

    dataBase->dataHost.resize( LENGTH_CELL_DATA * dataBase->numberOfCells );
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

    std::vector<CellProperties> cellPropertiesBuffer ( dataBase->numberOfCells );

    std::vector<char> faceOrientationBuffer( dataBase->numberOfFaces );

    std::vector<uint> fineToCoarseBuffer ( LENGTH_FINE_TO_COARSE * dataBase->numberOfCoarseGhostCells );
    std::vector<uint> coarseToFineBuffer ( LENGTH_COARSE_TO_FINE * dataBase->numberOfFineGhostCells   );

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

        cellPropertiesBuffer[ cellIdx ] = CELL_PROPERTIES_DEFAULT;
        if( adapter.cells[ cellIdx ].isWall )
            setCellProperties( cellPropertiesBuffer[ cellIdx ], CELL_PROPERTIES_WALL ); 
        if( adapter.cells[ cellIdx ].isGhostCell )
            setCellProperties( cellPropertiesBuffer[ cellIdx ], CELL_PROPERTIES_GHOST ); 

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

        faceOrientationBuffer[ faceIdx ] = adapter.faces[ faceIdx ].orientation;
    }

    //////////////////////////////////////////////////////////////////////////

    for( uint cellIdx = 0; cellIdx < dataBase->numberOfCoarseGhostCells; cellIdx++ ){
        for( uint connectivityIdx = 0; connectivityIdx < LENGTH_FINE_TO_COARSE; connectivityIdx++ ){
            fineToCoarseBuffer[ FINE_TO_COARSE( cellIdx, connectivityIdx, dataBase->numberOfCoarseGhostCells ) ]
                = adapter.fineToCoarse[cellIdx][connectivityIdx];
        }
    }

    for( uint cellIdx = 0; cellIdx < dataBase->numberOfFineGhostCells; cellIdx++ ){
        for( uint connectivityIdx = 0; connectivityIdx < LENGTH_COARSE_TO_FINE; connectivityIdx++ ){
            coarseToFineBuffer[ COARSE_TO_FINE( cellIdx, connectivityIdx, dataBase->numberOfFineGhostCells ) ]
                = adapter.coarseToFine[cellIdx][connectivityIdx];
        }
    }

    //////////////////////////////////////////////////////////////////////////

    checkCudaErrors( cudaMemcpy ( dataBase->cellToCell,     cellToCellBuffer.data(),     sizeof(uint) * LENGTH_CELL_TO_CELL * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    
    checkCudaErrors( cudaMemcpy ( dataBase->faceToCell,     faceToCellBuffer.data(),     sizeof(uint) * LENGTH_FACE_TO_CELL * dataBase->numberOfFaces, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->parentCell,     parentCellBuffer.data(),     sizeof(uint) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->faceCenter,     faceCenterBuffer.data(),     sizeof(real) * LENGTH_VECTOR * dataBase->numberOfFaces, cudaMemcpyHostToDevice ) );
    checkCudaErrors( cudaMemcpy ( dataBase->cellCenter,     cellCenterBuffer.data(),     sizeof(real) * LENGTH_VECTOR * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->cellProperties, cellPropertiesBuffer.data(), sizeof(CellProperties) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->faceOrientation, faceOrientationBuffer.data(), sizeof(char) * dataBase->numberOfFaces, cudaMemcpyHostToDevice ) );

    checkCudaErrors( cudaMemcpy ( dataBase->fineToCoarse,   fineToCoarseBuffer.data(),   sizeof(uint) * LENGTH_FINE_TO_COARSE * dataBase->numberOfCoarseGhostCells, cudaMemcpyHostToDevice ) );
    checkCudaErrors( cudaMemcpy ( dataBase->coarseToFine,   coarseToFineBuffer.data(),   sizeof(uint) * LENGTH_COARSE_TO_FINE * dataBase->numberOfFineGhostCells  , cudaMemcpyHostToDevice ) );

    //////////////////////////////////////////////////////////////////////////
}

void DataBaseAllocatorGPU::copyDataHostToDevice(SPtr<DataBase> dataBase)
{
    checkCudaErrors( cudaMemcpy( dataBase->data, dataBase->dataHost.data(), sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
}

void DataBaseAllocatorGPU::copyDataDeviceToHost(SPtr<DataBase> dataBase,  real* hostData )
{
    checkCudaErrors( cudaMemcpy( hostData, dataBase->data, sizeof(real) * LENGTH_CELL_DATA * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
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

void DataBaseAllocatorGPU::freeMemory(Communicator & communicator)
{
    checkCudaErrors( cudaFree ( communicator.sendIndices ) );
    checkCudaErrors( cudaFree ( communicator.recvIndices ) );

    checkCudaErrors( cudaFree ( communicator.sendBuffer  ) );
    checkCudaErrors( cudaFree ( communicator.recvBuffer  ) );
}

void DataBaseAllocatorGPU::allocateMemory(Communicator & communicator, std::vector<uint>& sendIndices, std::vector<uint>& recvIndices)
{
    checkCudaErrors( cudaMalloc ( &communicator.sendIndices , sizeof(uint) * communicator.numberOfSendNodes ) );
    checkCudaErrors( cudaMalloc ( &communicator.recvIndices , sizeof(uint) * communicator.numberOfRecvNodes ) );
    
    checkCudaErrors( cudaMalloc ( &communicator.sendBuffer , LENGTH_CELL_DATA * sizeof(uint) * communicator.numberOfSendNodes ) );
    checkCudaErrors( cudaMalloc ( &communicator.recvBuffer , LENGTH_CELL_DATA * sizeof(uint) * communicator.numberOfRecvNodes ) );

    checkCudaErrors( cudaMemcpy ( communicator.sendIndices , sendIndices.data() , sizeof(uint) * communicator.numberOfSendNodes, cudaMemcpyHostToDevice ) );
    checkCudaErrors( cudaMemcpy ( communicator.recvIndices , recvIndices.data() , sizeof(uint) * communicator.numberOfRecvNodes, cudaMemcpyHostToDevice ) );
}

void DataBaseAllocatorGPU::copyDataDeviceToDevice(SPtr<Communicator> dst, SPtr<Communicator> src)
{
    checkCudaErrors( cudaMemcpy ( dst->recvBuffer, src->sendBuffer, LENGTH_CELL_DATA * sizeof(real) * src->numberOfSendNodes, cudaMemcpyDefault ) );
}

std::string DataBaseAllocatorGPU::getDeviceType()
{
    return std::string("GPU");
}
