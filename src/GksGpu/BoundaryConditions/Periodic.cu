#include "Periodic.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "GksMeshAdapter/GksMeshAdapter.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseStruct.h"
#include "DataBase/DataBaseAllocator.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const BoundaryConditionStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const BoundaryConditionStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Periodic::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
                                          const Parameters parameters, 
                                          const uint level)
{    
    CudaUtility::CudaGrid grid( this->numberOfCellsPerLevel[ level ], 32 );

    runKernel( boundaryConditionKernel,
               boundaryConditionFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               this->toStruct(),
               parameters,
               this->startOfCellsPerLevel[ level ] );

    cudaDeviceSynchronize();

    getLastCudaError("Periodic::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const BoundaryConditionStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const BoundaryConditionStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];
    
    ConservedVariables domainCellData;
    readCellData ( domainCellIdx, dataBase, domainCellData );
    writeCellData( ghostCellIdx , dataBase, domainCellData );
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Periodic::findBoundaryCells(GksMeshAdapter & adapter, bool allowGhostCells, std::function<bool(Vec3)> boundaryFinder)
{
    this->myAllocator->freeMemory( *this );

    std::vector<uint> ghostCells;
    std::vector<uint> domainCells;
    std::vector<uint> secondCells;

    numberOfCellsPerLevel.resize( adapter.numberOfLevels );
    startOfCellsPerLevel.resize ( adapter.numberOfLevels );

    for( auto& n : numberOfCellsPerLevel ) n = 0;

    for( uint level = 0; level < adapter.numberOfLevels; level++ )
    {
        uint startIdx = adapter.startOfCellsPerLevel[level] 
                      + adapter.numberOfBulkCellsPerLevel[level];

        uint endIdx   = adapter.startOfCellsPerLevel[level] 
                      + adapter.numberOfCellsPerLevel[level];

        for( uint_2 candidate : adapter.periodicBoundaryNeighbors )
        {
            MeshCell& cell = adapter.cells[ candidate[0] ];

            if( !boundaryFinder( cell.cellCenter ) ) continue;
         
            if( candidate[1] == INVALID_INDEX ) continue;
            
            ghostCells.push_back ( candidate[0] );
            domainCells.push_back( candidate[1] );
                
            this->numberOfCellsPerLevel[ level ]++;
        }
    }

    startOfCellsPerLevel[ 0 ] = 0;

    for( uint level = 1; level < adapter.numberOfLevels; level++ )
    {
        startOfCellsPerLevel[ level ] = startOfCellsPerLevel [ level - 1 ]
                                      + numberOfCellsPerLevel[ level - 1 ];
    }

    this->numberOfCells = ghostCells.size();

    this->myAllocator->allocateMemory( shared_from_this(), ghostCells, domainCells, secondCells );

}

bool Periodic::isWall()
{
    return false;
}
