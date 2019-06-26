#include "Symmetry.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const SymmetryStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const SymmetryStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Symmetry::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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

    getLastCudaError("IsothermalWall::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const SymmetryStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const SymmetryStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];

    //////////////////////////////////////////////////////////////////////////

    ConservedVariables domainCellData;
    readCellData ( domainCellIdx, dataBase, domainCellData );

    //////////////////////////////////////////////////////////////////////////
    
    PrimitiveVariables domainCellPrim = toPrimitiveVariables( domainCellData, parameters.K );
    PrimitiveVariables ghostCellPrim  = toPrimitiveVariables( domainCellData, parameters.K );

    //////////////////////////////////////////////////////////////////////////

    if( boundaryCondition.direction == 'x' ) ghostCellPrim.U = - domainCellPrim.U;
    if( boundaryCondition.direction == 'y' ) ghostCellPrim.V = - domainCellPrim.V;
    if( boundaryCondition.direction == 'z' ) ghostCellPrim.W = - domainCellPrim.W;

    //////////////////////////////////////////////////////////////////////////

    ConservedVariables ghostCellData = toConservedVariables(ghostCellPrim, parameters.K);

    writeCellData( ghostCellIdx , dataBase, ghostCellData );
}

Symmetry::Symmetry(SPtr<DataBase> dataBase, char direction)
    : BoundaryCondition( dataBase )
{
    this->direction = direction;
}

bool Symmetry::isWall()
{
    return true;
}

bool Symmetry::secondCellsNeeded()
{
    return false;
}

