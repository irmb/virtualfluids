#include "IsothermalWall.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu{

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const IsothermalWallStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const IsothermalWallStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void IsothermalWall::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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

    getLastCudaError("IsothermalWall::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const IsothermalWallStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const IsothermalWallStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];
    uint secondCellIdx = boundaryCondition.secondCells[ startIndex + index ];

    PrimitiveVariables ghostCellPrim;
    {
        PrimitiveVariables domainCellPrim;
        PrimitiveVariables secondCellPrim;

        {
            ConservedVariables domainCellData;
            readCellData( domainCellIdx, dataBase, domainCellData );
            domainCellPrim = toPrimitiveVariables( domainCellData, parameters.K );

            ConservedVariables secondCellData;
            if( secondCellIdx != INVALID_INDEX ){
                readCellData( secondCellIdx, dataBase, secondCellData );
                secondCellPrim = toPrimitiveVariables( secondCellData, parameters.K );
            }
        }

        ghostCellPrim.U      = c2o1 * boundaryCondition.velocity.x - domainCellPrim.U;
        ghostCellPrim.V      = c2o1 * boundaryCondition.velocity.y - domainCellPrim.V;
        ghostCellPrim.W      = c2o1 * boundaryCondition.velocity.z - domainCellPrim.W;
        ghostCellPrim.lambda = c2o1 * boundaryCondition.lambda     - domainCellPrim.lambda;
    #ifdef USE_PASSIVE_SCALAR
        ghostCellPrim.S_1    = /*two * boundaryCondition.S_1 -*/ domainCellPrim.S_1;
        ghostCellPrim.S_2    = /*two * boundaryCondition.S_2 -*/ domainCellPrim.S_2;
    #endif // USE_PASSIVE_SCALAR


        if( boundaryCondition.useSecondCells && secondCellIdx != INVALID_INDEX ){
            real p1 = c1o2 * domainCellPrim.rho / domainCellPrim.lambda;
            real p2 = c1o2 * secondCellPrim.rho / secondCellPrim.lambda;

            ghostCellPrim.rho = c2o1 * ( c2o1 * p1 - p2 ) * ghostCellPrim.lambda;
        }
        else{
            real p = c1o2 * domainCellPrim.rho / domainCellPrim.lambda;

            ghostCellPrim.rho = c2o1 * p * ghostCellPrim.lambda;
        }
    }

    {
        ConservedVariables ghostCons = toConservedVariables( ghostCellPrim, parameters.K );

        writeCellData( ghostCellIdx, dataBase, ghostCons );
    }
}

IsothermalWall::IsothermalWall(SPtr<DataBase> dataBase, Vec3 velocity, real lambda, bool useSecondCells, real S_1, real S_2)
    : BoundaryCondition( dataBase )
{
    this->velocity       = velocity;
    this->lambda         = lambda;
    this->S_1            = S_1;
    this->S_2            = S_2;
    this->useSecondCells = useSecondCells;
}

bool IsothermalWall::isWall()
{
    return true;
}

bool IsothermalWall::secondCellsNeeded()
{
    return true;
}

} // namespace GksGpu

