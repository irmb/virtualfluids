#include "SalinasVazquez.h"

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
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const SalinasVazquezStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const SalinasVazquezStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void SalinasVazquez::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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
                                        const SalinasVazquezStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const SalinasVazquezStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];
    uint secondCellIdx = boundaryCondition.secondCells[ startIndex + index ];

    real lambda;
    {
        real x = dataBase.cellCenter[ VEC_X(ghostCellIdx, dataBase.numberOfCells) ];

        real TMX = one / boundaryCondition.lambdaMX;
        real TPX = one / boundaryCondition.lambdaPX;

        real T = TPX + ( TMX - TPX ) * ( boundaryCondition.a0 
                                       + boundaryCondition.a1*x 
                                       + boundaryCondition.a2*x*x 
                                       + boundaryCondition.a3*x*x*x );

        lambda = one / T;
    }

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

        ghostCellPrim.U      =              - domainCellPrim.U;
        ghostCellPrim.V      =              - domainCellPrim.V;
        ghostCellPrim.W      =              - domainCellPrim.W;
        ghostCellPrim.lambda = two * lambda - domainCellPrim.lambda;
    #ifdef USE_PASSIVE_SCALAR
        ghostCellPrim.S_1    =                domainCellPrim.S_1;
        ghostCellPrim.S_2    =                domainCellPrim.S_2;
    #endif // USE_PASSIVE_SCALAR


        if( boundaryCondition.useSecondCells && secondCellIdx != INVALID_INDEX ){
            real p1 = c1o2 * domainCellPrim.rho / domainCellPrim.lambda;
            real p2 = c1o2 * secondCellPrim.rho / secondCellPrim.lambda;

            ghostCellPrim.rho = two * ( two * p1 - p2 ) * ghostCellPrim.lambda;
        }
        else{
            real p = c1o2 * domainCellPrim.rho / domainCellPrim.lambda;

            ghostCellPrim.rho = two * p * ghostCellPrim.lambda;
        }
    }

    {
        ConservedVariables ghostCons = toConservedVariables( ghostCellPrim, parameters.K );

        writeCellData( ghostCellIdx, dataBase, ghostCons );
    }
}

SalinasVazquez::SalinasVazquez(SPtr<DataBase> dataBase, real lambdaMX, real lambdaPX, real a0, real a1, real a2, real a3, bool useSecondCells)
    : BoundaryCondition( dataBase )
{
    this->lambdaMX       = lambdaMX;
    this->lambdaPX       = lambdaPX;

    this->a0             = a0;
    this->a1             = a1;
    this->a2             = a2;
    this->a3             = a3;

    this->useSecondCells = useSecondCells;
}

bool SalinasVazquez::isWall()
{
    return true;
}

bool SalinasVazquez::secondCellsNeeded()
{
    return true;
}

