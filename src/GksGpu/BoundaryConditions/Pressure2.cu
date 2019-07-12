#include "Pressure2.h"

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
                                                           const Pressure2Struct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const Pressure2Struct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Pressure2::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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

    getLastCudaError("Pressure::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const Pressure2Struct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const Pressure2Struct& boundaryCondition, 
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

        //ghostCellPrim.rho    = two * domainCellPrim.rho    - secondCellPrim.rho;
        ghostCellPrim.U      = c2o1 * domainCellPrim.U      - secondCellPrim.U;
        ghostCellPrim.V      = c2o1 * domainCellPrim.V      - secondCellPrim.V;
        ghostCellPrim.W      = c2o1 * domainCellPrim.W      - secondCellPrim.W;
        //ghostCellPrim.lambda = two * domainCellPrim.lambda - secondCellPrim.lambda;
        ghostCellPrim.lambda = domainCellPrim.lambda;
    #ifdef USE_PASSIVE_SCALAR
        //ghostCellPrim.S_1    = two * domainCellPrim.S_1    - secondCellPrim.S_1;
        //ghostCellPrim.S_2    = two * domainCellPrim.S_2    - secondCellPrim.S_2;
        ghostCellPrim.S_1    = domainCellPrim.S_1;
        ghostCellPrim.S_2    = domainCellPrim.S_2;
        //ghostCellPrim.S_1    = zero;
        //ghostCellPrim.S_2    = zero;
    #endif // USE_PASSIVE_SCALAR


        real rho0 = ( c2o1 * boundaryCondition.p0 * c1o2 * ( domainCellPrim.lambda + ghostCellPrim.lambda ) );
        ghostCellPrim.rho = c2o1 * rho0 - domainCellPrim.rho;

        //real lambda0 = ( c1o2 * ( domainCellPrim.rho + ghostCellPrim.rho  ) * c1o2 / boundaryCondition.p0 );
        //ghostCellPrim.lambda = two * lambda0 - domainCellPrim.lambda;
    
        //////////////////////////////////////////////////////////////////////////

        real xGhostCell = dataBase.cellCenter[VEC_X(ghostCellIdx, dataBase.numberOfCells)];
        real yGhostCell = dataBase.cellCenter[VEC_Y(ghostCellIdx, dataBase.numberOfCells)];
        real zGhostCell = dataBase.cellCenter[VEC_Z(ghostCellIdx, dataBase.numberOfCells)];

        real xDomainCell = dataBase.cellCenter[VEC_X(domainCellIdx, dataBase.numberOfCells)];
        real yDomainCell = dataBase.cellCenter[VEC_Y(domainCellIdx, dataBase.numberOfCells)];
        real zDomainCell = dataBase.cellCenter[VEC_Z(domainCellIdx, dataBase.numberOfCells)];

        real dx = xGhostCell - xDomainCell;
        real dy = yGhostCell - yDomainCell;
        real dz = zGhostCell - zDomainCell;

        real sign = domainCellPrim.U * dx
                  + domainCellPrim.V * dy
                  + domainCellPrim.W * dz;

        //////////////////////////////////////////////////////////////////////////

        if( sign < c0o1 )
        {
            //ghostCellPrim.U = - domainCellPrim.U;
            //ghostCellPrim.V = - domainCellPrim.V;
            //ghostCellPrim.W = - domainCellPrim.W;
            ghostCellPrim.U = c0o1;
            ghostCellPrim.V = c0o1;
            ghostCellPrim.W = c0o1;
        }
    }

    //////////////////////////////////////////////////////////////////////////


    {
        ConservedVariables ghostCons = toConservedVariables( ghostCellPrim, parameters.K );

        writeCellData( ghostCellIdx, dataBase, ghostCons );
    }
}

Pressure2::Pressure2(SPtr<DataBase> dataBase, real p0)
    : BoundaryCondition( dataBase )
{
    this->p0 = p0;
}

bool Pressure2::isWall()
{
    return false;
}

bool Pressure2::secondCellsNeeded()
{
    return true;
}

