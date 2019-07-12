#include "Open.h"

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
                                                           const OpenStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const OpenStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Open::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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
                                        const OpenStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const OpenStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];
    uint secondCellIdx = boundaryCondition.secondCells[ startIndex + index ];

    ConservedVariables domainCellData, secondCellData, ghostCellData;
    readCellData ( domainCellIdx, dataBase, domainCellData );
    readCellData ( secondCellIdx, dataBase, secondCellData );

    PrimitiveVariables domainCellPrim = toPrimitiveVariables( domainCellData, parameters.K );
    PrimitiveVariables secondCellPrim = toPrimitiveVariables( secondCellData, parameters.K );
    
    //////////////////////////////////////////////////////////////////////////

    real xGhostCell  = dataBase.cellCenter[ VEC_X(ghostCellIdx, dataBase.numberOfCells) ];
    real yGhostCell  = dataBase.cellCenter[ VEC_Y(ghostCellIdx, dataBase.numberOfCells) ];
    real zGhostCell  = dataBase.cellCenter[ VEC_Z(ghostCellIdx, dataBase.numberOfCells) ];
    
    real xDomainCell = dataBase.cellCenter[ VEC_X(domainCellIdx, dataBase.numberOfCells) ];
    real yDomainCell = dataBase.cellCenter[ VEC_Y(domainCellIdx, dataBase.numberOfCells) ];
    real zDomainCell = dataBase.cellCenter[ VEC_Z(domainCellIdx, dataBase.numberOfCells) ];

    real dx = xGhostCell - xDomainCell;
    real dy = yGhostCell - yDomainCell;
    real dz = zGhostCell - zDomainCell;

    real sign = domainCellPrim.U * dx 
              + domainCellPrim.V * dy 
              + domainCellPrim.W * dz;

    //////////////////////////////////////////////////////////////////////////

    real p1 = c1o2 * domainCellPrim.rho / domainCellPrim.lambda;
    real p2 = c1o2 * secondCellPrim.rho / secondCellPrim.lambda;

    real p0 = c1o2 * boundaryCondition.prim.rho / boundaryCondition.prim.lambda;

    //////////////////////////////////////////////////////////////////////////

    //if( sign > zero )
    //if( p2 > p1 )
    if( p1 > p0 )
    {
        ghostCellData = domainCellData;
        //ghostCellData = two * domainCellData + ( - one ) * secondCellData;

        
        ghostCellData.rhoU = c0o1;
        ghostCellData.rhoV = c0o1;
        ghostCellData.rhoW = c0o1;
    }
    else
    {
        PrimitiveVariables ghostCellPrim  = boundaryCondition.prim;

        ghostCellPrim.U = domainCellPrim.U;
        ghostCellPrim.V = domainCellPrim.V;
        ghostCellPrim.W = domainCellPrim.W;

        //ghostCellPrim.U = p0/p1;
        //ghostCellPrim.V = p0/p1;
        //ghostCellPrim.W = p0/p1;

        //ghostCellPrim.U = two * domainCellPrim.U - secondCellPrim.U;
        //ghostCellPrim.V = two * domainCellPrim.V - secondCellPrim.V;
        //ghostCellPrim.W = two * domainCellPrim.W - secondCellPrim.W;

        real velocity = sqrt( ghostCellPrim.U * ghostCellPrim.U + ghostCellPrim.V * ghostCellPrim.V + ghostCellPrim.W * ghostCellPrim.W );

        if( velocity > boundaryCondition.velocityLimiter  )
        {
            ghostCellPrim.U *= boundaryCondition.velocityLimiter / velocity;
            ghostCellPrim.V *= boundaryCondition.velocityLimiter / velocity;
            ghostCellPrim.W *= boundaryCondition.velocityLimiter / velocity;
        }

        ghostCellData = toConservedVariables(ghostCellPrim, parameters.K);
    }

    //////////////////////////////////////////////////////////////////////////

    //ghostCellData = two * domainCellData + ( - one ) * secondCellData;

    //ghostCellData = domainCellData;

    //////////////////////////////////////////////////////////////////////////

    writeCellData(ghostCellIdx, dataBase, ghostCellData);
}

Open::Open(SPtr<DataBase> dataBase, PrimitiveVariables prim, real velocityLimiter)
    : BoundaryCondition( dataBase )
{
    this->prim = prim;

    this->velocityLimiter = velocityLimiter;
}

bool Open::isWall()
{
    return false;
}

bool Open::secondCellsNeeded()
{
    return true;
}

