#include "PassiveScalarDiriclet.h"

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
                                                           const PassiveScalarDiricletStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const PassiveScalarDiricletStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void PassiveScalarDiriclet::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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
                                        const PassiveScalarDiricletStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const PassiveScalarDiricletStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
#ifdef USE_PASSIVE_SCALAR
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //real ghostCell_rhoS_1;
    //real ghostCell_rhoS_2;

    //PrimitiveVariables ghostCellPrim;
    //{
    //    PrimitiveVariables domainCellPrim;
    //    PrimitiveVariables secondCellPrim;

    //    {
    //        ConservedVariables domainCellData;
    //        readCellData( domainCellIdx, dataBase, domainCellData );
    //        domainCellPrim = toPrimitiveVariables( domainCellData, parameters.K );
    //    }

    //    real rho = dataBase.data[RHO__(ghostCellIdx, dataBase.numberOfCells)];

    //    ghostCell_rhoS_1 = rho * ( two * boundaryCondition.S_1 - domainCellPrim.S_1 );
    //    ghostCell_rhoS_2 = rho * ( two * boundaryCondition.S_2 - domainCellPrim.S_2 );

    //    //ghostCell_rhoS_1 = rho * domainCellPrim.S_1;
    //    //ghostCell_rhoS_2 = rho * domainCellPrim.S_2;
    //}

    //{
    //    real x = dataBase.cellCenter[VEC_X(ghostCellIdx, dataBase.numberOfCells)];
    //    real y = dataBase.cellCenter[VEC_Y(ghostCellIdx, dataBase.numberOfCells)];

    //    real r = sqrt( x * x + y * y );

    //    if( r > 0.25 )
    //    {
    //        ghostCell_rhoS_1 *= four * (c1o2 - r);
    //        ghostCell_rhoS_2 *= four * (c1o2 - r);
    //    }
    //}

    //{
    //    dataBase.data[RHO_S_1(ghostCellIdx, dataBase.numberOfCells)] = ghostCell_rhoS_1;
    //    dataBase.data[RHO_S_2(ghostCellIdx, dataBase.numberOfCells)] = ghostCell_rhoS_2;
    //}

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    {
        PrimitiveVariables domainCellPrim;

        ConservedVariables domainCellData;
        readCellData( domainCellIdx, dataBase, domainCellData );
        domainCellPrim = toPrimitiveVariables( domainCellData, parameters.K );

        real dS_1 = ( boundaryCondition.S_1 * ( 1.0 - domainCellPrim.S_1 ) ) * parameters.dt;
        
        real x = dataBase.cellCenter[VEC_X(ghostCellIdx, dataBase.numberOfCells)];
        real y = dataBase.cellCenter[VEC_Y(ghostCellIdx, dataBase.numberOfCells)];

        real r = sqrt( x * x + y * y );

        //if( r > 0.25 ) dS_1 *= four * (c1o2 - r);

        domainCellPrim.S_1 += dS_1;

        domainCellPrim.S_2 = 1.0 - domainCellPrim.S_1;

        domainCellData = toConservedVariables( domainCellPrim, parameters.K );

        writeCellData(domainCellIdx, dataBase, domainCellData);
    }


#endif // USE_PASSIVE_SCALAR
}

PassiveScalarDiriclet::PassiveScalarDiriclet(SPtr<DataBase> dataBase, real S_1, real S_2)
    : BoundaryCondition( dataBase )
{
    this->S_1 = S_1;
    this->S_2 = S_2;
}

bool PassiveScalarDiriclet::isWall()
{
    return true;
}

bool PassiveScalarDiriclet::secondCellsNeeded()
{
    return false;
}

