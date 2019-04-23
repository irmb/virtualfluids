#include "InflowComplete.h"

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
                                                           const InflowCompleteStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const InflowCompleteStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void InflowComplete::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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

    getLastCudaError("Inflow::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const InflowCompleteStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const InflowCompleteStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];

    PrimitiveVariables ghostCellPrim;
    {
        PrimitiveVariables domainCellPrim;

        {
            ConservedVariables domainCellData;
            readCellData( domainCellIdx, dataBase, domainCellData );
            domainCellPrim = toPrimitiveVariables( domainCellData, parameters.K );
        }

    //    ghostCellPrim.rho    = two * boundaryCondition.prim.rho    - domainCellPrim.rho;
    //    ghostCellPrim.U      = two * boundaryCondition.prim.U      - domainCellPrim.U;
    //    ghostCellPrim.V      = two * boundaryCondition.prim.V      - domainCellPrim.V;
    //    ghostCellPrim.W      = two * boundaryCondition.prim.W      - domainCellPrim.W;
    //    ghostCellPrim.lambda = /*two * boundaryCondition.prim.lambda -*/ domainCellPrim.lambda;
    //#ifdef USE_PASSIVE_SCALAR
    //    ghostCellPrim.S_1    = two * boundaryCondition.prim.S_1    - domainCellPrim.S_1;
    //    ghostCellPrim.S_2    = two * boundaryCondition.prim.S_2    - domainCellPrim.S_2;
    //#endif // USE_PASSIVE_SCALAR

        ghostCellPrim.rho    = boundaryCondition.prim.rho;
        ghostCellPrim.U      = two * boundaryCondition.prim.U - domainCellPrim.U;
        ghostCellPrim.V      = two * boundaryCondition.prim.V - domainCellPrim.V;
        ghostCellPrim.W      = boundaryCondition.prim.W;
        ghostCellPrim.lambda = boundaryCondition.prim.lambda;
    #ifdef USE_PASSIVE_SCALAR
        ghostCellPrim.S_1    = boundaryCondition.prim.S_1;
        ghostCellPrim.S_2    = boundaryCondition.prim.S_2;
    #endif // USE_PASSIVE_SCALAR

        real y = dataBase.cellCenter[ VEC_Y(ghostCellIdx, dataBase.numberOfCells) ];
        real x = dataBase.cellCenter[ VEC_X(ghostCellIdx, dataBase.numberOfCells) ];

        real r = sqrt( y*y + x*x );

        ghostCellPrim.W   *= (one - four*r*r);
        ghostCellPrim.S_1 *=       (one - four*r*r);
        ghostCellPrim.S_2 = one - ghostCellPrim.S_1;
    }

    {
        ConservedVariables ghostCons = toConservedVariables( ghostCellPrim, parameters.K );

        writeCellData( ghostCellIdx, dataBase, ghostCons );
    }
}

InflowComplete::InflowComplete(SPtr<DataBase> dataBase, PrimitiveVariables prim)
    : BoundaryCondition( dataBase )
{
    this->prim = prim;
}

bool InflowComplete::isWall()
{
    return false;
}

bool InflowComplete::secondCellsNeeded()
{
    return false;
}

