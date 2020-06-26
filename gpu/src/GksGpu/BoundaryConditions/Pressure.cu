#include "Pressure.h"

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

namespace GksGpu{

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const PressureStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const PressureStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Pressure::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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
                                        const PressureStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const PressureStruct& boundaryCondition, 
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
        ghostCellPrim.lambda = c2o1 * domainCellPrim.lambda - secondCellPrim.lambda;
    #ifdef USE_PASSIVE_SCALAR
        ghostCellPrim.S_1    = c2o1 * domainCellPrim.S_1    - secondCellPrim.S_1;
        ghostCellPrim.S_2    = c2o1 * domainCellPrim.S_2    - secondCellPrim.S_2;
    #endif // USE_PASSIVE_SCALAR


        real rho0 = ( c2o1 * boundaryCondition.p0 * c1o2 * ( domainCellPrim.lambda + ghostCellPrim.lambda ) );
        ghostCellPrim.rho = c2o1 * rho0 - domainCellPrim.rho;

        //real lambda0 = ( c1o2 * ( domainCellPrim.rho + ghostCellPrim.rho  ) * c1o2 / boundaryCondition.p0 );
        //ghostCellPrim.lambda = two * lambda0 - domainCellPrim.lambda;
    }

    {
        ConservedVariables ghostCons = toConservedVariables( ghostCellPrim, parameters.K );

        writeCellData( ghostCellIdx, dataBase, ghostCons );
    }
}

Pressure::Pressure(SPtr<DataBase> dataBase, real p0)
    : BoundaryCondition( dataBase )
{
    this->p0 = p0;
}

bool Pressure::isWall()
{
    return false;
}

bool Pressure::secondCellsNeeded()
{
    return true;
}

} // namespace GksGpu
