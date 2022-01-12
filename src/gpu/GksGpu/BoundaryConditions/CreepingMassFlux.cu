#include "CreepingMassFlux.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

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

#include "FluxComputation/Moments.cuh"
#include "FluxComputation/ApplyFlux.cuh"
#include "FluxComputation/Transformation.cuh"
#include "FluxComputation/AssembleFlux.cuh"
#include "FluxComputation/ExpansionCoefficients.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu{

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const CreepingMassFluxStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const CreepingMassFluxStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void CreepingMassFlux::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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

    getLastCudaError("HeatFlux::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const CreepingMassFluxStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const CreepingMassFluxStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    // uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    PrimitiveVariables domainCellPrim;
    {
        ConservedVariables domainCellData;
        readCellData(domainCellIdx, dataBase, domainCellData);
        domainCellPrim = toPrimitiveVariables(domainCellData, parameters.K);
    }

    ConservedVariables flux;

    if( boundaryCondition.velocity > 0.0 )
    {
        flux.rho = boundaryCondition.velocity * boundaryCondition.rho;

        flux.rhoE = ( parameters.K + c3o1 ) / ( c4o1 * boundaryCondition.lambda ) * flux.rho;
        //flux.rhoE = (parameters.K + three) / (four * domainCellPrim.lambda) * flux.rho;

    #ifdef USE_PASSIVE_SCALAR
        flux.rhoS_1 = flux.rho;
    #endif // USE_PASSIVE_SCALAR
        flux = (parameters.dt * parameters.dx * parameters.dx) * flux;

        applyFluxToPosCell(dataBase, domainCellIdx, flux, 'z', parameters);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

CreepingMassFlux::CreepingMassFlux(SPtr<DataBase> dataBase, real rho, real velocity, real lambda)
    : BoundaryCondition( dataBase )
{
    this->rho      = rho;
    this->velocity = velocity;
    this->lambda   = lambda;
}

bool CreepingMassFlux::isWall()
{
    return true;
}

bool CreepingMassFlux::isFluxBC()
{
    return true;
}

bool CreepingMassFlux::secondCellsNeeded()
{
    return false;
}

} // namespace GksGpu

