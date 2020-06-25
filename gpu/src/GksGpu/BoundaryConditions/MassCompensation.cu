#include "MassCompensation.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

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

#include "FluxComputation/Moments.cuh"
#include "FluxComputation/ApplyFlux.cuh"
#include "FluxComputation/Transformation.cuh"
#include "FluxComputation/AssembleFlux.cuh"
#include "FluxComputation/ExpansionCoefficients.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu{

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const MassCompensationStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const MassCompensationStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void MassCompensation::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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
                                        const MassCompensationStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const MassCompensationStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    PrimitiveVariables domainCellPrim;
    {
        ConservedVariables domainCellData;
        readCellData(domainCellIdx, dataBase, domainCellData);
        domainCellPrim = toPrimitiveVariables(domainCellData, parameters.K);
    }

    //////////////////////////////////////////////////////////////////////////

    real p0 = c1o2 * boundaryCondition.rho / boundaryCondition.lambda;
    real p1 = c1o2 * domainCellPrim.rho / domainCellPrim.lambda;

    //////////////////////////////////////////////////////////////////////////

    if( p1 > p0 )
    {
        ConservedVariables flux;

        flux.rho = c2o1 * p0 * domainCellPrim.lambda - domainCellPrim.rho;

        //flux.rhoE = ( parameters.K + three ) / ( four * boundaryCondition.lambda ) * flux.rho;
        flux.rhoE = (parameters.K + c3o1) / (c4o1 * domainCellPrim.lambda) * flux.rho;

        flux = (parameters.dt * parameters.dx * parameters.dx) * flux;

        applyFluxToPosCell(dataBase, domainCellIdx, flux, 'z', parameters);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

MassCompensation::MassCompensation(SPtr<DataBase> dataBase, real rho, real velocity, real lambda)
    : BoundaryCondition( dataBase )
{
    this->rho      = rho;
    this->velocity = velocity;
    this->lambda   = lambda;
}

bool MassCompensation::isWall()
{
    return false;
}

bool MassCompensation::isFluxBC()
{
    return false;
}

bool MassCompensation::secondCellsNeeded()
{
    return false;
}

} // namespace GksGpu

