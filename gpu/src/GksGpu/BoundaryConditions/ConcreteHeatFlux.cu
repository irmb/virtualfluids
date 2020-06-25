#include "ConcreteHeatFlux.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#include <thrust/host_vector.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"
#include "Core/Logger/Logger.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"
#include "FlowStateData/ThermalDependencies.cuh"

#include "FluxComputation/Moments.cuh"
#include "FluxComputation/ApplyFlux.cuh"
#include "FluxComputation/Transformation.cuh"
#include "FluxComputation/AssembleFlux.cuh"
#include "FluxComputation/ExpansionCoefficients.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu{

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const ConcreteHeatFluxStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const ConcreteHeatFluxStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void ConcreteHeatFlux::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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
                                        const ConcreteHeatFluxStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const ConcreteHeatFluxStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
#ifdef USE_PASSIVE_SCALAR
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];
    uint secondCellIdx = boundaryCondition.secondCells[ startIndex + index ];

    real dx = boundaryCondition.L / real(boundaryCondition.numberOfPoints + 1);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    PrimitiveVariables domainCellPrim;
    {
        ConservedVariables domainCellData;
        readCellData(domainCellIdx, dataBase, domainCellData);
        domainCellPrim = toPrimitiveVariables(domainCellData, parameters.K);
    }
    PrimitiveVariables secondCellPrim;
    {
        ConservedVariables secondCellData;
        readCellData(secondCellIdx, dataBase, secondCellData);
        secondCellPrim = toPrimitiveVariables(secondCellData, parameters.K);
    }

    real TF = c3o2 * getT(domainCellPrim) - c1o2 * getT(secondCellPrim);
    //real TF = getT(domainCellPrim);

    for( uint i = 0; i < boundaryCondition.numberOfPoints; i++ )
    {
        uint finiteDifferenceIndex = ( startIndex + index ) * boundaryCondition.numberOfPoints + i;

        real T0 = boundaryCondition.temperatures[ finiteDifferenceIndex ];

        real Tn;
        if( i == 0 )
            Tn = TF;
        else
            Tn = boundaryCondition.temperatures[ finiteDifferenceIndex - 1 ];

        real Tp;
        if( i == boundaryCondition.numberOfPoints - 1 )
            Tp = boundaryCondition.ambientTemperature;
        else
            Tp = boundaryCondition.temperatures[ finiteDifferenceIndex + 1 ];

        real dTdxx = ( Tp + Tn - c2o1 * T0 ) / ( dx * dx );

        boundaryCondition.temperatures[ finiteDifferenceIndex ] += parameters.dt * boundaryCondition.temperatureConductivity * dTdxx;
    }

    ConservedVariables flux;

    {
        real T0 = TF;
        real T1 = boundaryCondition.temperatures[ ( startIndex + index ) * boundaryCondition.numberOfPoints     ];
        real T2 = boundaryCondition.temperatures[ ( startIndex + index ) * boundaryCondition.numberOfPoints + 1 ];


        real k = boundaryCondition.temperatureConductivity * boundaryCondition.density * boundaryCondition.specificHeatCapacity;

        flux.rhoE = - k * ( - c3o2 * T0 + c2o1 * T1 - c1o2 * T2 ) / dx;
    }

    flux = (parameters.dt * parameters.dx * parameters.dx) * flux;

    applyFluxToNegCell(dataBase, domainCellIdx, flux, 'a', parameters);

#endif
}

ConcreteHeatFlux::~ConcreteHeatFlux()
{
    checkCudaErrors( cudaFree( this->temperatures ) );
}

ConcreteHeatFlux::ConcreteHeatFlux(SPtr<DataBase> dataBase, uint numberOfPoints, real temperatureConductivity, real density, real specificHeatCapacity, real L, real ambientTemperature)
    : BoundaryCondition( dataBase )
{
    this->numberOfPoints = numberOfPoints;

    this->temperatureConductivity = temperatureConductivity;
    this->density                 = density;
    this->specificHeatCapacity    = specificHeatCapacity;

    this->L = L;
    this->ambientTemperature = ambientTemperature;

    this->temperatures = nullptr;
}

void ConcreteHeatFlux::init()
{
    checkCudaErrors( cudaMalloc( &this->temperatures, sizeof(real) * numberOfPoints * this->numberOfCells ) );

    // initialize values
    thrust::device_ptr<real> dev_ptr(this->temperatures);
    thrust::fill(dev_ptr, dev_ptr + numberOfPoints * this->numberOfCells, this->ambientTemperature);

    this->ghostCellsHost.resize(this->numberOfCells);
    this->domainCellsHost.resize(this->numberOfCells);
    this->secondCellsHost.resize(this->numberOfCells);

    this->temperaturesHost.resize(this->numberOfPoints * this->numberOfCells);
}

void ConcreteHeatFlux::download()
{
    checkCudaErrors( cudaMemcpy(this->ghostCellsHost.data() , this->ghostCells , sizeof(uint) * this->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy(this->domainCellsHost.data(), this->domainCells, sizeof(uint) * this->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy(this->secondCellsHost.data(), this->secondCells, sizeof(uint) * this->numberOfCells, cudaMemcpyDeviceToHost ) );

    checkCudaErrors( cudaMemcpy(this->temperaturesHost.data(), this->temperatures, sizeof(real) * this->numberOfCells * this->numberOfPoints, cudaMemcpyDeviceToHost ) );
}

bool ConcreteHeatFlux::isWall()
{
    return true;
}

bool ConcreteHeatFlux::isInsulated()
{
    return true;
}

bool ConcreteHeatFlux::isFluxBC()
{
    return true;
}

bool ConcreteHeatFlux::secondCellsNeeded()
{
    return true;
}

} // namespace GksGpu

