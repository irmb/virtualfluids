#include "CellUpdate.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/ThermalDependencies.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CellUpdate/Reaction.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void cellUpdateKernel  ( DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void cellUpdateFunction( DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CellUpdate::run( SPtr<DataBase> dataBase, Parameters parameters, uint level )
{
    CudaUtility::CudaGrid grid( dataBase->perLevelCount[ level ].numberOfBulkCells, 32 );

    runKernel( cellUpdateKernel,
               cellUpdateFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               parameters,
               dataBase->perLevelCount[ level ].startOfCells );

    cudaDeviceSynchronize();

    getLastCudaError("CellUpdate::run( SPtr<DataBase> dataBase, Parameters parameters, uint level )");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void cellUpdateKernel(DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    cellUpdateFunction( dataBase, parameters, startIndex, index );
}

__host__ __device__ inline void cellUpdateFunction(DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint index)
{
    uint cellIndex = startIndex + index;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real cellVolume = parameters.dx * parameters.dx * parameters.dx;

    ConservedVariables cons;

    readCellData      (cellIndex, dataBase, cons);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    {
        ConservedVariables update, zeroCons;
        readCellDataUpdate(cellIndex, dataBase, update);
        writeCellDataUpdate(cellIndex, dataBase, zeroCons);

        cons = cons + (one / cellVolume) * update;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(true)
    {
        // consistent source term treatment of Tian et al. (2007)
        cons.rhoU += parameters.force.x * parameters.dt * cons.rho;
        cons.rhoV += parameters.force.y * parameters.dt * cons.rho;
        cons.rhoW += parameters.force.z * parameters.dt * cons.rho;
        cons.rhoE += parameters.force.x * dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx)
                   + parameters.force.y * dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx)
                   + parameters.force.z * dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx);

        dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] = zero;
        dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] = zero;
        dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] = zero;
    }

    if(false)
    {
        // forcing only on density variation
        cons.rhoU += parameters.force.x * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoV += parameters.force.y * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoW += parameters.force.z * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoE += parameters.force.x * dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx)
                   + parameters.force.y * dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx)
                   + parameters.force.z * dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx);

        dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] = zero;
        dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] = zero;
        dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] = zero;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    chemicalReaction(dataBase, parameters, cellIndex, cons);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    writeCellData(cellIndex, dataBase, cons);
}
