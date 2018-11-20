#include "CellUpdate.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void cellUpdateKernel  ( DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void cellUpdateFunction( DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CellUpdate::updateCells( SPtr<DataBase> dataBase, Parameters parameters, uint level )
{
    CudaUtility::CudaGrid grid( dataBase->perLevelCount[ level ].numberOfBulkCells, 32 );

    runKernel( cellUpdateKernel,
               cellUpdateFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               parameters,
               dataBase->perLevelCount[ level ].startOfCells );

    getLastCudaError("CellUpdate::updateCells( SPtr<DataBase> dataBase, Parameters parameters, uint level )");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void cellUpdateKernel(DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index > numberOfEntities ) return;

    cellUpdateFunction( dataBase, parameters, startIndex, index );
}

__host__ __device__ inline void cellUpdateFunction(DataBaseStruct dataBase, Parameters parameters, uint startIndex, uint index)
{
    uint cellIndex = startIndex + index;

    //////////////////////////////////////////////////////////////////////////

    real cellVolume = parameters.dx * parameters.dx;

    ConservedVariables update;

    update.rho  = dataBase.dataUpdate[ RHO__(cellIndex, dataBase.numberOfCells) ] / cellVolume;
    update.rhoU = dataBase.dataUpdate[ RHO_U(cellIndex, dataBase.numberOfCells) ] / cellVolume;
    update.rhoV = dataBase.dataUpdate[ RHO_V(cellIndex, dataBase.numberOfCells) ] / cellVolume;
    update.rhoV = dataBase.dataUpdate[ RHO_W(cellIndex, dataBase.numberOfCells) ] / cellVolume;
    update.rhoE = dataBase.dataUpdate[ RHO_E(cellIndex, dataBase.numberOfCells) ] / cellVolume;

    dataBase.dataUpdate[ RHO__(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.dataUpdate[ RHO_U(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.dataUpdate[ RHO_V(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.dataUpdate[ RHO_W(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.dataUpdate[ RHO_E(cellIndex, dataBase.numberOfCells) ] = zero;

    //////////////////////////////////////////////////////////////////////////

    real rho = dataBase.data[ RHO__(cellIndex, dataBase.numberOfCells) ] + update.rho;

    Vec3 force = parameters.force;

    update.rhoU += force.x * parameters.dt * rho ;
    update.rhoV += force.y * parameters.dt * rho ;
    update.rhoW += force.z * parameters.dt * rho ;
    update.rhoE += force.x * dataBase.massFlux[ VEC_X(cellIndex, dataBase.numberOfCells) ] / ( four * parameters.dx )
                 + force.y * dataBase.massFlux[ VEC_Y(cellIndex, dataBase.numberOfCells) ] / ( four * parameters.dx ) 
                 + force.z * dataBase.massFlux[ VEC_Z(cellIndex, dataBase.numberOfCells) ] / ( four * parameters.dx ) ;

    dataBase.massFlux[ VEC_X(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.massFlux[ VEC_Y(cellIndex, dataBase.numberOfCells) ] = zero;

    //////////////////////////////////////////////////////////////////////////

    dataBase.data[ RHO__(cellIndex, dataBase.numberOfCells) ] += update.rho ;
    dataBase.data[ RHO_U(cellIndex, dataBase.numberOfCells) ] += update.rhoU;
    dataBase.data[ RHO_V(cellIndex, dataBase.numberOfCells) ] += update.rhoV;
    dataBase.data[ RHO_W(cellIndex, dataBase.numberOfCells) ] += update.rhoW;
    dataBase.data[ RHO_E(cellIndex, dataBase.numberOfCells) ] += update.rhoE;

#ifdef USE_PASSIVE_SCALAR
	update.rhoS = dataBase.dataUpdate[ RHO_S(cellIndex, dataBase.numberOfCells) ] / cellVolume;

    dataBase.dataUpdate[ RHO_S(cellIndex, dataBase.numberOfCells) ] = zero;

    dataBase.data[ RHO_S(cellIndex, dataBase.numberOfCells) ] += update.rhoS;
#endif // USE_PASSIVE_SCALAR
}
