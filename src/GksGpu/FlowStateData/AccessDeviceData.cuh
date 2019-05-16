#ifndef AccessDeviceData_CUH
#define AccessDeviceData_CUH

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

#include "Core/DataTypes.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void readCellData(const uint cellIdx, const DataBaseStruct& dataBase, ConservedVariables& cellCons)
{
    cellCons.rho  = dataBase.data[ RHO__( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoU = dataBase.data[ RHO_U( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoV = dataBase.data[ RHO_V( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoW = dataBase.data[ RHO_W( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoE = dataBase.data[ RHO_E( cellIdx, dataBase.numberOfCells ) ];
#ifdef USE_PASSIVE_SCALAR
	cellCons.rhoS_1 = dataBase.data[ RHO_S_1( cellIdx, dataBase.numberOfCells ) ];
	cellCons.rhoS_2 = dataBase.data[ RHO_S_2( cellIdx, dataBase.numberOfCells ) ];
#endif // USE_PASSIVE_SCALAR
}

__host__ __device__ inline void writeCellData(const uint cellIdx, const DataBaseStruct& dataBase, ConservedVariables& cellCons)
{
    dataBase.data[ RHO__( cellIdx, dataBase.numberOfCells ) ] = cellCons.rho ;
    dataBase.data[ RHO_U( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoU;
    dataBase.data[ RHO_V( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoV;
    dataBase.data[ RHO_W( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoW;
    dataBase.data[ RHO_E( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoE;
#ifdef USE_PASSIVE_SCALAR
	dataBase.data[ RHO_S_1( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoS_1;
	dataBase.data[ RHO_S_2( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoS_2;
#endif // USE_PASSIVE_SCALAR
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void readCellDataUpdate(const uint cellIdx, const DataBaseStruct& dataBase, ConservedVariables& cellCons)
{
    cellCons.rho  = dataBase.dataUpdate[ RHO__( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoU = dataBase.dataUpdate[ RHO_U( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoV = dataBase.dataUpdate[ RHO_V( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoW = dataBase.dataUpdate[ RHO_W( cellIdx, dataBase.numberOfCells ) ];
    cellCons.rhoE = dataBase.dataUpdate[ RHO_E( cellIdx, dataBase.numberOfCells ) ];
#ifdef USE_PASSIVE_SCALAR
	cellCons.rhoS_1 = dataBase.dataUpdate[ RHO_S_1( cellIdx, dataBase.numberOfCells ) ];
	cellCons.rhoS_2 = dataBase.dataUpdate[ RHO_S_2( cellIdx, dataBase.numberOfCells ) ];
#endif // USE_PASSIVE_SCALAR
}

__host__ __device__ inline void writeCellDataUpdate(const uint cellIdx, const DataBaseStruct& dataBase, ConservedVariables& cellCons)
{
    dataBase.dataUpdate[ RHO__( cellIdx, dataBase.numberOfCells ) ] = cellCons.rho ;
    dataBase.dataUpdate[ RHO_U( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoU;
    dataBase.dataUpdate[ RHO_V( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoV;
    dataBase.dataUpdate[ RHO_W( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoW;
    dataBase.dataUpdate[ RHO_E( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoE;
#ifdef USE_PASSIVE_SCALAR
	dataBase.dataUpdate[ RHO_S_1( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoS_1;
	dataBase.dataUpdate[ RHO_S_2( cellIdx, dataBase.numberOfCells ) ] = cellCons.rhoS_2;
#endif // USE_PASSIVE_SCALAR
}

#endif