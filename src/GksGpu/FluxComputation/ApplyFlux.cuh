#ifndef ApplyFlux_CUH
#define ApplyFlux_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void applyFluxToNegCell( const DataBaseStruct& dataBase,
                                                    const uint& negCellIdx,
                                                    const ConservedVariables& flux,
                                                    const char direction,
                                                    const real& dt)
{

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
    atomicAdd( &( dataBase.dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] ), - flux.rho  );
    atomicAdd( &( dataBase.dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoU );
    atomicAdd( &( dataBase.dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoV );
    atomicAdd( &( dataBase.dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoW );
    atomicAdd( &( dataBase.dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoE );
#ifdef USE_PASSIVE_SCALAR
	atomicAdd( &( dataBase.dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_1 );
	atomicAdd( &( dataBase.dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_2 );
#endif // USE_PASSIVE_SCALAR
    
    if( direction == 'x' )
        atomicAdd( &( dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'y' )
        atomicAdd( &( dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'z' )
        atomicAdd( &( dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
#else
#pragma omp atomic
    dataBase.dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] -= flux.rho ;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoU;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoV;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoW;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoE;
#ifdef USE_PASSIVE_SCALAR
#pragma omp atomic
	dataBase.dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_1;
	dataBase.dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_2;
#endif // USE_PASSIVE_SCALAR
    
    if( direction == 'x' )
#pragma omp atomic
        dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
    if( direction == 'y' )
#pragma omp atomic
        dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
    if( direction == 'z' )
#pragma omp atomic
        dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
#endif

}

__host__ __device__ inline void applyFluxToPosCell( const DataBaseStruct& dataBase,
                                                    const uint& posCellIdx,
                                                    const ConservedVariables& flux,
                                                    const char& direction,
                                                    const real& dt )
{

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
    atomicAdd( &( dataBase.dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] ),   flux.rho  );
    atomicAdd( &( dataBase.dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoU );
    atomicAdd( &( dataBase.dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoV );
    atomicAdd( &( dataBase.dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoW );
    atomicAdd( &( dataBase.dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoE );
#ifdef USE_PASSIVE_SCALAR
	atomicAdd( &( dataBase.dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_1 );
	atomicAdd( &( dataBase.dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_2 );
#endif // USE_PASSIVE_SCALAR
    
    if( direction == 'x' )
        atomicAdd( &( dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'y' )
        atomicAdd( &( dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'z' )
        atomicAdd( &( dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
#else
#pragma omp atomic
    dataBase.dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] += flux.rho ;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] += flux.rhoU;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] += flux.rhoV;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] += flux.rhoW;
#pragma omp atomic
    dataBase.dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] += flux.rhoE;
#ifdef USE_PASSIVE_SCALAR
#pragma omp atomic
	dataBase.dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_1;
	dataBase.dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_2;
#endif // USE_PASSIVE_SCALAR
    
    if( direction == 'x' )
#pragma omp atomic
        dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
    if( direction == 'y' )
#pragma omp atomic
        dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
    if( direction == 'z' )
#pragma omp atomic
        dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
#endif
}

#endif