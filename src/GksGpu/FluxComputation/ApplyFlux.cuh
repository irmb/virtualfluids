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
    double* dataUpdate = dataBase.dataUpdate;

    //real* dataUpdate = dataBase.dataUpdatePX;

    //if(direction == 'x') dataUpdate = dataBase.dataUpdatePX;
    //if(direction == 'y') dataUpdate = dataBase.dataUpdatePY;
    //if(direction == 'z') dataUpdate = dataBase.dataUpdatePZ;

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
    atomicAdd( &( dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] ), - (double)flux.rho  );
    atomicAdd( &( dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] ), - (double)flux.rhoU );
    atomicAdd( &( dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] ), - (double)flux.rhoV );
    atomicAdd( &( dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] ), - (double)flux.rhoW );
    atomicAdd( &( dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] ), - (double)flux.rhoE );
#ifdef USE_PASSIVE_SCALAR
	atomicAdd( &( dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] ), - (double)flux.rhoS_1 );
	atomicAdd( &( dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] ), - (double)flux.rhoS_2 );
#endif // USE_PASSIVE_SCALAR
    
    if( direction == 'x' )
        atomicAdd( &( dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'y' )
        atomicAdd( &( dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'z' )
        atomicAdd( &( dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
#else
#pragma omp atomic
    dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] -= flux.rho ;
#pragma omp atomic
    dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoU;
#pragma omp atomic
    dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoV;
#pragma omp atomic
    dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoW;
#pragma omp atomic
    dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoE;
#ifdef USE_PASSIVE_SCALAR
#pragma omp atomic
	dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_1;
	dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_2;
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
    double* dataUpdate = dataBase.dataUpdate;

    //real* dataUpdate = dataBase.dataUpdateMX;

    //if(direction == 'x') dataUpdate = dataBase.dataUpdateMX;
    //if(direction == 'y') dataUpdate = dataBase.dataUpdateMY;
    //if(direction == 'z') dataUpdate = dataBase.dataUpdateMZ;

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
    atomicAdd( &( dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] ),   (double)flux.rho  );
    atomicAdd( &( dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] ),   (double)flux.rhoU );
    atomicAdd( &( dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] ),   (double)flux.rhoV );
    atomicAdd( &( dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] ),   (double)flux.rhoW );
    atomicAdd( &( dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] ),   (double)flux.rhoE );
#ifdef USE_PASSIVE_SCALAR
	atomicAdd( &( dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] ),   (double)flux.rhoS_1 );
	atomicAdd( &( dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] ),   (double)flux.rhoS_2 );
#endif // USE_PASSIVE_SCALAR
    
    if( direction == 'x' )
        atomicAdd( &( dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'y' )
        atomicAdd( &( dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
    if( direction == 'z' )
        atomicAdd( &( dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
#else
#pragma omp atomic
    dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] += flux.rho ;
#pragma omp atomic
    dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] += flux.rhoU;
#pragma omp atomic
    dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] += flux.rhoV;
#pragma omp atomic
    dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] += flux.rhoW;
#pragma omp atomic
    dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] += flux.rhoE;
#ifdef USE_PASSIVE_SCALAR
#pragma omp atomic
	dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_1;
	dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_2;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//__host__ __device__ inline void applyFluxToNegCell_X( const DataBaseStruct& dataBase,
//                                                    const uint& negCellIdx,
//                                                    const ConservedVariables& flux,
//                                                    const char direction,
//                                                    const real& dt)
//{
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//    atomicAdd( &( dataBase.dataUpdatePX[ RHO__(negCellIdx, dataBase.numberOfCells) ] ), - flux.rho  );
//    atomicAdd( &( dataBase.dataUpdatePX[ RHO_U(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoU );
//    atomicAdd( &( dataBase.dataUpdatePX[ RHO_V(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoV );
//    atomicAdd( &( dataBase.dataUpdatePX[ RHO_W(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoW );
//    atomicAdd( &( dataBase.dataUpdatePX[ RHO_E(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoE );
//#ifdef USE_PASSIVE_SCALAR
//	atomicAdd( &( dataBase.dataUpdatePX[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_1 );
//	atomicAdd( &( dataBase.dataUpdatePX[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_2 );
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//        atomicAdd( &( dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'y' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'z' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//#else
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] -= flux.rho ;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoU;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoV;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoW;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoE;
//#ifdef USE_PASSIVE_SCALAR
//#pragma omp atomic
//	dataBase.dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_1;
//	dataBase.dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_2;
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'y' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'z' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//#endif
//
//}
//
//__host__ __device__ inline void applyFluxToPosCell_X( const DataBaseStruct& dataBase,
//                                                    const uint& posCellIdx,
//                                                    const ConservedVariables& flux,
//                                                    const char& direction,
//                                                    const real& dt )
//{
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//    atomicAdd( &( dataBase.dataUpdateMX[ RHO__(posCellIdx, dataBase.numberOfCells) ] ),   flux.rho  );
//    atomicAdd( &( dataBase.dataUpdateMX[ RHO_U(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoU );
//    atomicAdd( &( dataBase.dataUpdateMX[ RHO_V(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoV );
//    atomicAdd( &( dataBase.dataUpdateMX[ RHO_W(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoW );
//    atomicAdd( &( dataBase.dataUpdateMX[ RHO_E(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoE );
//#ifdef USE_PASSIVE_SCALAR
//	atomicAdd( &( dataBase.dataUpdateMX[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_1 );
//	atomicAdd( &( dataBase.dataUpdateMX[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_2 );
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//        atomicAdd( &( dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'y' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'z' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//#else
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] += flux.rho ;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] += flux.rhoU;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] += flux.rhoV;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] += flux.rhoW;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] += flux.rhoE;
//#ifdef USE_PASSIVE_SCALAR
//#pragma omp atomic
//	dataBase.dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_1;
//	dataBase.dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_2;
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'y' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'z' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//#endif
//}
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//__host__ __device__ inline void applyFluxToNegCell_Y( const DataBaseStruct& dataBase,
//                                                    const uint& negCellIdx,
//                                                    const ConservedVariables& flux,
//                                                    const char direction,
//                                                    const real& dt)
//{
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//    atomicAdd( &( dataBase.dataUpdatePY[ RHO__(negCellIdx, dataBase.numberOfCells) ] ), - flux.rho  );
//    atomicAdd( &( dataBase.dataUpdatePY[ RHO_U(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoU );
//    atomicAdd( &( dataBase.dataUpdatePY[ RHO_V(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoV );
//    atomicAdd( &( dataBase.dataUpdatePY[ RHO_W(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoW );
//    atomicAdd( &( dataBase.dataUpdatePY[ RHO_E(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoE );
//#ifdef USE_PASSIVE_SCALAR
//	atomicAdd( &( dataBase.dataUpdatePY[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_1 );
//	atomicAdd( &( dataBase.dataUpdatePY[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_2 );
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//        atomicAdd( &( dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'y' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'z' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//#else
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] -= flux.rho ;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoU;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoV;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoW;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoE;
//#ifdef USE_PASSIVE_SCALAR
//#pragma omp atomic
//	dataBase.dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_1;
//	dataBase.dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_2;
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'y' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'z' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//#endif
//
//}
//
//__host__ __device__ inline void applyFluxToPosCell_Y( const DataBaseStruct& dataBase,
//                                                    const uint& posCellIdx,
//                                                    const ConservedVariables& flux,
//                                                    const char& direction,
//                                                    const real& dt )
//{
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//    atomicAdd( &( dataBase.dataUpdateMY[ RHO__(posCellIdx, dataBase.numberOfCells) ] ),   flux.rho  );
//    atomicAdd( &( dataBase.dataUpdateMY[ RHO_U(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoU );
//    atomicAdd( &( dataBase.dataUpdateMY[ RHO_V(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoV );
//    atomicAdd( &( dataBase.dataUpdateMY[ RHO_W(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoW );
//    atomicAdd( &( dataBase.dataUpdateMY[ RHO_E(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoE );
//#ifdef USE_PASSIVE_SCALAR
//	atomicAdd( &( dataBase.dataUpdateMY[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_1 );
//	atomicAdd( &( dataBase.dataUpdateMY[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_2 );
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//        atomicAdd( &( dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'y' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'z' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//#else
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] += flux.rho ;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] += flux.rhoU;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] += flux.rhoV;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] += flux.rhoW;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] += flux.rhoE;
//#ifdef USE_PASSIVE_SCALAR
//#pragma omp atomic
//	dataBase.dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_1;
//	dataBase.dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_2;
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'y' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'z' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//#endif
//}
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//__host__ __device__ inline void applyFluxToNegCell_Z( const DataBaseStruct& dataBase,
//                                                    const uint& negCellIdx,
//                                                    const ConservedVariables& flux,
//                                                    const char direction,
//                                                    const real& dt)
//{
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//    atomicAdd( &( dataBase.dataUpdatePZ[ RHO__(negCellIdx, dataBase.numberOfCells) ] ), - flux.rho  );
//    atomicAdd( &( dataBase.dataUpdatePZ[ RHO_U(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoU );
//    atomicAdd( &( dataBase.dataUpdatePZ[ RHO_V(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoV );
//    atomicAdd( &( dataBase.dataUpdatePZ[ RHO_W(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoW );
//    atomicAdd( &( dataBase.dataUpdatePZ[ RHO_E(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoE );
//#ifdef USE_PASSIVE_SCALAR
//	atomicAdd( &( dataBase.dataUpdatePZ[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_1 );
//	atomicAdd( &( dataBase.dataUpdatePZ[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] ), - flux.rhoS_2 );
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//        atomicAdd( &( dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'y' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'z' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//#else
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] -= flux.rho ;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoU;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoV;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoW;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoE;
//#ifdef USE_PASSIVE_SCALAR
//#pragma omp atomic
//	dataBase.dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_1;
//	dataBase.dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] -= flux.rhoS_2;
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'y' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'z' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
//#endif
//
//}
//
//__host__ __device__ inline void applyFluxToPosCell_Z( const DataBaseStruct& dataBase,
//                                                    const uint& posCellIdx,
//                                                    const ConservedVariables& flux,
//                                                    const char& direction,
//                                                    const real& dt )
//{
//
//#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
//    atomicAdd( &( dataBase.dataUpdateMZ[ RHO__(posCellIdx, dataBase.numberOfCells) ] ),   flux.rho  );
//    atomicAdd( &( dataBase.dataUpdateMZ[ RHO_U(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoU );
//    atomicAdd( &( dataBase.dataUpdateMZ[ RHO_V(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoV );
//    atomicAdd( &( dataBase.dataUpdateMZ[ RHO_W(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoW );
//    atomicAdd( &( dataBase.dataUpdateMZ[ RHO_E(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoE );
//#ifdef USE_PASSIVE_SCALAR
//	atomicAdd( &( dataBase.dataUpdateMZ[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_1 );
//	atomicAdd( &( dataBase.dataUpdateMZ[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] ),   flux.rhoS_2 );
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//        atomicAdd( &( dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'y' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//    if( direction == 'z' )
//        atomicAdd( &( dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] ), flux.rho );
//#else
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] += flux.rho ;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] += flux.rhoU;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] += flux.rhoV;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] += flux.rhoW;
//#pragma omp atomic
//    dataBase.dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] += flux.rhoE;
//#ifdef USE_PASSIVE_SCALAR
//#pragma omp atomic
//	dataBase.dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_1;
//	dataBase.dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] += flux.rhoS_2;
//#endif // USE_PASSIVE_SCALAR
//    
//    if( direction == 'x' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_X(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'y' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Y(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//    if( direction == 'z' )
//#pragma omp atomic
//        dataBase.massFlux[ VEC_Z(posCellIdx, dataBase.numberOfCells) ] += flux.rho;
//#endif
//}

#endif