//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file ApplyFlux.cuh
//! \ingroup FluxComputation
//! \author Stephan Lenz
//=======================================================================================
#ifndef ApplyFlux_CUH
#define ApplyFlux_CUH

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! \brief This function applies the negative flux to the \ref DataBase::dataUpdate variable of the cell with index negCellIdx
//! 
//! \param dataBase     \ref DataBaseStruct that holds the memory pointers
//! \param negCellIdx   index of the negative cell
//! \param flux         flux that goes from negative to positive cell
//! \param direction    char with 'x', 'y' or 'z', used for mass flux computation
//! \param parameters   \ref Parameters struct
__host__ __device__ inline void applyFluxToNegCell( const DataBaseStruct& dataBase,
                                                    const uint& negCellIdx,
                                                    const ConservedVariables& flux,
                                                    const char direction,
                                                    const Parameters& parameters)
{
    realAccumulator* dataUpdate = dataBase.dataUpdate;

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
    atomicAdd( &( dataUpdate[ RHO__(negCellIdx, dataBase.numberOfCells) ] ), - (realAccumulator)flux.rho  );
    atomicAdd( &( dataUpdate[ RHO_U(negCellIdx, dataBase.numberOfCells) ] ), - (realAccumulator)flux.rhoU );
    atomicAdd( &( dataUpdate[ RHO_V(negCellIdx, dataBase.numberOfCells) ] ), - (realAccumulator)flux.rhoV );
    atomicAdd( &( dataUpdate[ RHO_W(negCellIdx, dataBase.numberOfCells) ] ), - (realAccumulator)flux.rhoW );
    atomicAdd( &( dataUpdate[ RHO_E(negCellIdx, dataBase.numberOfCells) ] ), - (realAccumulator)flux.rhoE );
#ifdef USE_PASSIVE_SCALAR
	atomicAdd( &( dataUpdate[ RHO_S_1(negCellIdx, dataBase.numberOfCells) ] ), - (realAccumulator)flux.rhoS_1 );
	atomicAdd( &( dataUpdate[ RHO_S_2(negCellIdx, dataBase.numberOfCells) ] ), - (realAccumulator)flux.rhoS_2 );
#endif // USE_PASSIVE_SCALAR
    
    if( parameters.forcingSchemeIdx == 0 || parameters.forcingSchemeIdx == 1 )
    {
        if (direction == 'x')
            atomicAdd(&(dataBase.massFlux[VEC_X(negCellIdx, dataBase.numberOfCells)]), flux.rho);
        if (direction == 'y')
            atomicAdd(&(dataBase.massFlux[VEC_Y(negCellIdx, dataBase.numberOfCells)]), flux.rho);
        if (direction == 'z')
            atomicAdd(&(dataBase.massFlux[VEC_Z(negCellIdx, dataBase.numberOfCells)]), flux.rho);
    }
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

    if( parameters.forcingSchemeIdx == 0 || parameters.forcingSchemeIdx == 1 )
    {
        if( direction == 'x' )
    #pragma omp atomic
            dataBase.massFlux[ VEC_X(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
        if( direction == 'y' )
    #pragma omp atomic
            dataBase.massFlux[ VEC_Y(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
        if( direction == 'z' )
    #pragma omp atomic
            dataBase.massFlux[ VEC_Z(negCellIdx, dataBase.numberOfCells) ] += flux.rho;
    }
#endif

}

//! \brief This function applies the positive flux to the \ref DataBase::dataUpdate variable of the cell with index posCellIdx
//! 
//! \param dataBase     \ref DataBaseStruct that holds the memory pointers
//! \param posCellIdx   index of the positive cell
//! \param flux         flux that goes from negative to positive cell
//! \param direction    char with 'x', 'y' or 'z', used for mass flux computation
//! \param parameters   \ref Parameters struct
__host__ __device__ inline void applyFluxToPosCell( const DataBaseStruct& dataBase,
                                                    const uint& posCellIdx,
                                                    const ConservedVariables& flux,
                                                    const char& direction,
                                                    const Parameters& parameters)
{
    realAccumulator* dataUpdate = dataBase.dataUpdate;

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
    atomicAdd( &( dataUpdate[ RHO__(posCellIdx, dataBase.numberOfCells) ] ),   (realAccumulator)flux.rho  );
    atomicAdd( &( dataUpdate[ RHO_U(posCellIdx, dataBase.numberOfCells) ] ),   (realAccumulator)flux.rhoU );
    atomicAdd( &( dataUpdate[ RHO_V(posCellIdx, dataBase.numberOfCells) ] ),   (realAccumulator)flux.rhoV );
    atomicAdd( &( dataUpdate[ RHO_W(posCellIdx, dataBase.numberOfCells) ] ),   (realAccumulator)flux.rhoW );
    atomicAdd( &( dataUpdate[ RHO_E(posCellIdx, dataBase.numberOfCells) ] ),   (realAccumulator)flux.rhoE );
#ifdef USE_PASSIVE_SCALAR
	atomicAdd( &( dataUpdate[ RHO_S_1(posCellIdx, dataBase.numberOfCells) ] ),   (realAccumulator)flux.rhoS_1 );
	atomicAdd( &( dataUpdate[ RHO_S_2(posCellIdx, dataBase.numberOfCells) ] ),   (realAccumulator)flux.rhoS_2 );
#endif // USE_PASSIVE_SCALAR
    
    if( parameters.forcingSchemeIdx == 0 || parameters.forcingSchemeIdx == 1 )
    {
        if (direction == 'x')
            atomicAdd(&(dataBase.massFlux[VEC_X(posCellIdx, dataBase.numberOfCells)]), flux.rho);
        if (direction == 'y')
            atomicAdd(&(dataBase.massFlux[VEC_Y(posCellIdx, dataBase.numberOfCells)]), flux.rho);
        if (direction == 'z')
            atomicAdd(&(dataBase.massFlux[VEC_Z(posCellIdx, dataBase.numberOfCells)]), flux.rho);
    }
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
    
    if( parameters.forcingSchemeIdx == 0 || parameters.forcingSchemeIdx == 1 )
    {
        if (direction == 'x')
    #pragma omp atomic
            dataBase.massFlux[VEC_X(posCellIdx, dataBase.numberOfCells)] += flux.rho;
        if (direction == 'y')
    #pragma omp atomic
            dataBase.massFlux[VEC_Y(posCellIdx, dataBase.numberOfCells)] += flux.rho;
        if (direction == 'z')
    #pragma omp atomic
            dataBase.massFlux[VEC_Z(posCellIdx, dataBase.numberOfCells)] += flux.rho;
    }
#endif
}

#endif