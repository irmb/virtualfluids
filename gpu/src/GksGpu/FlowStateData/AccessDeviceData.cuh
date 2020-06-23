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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
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
//! \file AccessDeviceData.cuh
//! \ingroup FlowStateData
//! \author Stephan Lenz
//=======================================================================================
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

//! \brief reads flow state from DataBase::data
//! \param[in]  cellIdx    index of the cell from which the data is read
//! \param[in]  dataBase   reference to a DataBaseStruct
//! \param[out] cellCons   reference to ConservedVariables object in which the data is stored
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

//! \brief writes flow state to DataBase::data
//! \param[in]  cellIdx    index of the cell to which the data is written
//! \param[in]  dataBase   reference to a DataBaseStruct
//! \param[out] cellCons   reference to ConservedVariables object with the data that is written
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

//! \brief reads flow state update from DataBase::dataUpdate
//! \param[in]  cellIdx    index of the cell from which the data is read
//! \param[in]  dataBase   reference to a DataBaseStruct
//! \param[out] cellCons   reference to ConservedVariables object in which the data is stored
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

//! \brief writes flow state update to DataBase::dataUpdate
//! \param[in]  cellIdx    index of the cell to which the data is written
//! \param[in]  dataBase   reference to a DataBaseStruct
//! \param[out] cellCons   reference to ConservedVariables object with the data that is written
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