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
//! \file CellUpdate.cu
//! \ingroup CellUpdate
//! \author Stephan Lenz
//=======================================================================================
#include "CellUpdate.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <math.h>

#include "PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//! \brief This is a CUDA Kernel that computes the cell index and calls \ref cellUpdateFunction for this index
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

        //////////////////////////////////////////////////////////////////////////

        cons = cons + (c1o1 / cellVolume) * update;
        
        //////////////////////////////////////////////////////////////////////////

        // check if simulation crashed
        if( isnan(cons.rho ) ||
            isnan(cons.rhoU) ||
            isnan(cons.rhoV) ||
            isnan(cons.rhoW) ||
            isnan(cons.rhoE) )
        {
            *dataBase.crashCellIndex = cellIndex;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(parameters.forcingSchemeIdx == 0)
    {
        // consistent source term treatment of Tian et al. (2007)
        cons.rhoU += parameters.force.x * parameters.dt * cons.rho;
        cons.rhoV += parameters.force.y * parameters.dt * cons.rho;
        cons.rhoW += parameters.force.z * parameters.dt * cons.rho;
        cons.rhoE += parameters.force.x * dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] / (c6o1 * parameters.dx * parameters.dx)
                   + parameters.force.y * dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] / (c6o1 * parameters.dx * parameters.dx)
                   + parameters.force.z * dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] / (c6o1 * parameters.dx * parameters.dx);

        dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] = c0o1;
        dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] = c0o1;
        dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] = c0o1;
    }

    if(parameters.forcingSchemeIdx == 1)
    {
        // forcing only on density variation
        cons.rhoU += parameters.force.x * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoV += parameters.force.y * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoW += parameters.force.z * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoE += parameters.force.x * dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] / (c6o1 * parameters.dx * parameters.dx)
                   + parameters.force.y * dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] / (c6o1 * parameters.dx * parameters.dx)
                   + parameters.force.z * dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] / (c6o1 * parameters.dx * parameters.dx);

        dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] = c0o1;
        dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] = c0o1;
        dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] = c0o1;
    }

    if(parameters.forcingSchemeIdx == 2)
    {
        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);
        real lambda = prim.lambda;

        // forcing only on density variation
        cons.rhoU += parameters.force.x * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoV += parameters.force.y * parameters.dt * ( cons.rho - parameters.rhoRef );
        cons.rhoW += parameters.force.z * parameters.dt * ( cons.rho - parameters.rhoRef );

        prim = toPrimitiveVariables(cons, parameters.K);
        prim.lambda = lambda;
        cons = toConservedVariables(prim, parameters.K);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    writeCellData(cellIndex, dataBase, cons);
}
