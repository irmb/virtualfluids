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

namespace GksGpu {

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

    //CellProperties cellProperties = dataBase.cellProperties[ cellIndex ];

    //if( isCellProperties( cellProperties, CELL_PROPERTIES_FINE_GHOST ) );

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
        // dirty fix to exclude viscous heating: Part 1
        //ConservedVariables testCons = cons;
        //testCons.rho  += update.rho / cellVolume;
        //testCons.rhoE += update.rhoE/ cellVolume;
        //PrimitiveVariables testPrim = toPrimitiveVariables(testCons, parameters.K);
        //////////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////////////////////
        //if( cellIndex == 415179 )
        //{
        //    //printf( "rho   = %14.4e  |  dRho   = %14.4e \n", cons.rho   , (one / cellVolume) * update.rho    );
        //    //printf( "rhoU  = %14.4e  |  dRhoU  = %14.4e \n", cons.rhoU  , (one / cellVolume) * update.rhoU   );
        //    //printf( "rhoV  = %14.4e  |  dRhoV  = %14.4e \n", cons.rhoV  , (one / cellVolume) * update.rhoV   );
        //    //printf( "rhoW  = %14.4e  |  dRhoW  = %14.4e \n", cons.rhoW  , (one / cellVolume) * update.rhoW   );
        //    printf( "rhoE  = %14.4e  |  dRhoE  = %14.4e \n", cons.rhoE  , (one / cellVolume) * update.rhoE   );
        //    //printf( "rhoS1 = %14.4e  |  dRhoS1 = %14.4e \n", cons.rhoS_1, (one / cellVolume) * update.rhoS_1 );
        //    //printf( "rhoS2 = %14.4e  |  dRhoS2 = %14.4e \n", cons.rhoS_2, (one / cellVolume) * update.rhoS_2 );
        //    printf( "=================================================================\n" );
        //}
        //////////////////////////////////////////////////////////////////////////

        cons = cons + (c1o1 / cellVolume) * update;
        
        //////////////////////////////////////////////////////////////////////////
        // dirty fix to exclude viscous heating: Part 2
        //PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);
        //prim.lambda = testPrim.lambda;
        //cons = toConservedVariables( prim, parameters.K );
        //////////////////////////////////////////////////////////////////////////

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
        //cons.rhoE += parameters.force.x * dataBase.massFlux[VEC_X(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx)
        //           + parameters.force.y * dataBase.massFlux[VEC_Y(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx)
        //           + parameters.force.z * dataBase.massFlux[VEC_Z(cellIndex, dataBase.numberOfCells)] / (six * parameters.dx * parameters.dx);

        prim = toPrimitiveVariables(cons, parameters.K);
        prim.lambda = lambda;
        cons = toConservedVariables(prim, parameters.K);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // mass conserving fix for out of bounds scalars

    //if( cons.rhoS_1 < zero ) { cons.rhoS_2 -= cons.rhoS_1; cons.rhoS_1 = zero; }
    //if( cons.rhoS_2 < zero ) { cons.rhoS_1 -= cons.rhoS_2; cons.rhoS_2 = zero; }

    //if( cons.rhoS_1 > cons.rho ) { cons.rhoS_2 += cons.rhoS_1 - cons.rho; cons.rhoS_1 = cons.rho; }
    //if( cons.rhoS_2 > cons.rho ) { cons.rhoS_1 += cons.rhoS_2 - cons.rho; cons.rhoS_2 = cons.rho; }

    //if( cons.rhoS_1 + cons.rhoS_2 > cons.rho )
    //{
    //    real faktor = (Z1 + Z2);

    //    Z1 /= faktor;
    //    Z2 /= faktor;
    //}
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    chemicalReaction(dataBase, parameters, cellIndex, cons);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Dirty fix that limits the velocity

    //PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

    //real velocity = sqrt( prim.U * prim.U + prim.V * prim.V + prim.W * prim.W );

    //if( velocity > five  )
    //{
    //    prim.U *= five / velocity;
    //    prim.V *= five / velocity;
    //    prim.W *= five / velocity;
    //}

    //cons = toConservedVariables(prim, parameters.K);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    writeCellData(cellIndex, dataBase, cons);
}

} // namespace GksGpu
