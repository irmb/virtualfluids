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

    //////////////////////////////////////////////////////////////////////////

    real cellVolume = parameters.dx * parameters.dx * parameters.dx;

    ConservedVariables update;

    update.rho  = dataBase.dataUpdate[ RHO__(cellIndex, dataBase.numberOfCells) ] / cellVolume;
    update.rhoU = dataBase.dataUpdate[ RHO_U(cellIndex, dataBase.numberOfCells) ] / cellVolume;
    update.rhoV = dataBase.dataUpdate[ RHO_V(cellIndex, dataBase.numberOfCells) ] / cellVolume;
    update.rhoW = dataBase.dataUpdate[ RHO_W(cellIndex, dataBase.numberOfCells) ] / cellVolume;
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
    update.rhoE += force.x * dataBase.massFlux[ VEC_X(cellIndex, dataBase.numberOfCells) ] / ( six * parameters.dx * parameters.dx )
                 + force.y * dataBase.massFlux[ VEC_Y(cellIndex, dataBase.numberOfCells) ] / ( six * parameters.dx * parameters.dx ) 
                 + force.z * dataBase.massFlux[ VEC_Z(cellIndex, dataBase.numberOfCells) ] / ( six * parameters.dx * parameters.dx );

    dataBase.massFlux[ VEC_X(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.massFlux[ VEC_Y(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.massFlux[ VEC_Z(cellIndex, dataBase.numberOfCells) ] = zero;

    //////////////////////////////////////////////////////////////////////////

    dataBase.data[ RHO__(cellIndex, dataBase.numberOfCells) ] += update.rho ;
    dataBase.data[ RHO_U(cellIndex, dataBase.numberOfCells) ] += update.rhoU;
    dataBase.data[ RHO_V(cellIndex, dataBase.numberOfCells) ] += update.rhoV;
    dataBase.data[ RHO_W(cellIndex, dataBase.numberOfCells) ] += update.rhoW;
    dataBase.data[ RHO_E(cellIndex, dataBase.numberOfCells) ] += update.rhoE;

#ifdef USE_PASSIVE_SCALAR
	update.rhoS_1 = dataBase.dataUpdate[ RHO_S_1(cellIndex, dataBase.numberOfCells) ] / cellVolume;
	update.rhoS_2 = dataBase.dataUpdate[ RHO_S_2(cellIndex, dataBase.numberOfCells) ] / cellVolume;

    dataBase.dataUpdate[ RHO_S_1(cellIndex, dataBase.numberOfCells) ] = zero;
    dataBase.dataUpdate[ RHO_S_2(cellIndex, dataBase.numberOfCells) ] = zero;

    dataBase.data[ RHO_S_1(cellIndex, dataBase.numberOfCells) ] += update.rhoS_1;
    dataBase.data[ RHO_S_2(cellIndex, dataBase.numberOfCells) ] += update.rhoS_2;
#endif // USE_PASSIVE_SCALAR

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_PASSIVE_SCALAR
    if (true)
    {
        CellProperties cellProperties = dataBase.cellProperties[ cellIndex ];

        if( !isCellProperties( cellProperties, CELL_PROPERTIES_FINE_GHOST ) )
        {
            PrimitiveVariables updatedPrimitive;
            ConservedVariables updatedConserved;

            updatedConserved.rho    = dataBase.data[RHO__(cellIndex, dataBase.numberOfCells)];
            updatedConserved.rhoU   = dataBase.data[RHO_U(cellIndex, dataBase.numberOfCells)];
            updatedConserved.rhoV   = dataBase.data[RHO_V(cellIndex, dataBase.numberOfCells)];
            updatedConserved.rhoW   = dataBase.data[RHO_W(cellIndex, dataBase.numberOfCells)];
            updatedConserved.rhoE   = dataBase.data[RHO_E(cellIndex, dataBase.numberOfCells)];
            updatedConserved.rhoS_1 = dataBase.data[RHO_S_1(cellIndex, dataBase.numberOfCells)];
            updatedConserved.rhoS_2 = dataBase.data[RHO_S_2(cellIndex, dataBase.numberOfCells)];

            updatedPrimitive = toPrimitiveVariables(updatedConserved, parameters.K);

            //////////////////////////////////////////////////////////////////////////

            real Y_F = updatedPrimitive.S_1;
            real Y_P = updatedPrimitive.S_2;

            real Y_A = one - Y_F - Y_P;

            real M = one / (Y_A / M_A
                          + Y_F / M_F
                          + Y_P / M_P);

            real X_A = Y_A * M / M_A;
            real X_F = Y_F * M / M_F;
            real X_P = Y_P * M / M_P;

            ///////////////////////////////////////////////////////////////////////////////

            real X_O2 = real(0.21) * X_A;

            ///////////////////////////////////////////////////////////////////////////////

            {
                real dX_F = fminf(X_F, c1o2 * X_O2);

                //if( dX_F < zero ) dX_F = zero;

                // Limit the combustion
                real dX_F_max = parameters.dt / real(0.0001) * real(0.0001);
                //real dX_F_max = real(0.0001);
                dX_F = fminf(dX_F_max, dX_F);

                //const real heatOfReaction = real(802310.0); // kJ / kmol
                const real heatOfReaction = real(800.0); // kJ / kmol

                real dn_F = updatedConserved.rho * dX_F / M;

                real releasedHeat = dn_F * heatOfReaction;

                ///////////////////////////////////////////////////////////////////////////////

                //real X_F_new = X_F - dX_F;
                //real X_P_new = X_P + dX_F;

                real X_A_new = X_A - two * dX_F / real(0.21);
                real X_F_new = X_F - dX_F;

                real X_P_new = one - X_A_new - X_F_new;

                real Z1 = X_F_new * M_F / M;
                real Z2 = X_P_new * M_P / M;

                //////////////////////////////////////////////////////////////////////////

                //if( Z1 < zero ) { Z2 -= Z1; Z1 = zero; }
                //if( Z2 < zero ) { Z1 -= Z2; Z2 = zero; }

                //if( Z1 > one  ) { Z2 += Z1 - one; Z1 = one; }
                //if( Z2 > one  ) { Z1 += Z2 - one; Z2 = one; }

                ///////////////////////////////////////////////////////////////////////////////

                dataBase.data[RHO_S_1(cellIndex, dataBase.numberOfCells)] = Z1 * updatedConserved.rho;
                dataBase.data[RHO_S_2(cellIndex, dataBase.numberOfCells)] = Z2 * updatedConserved.rho;

                dataBase.data[RHO_E(cellIndex, dataBase.numberOfCells)] = updatedConserved.rhoE + releasedHeat;
            }
        }
    }

#endif // USE_PASSIVE_SCALAR
}
