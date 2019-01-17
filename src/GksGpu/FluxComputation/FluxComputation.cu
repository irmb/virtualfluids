#include "FluxComputation.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"

#include "DataBase/DataBaseStruct.h"

#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"

#include "CellProperties/CellProperties.cuh"

#include "FluxComputation/Moments.cuh"
#include "FluxComputation/Reconstruction.cuh"
#include "FluxComputation/Transformation.cuh"
#include "FluxComputation/ExpansionCoefficients.cuh"
#include "FluxComputation/AssembleFlux.cuh"
#include "FluxComputation/ApplyFlux.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void fluxKernel  ( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void fluxFunction( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, uint level )
{
    //{
    //    CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFacesX, 128);

    //    runKernel(fluxKernel,
    //              fluxFunction,
    //              dataBase->getDeviceType(), grid,
    //              dataBase->toStruct(),
    //              parameters,
    //              'x',
    //              dataBase->perLevelCount[level].startOfFacesX);

    //    cudaDeviceSynchronize();

    //    getLastCudaError("FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, 'x', uint level )");
    //}
    //{
    //    CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFacesY, 128);

    //    runKernel(fluxKernel,
    //              fluxFunction,
    //              dataBase->getDeviceType(), grid,
    //              dataBase->toStruct(),
    //              parameters,
    //              'y',
    //              dataBase->perLevelCount[level].startOfFacesY);

    //    cudaDeviceSynchronize();

    //    getLastCudaError("FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, 'y', uint level )");
    //}
    //{
    //    CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFacesZ, 128);

    //    runKernel(fluxKernel,
    //              fluxFunction,
    //              dataBase->getDeviceType(), grid,
    //              dataBase->toStruct(),
    //              parameters,
    //              'z',
    //              dataBase->perLevelCount[level].startOfFacesZ);

    //    cudaDeviceSynchronize();

    //    getLastCudaError("FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, 'z', uint level )");
    //}
    //////////////////////////////////////////////////////////////////////////
    {
        CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFaces, 64);

        runKernel(fluxKernel,
                  fluxFunction,
                  dataBase->getDeviceType(), grid,
                  dataBase->toStruct(),
                  parameters,
                  'x',
                  dataBase->perLevelCount[level].startOfFacesX);

        cudaDeviceSynchronize();

        getLastCudaError("FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, 'x', uint level )");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void fluxKernel(DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    fluxFunction( dataBase, parameters, direction, startIndex, index );
}

__host__ __device__ inline void fluxFunction(DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint index)
{
    uint faceIndex = startIndex + index;

    direction = dataBase.faceOrientation[ faceIndex ];

    //////////////////////////////////////////////////////////////////////////

    //{   // SpongeLayer
    //    real x = dataBase.faceCenter[ VEC_X(faceIndex, dataBase.numberOfFaces) ];
    //    real y = dataBase.faceCenter[ VEC_Y(faceIndex, dataBase.numberOfFaces) ];

    //    real muNew = parameters.mu;

    //    if( x > two )
    //    {
    //        muNew += ( x - two ) * ten * ten * parameters.mu;
    //    }

    //    if( y > c1o2 )
    //    {
    //        muNew += ( y - c1o2 ) * ten * ten * ten * parameters.mu;
    //    }

    //    parameters.Pr = muNew / parameters.mu;

    //    parameters.mu = muNew;


    //}

    //////////////////////////////////////////////////////////////////////////

    PrimitiveVariables facePrim;

    //////////////////////////////////////////////////////////////////////////

    real ax[LENGTH_CELL_DATA];
    real ay[LENGTH_CELL_DATA];
    real az[LENGTH_CELL_DATA];
    real at[LENGTH_CELL_DATA];

#pragma unroll
    for( uint i = 0; i < LENGTH_CELL_DATA; i++ )
    { 
        ax[i] = zero; 
        ay[i] = zero; 
        az[i] = zero; 
        at[i] = zero;
    }

    //////////////////////////////////////////////////////////////////////////

    {
        ConservedVariables gradN, gradT1, gradT2;

        reconstructFiniteDifferences(faceIndex,
                                     dataBase,
                                     parameters,
                                     direction,
                                     gradN,
                                     gradT1,
                                     gradT2,
                                     facePrim);

        transformGlobalToLocal( gradN , direction );
        transformGlobalToLocal( gradT1, direction );
        transformGlobalToLocal( gradT2, direction );

        transformGlobalToLocal( facePrim, direction );

        computeExpansionCoefficients(facePrim, gradN , parameters.K, ax);
        computeExpansionCoefficients(facePrim, gradT1, parameters.K, ay);
        computeExpansionCoefficients(facePrim, gradT2, parameters.K, az);
    }

    //////////////////////////////////////////////////////////////////////////

    {
        ConservedVariables flux;

        ConservedVariables faceCons;

        {
            real momentU [ NUMBER_OF_MOMENTS ]; 
            real momentV [ NUMBER_OF_MOMENTS ]; 
            real momentW [ NUMBER_OF_MOMENTS ]; 
            real momentXi[ NUMBER_OF_MOMENTS ];

            computeMoments( facePrim, parameters.K, momentU, momentV, momentW, momentXi );

            Vec3 force = parameters.force;

            transformGlobalToLocal(force, direction);

            {
                ConservedVariables timeGrad;
                computeTimeDerivative( facePrim, 
                                       momentU, 
                                       momentV, 
                                       momentW, 
                                       momentXi, 
                                       ax, ay, az,
                                       force,
                                       timeGrad );

                computeExpansionCoefficients( facePrim, timeGrad, parameters.K, at );
            }
            {
                real timeCoefficients[4];
                computeTimeCoefficients( facePrim, parameters, timeCoefficients );

                real heatFlux;
                assembleFlux( facePrim, 
                              momentU, momentV, momentW, momentXi,
                              ax, ay, az, at, 
                              timeCoefficients, 
                              parameters,
                              force,
                              flux,
                              heatFlux );

                transformLocalToGlobal( flux, direction );
            }
        }

        //////////////////////////////////////////////////////////////////////////

        {
            uint negCellIdx = dataBase.faceToCell[ NEG_CELL(faceIndex, dataBase.numberOfFaces) ];
            uint posCellIdx = dataBase.faceToCell[ POS_CELL(faceIndex, dataBase.numberOfFaces) ];

            CellProperties negCellProperties = dataBase.cellProperties[ negCellIdx ];
            CellProperties posCellProperties = dataBase.cellProperties[ posCellIdx ];

            if( isCellProperties( negCellProperties, CELL_PROPERTIES_WALL ) || 
                isCellProperties( posCellProperties, CELL_PROPERTIES_WALL ) )
            {
                flux.rho = zero;
            }

            uint negCellParentIdx = dataBase.parentCell[ negCellIdx ];
            uint posCellParentIdx = dataBase.parentCell[ posCellIdx ];

            if( !( negCellParentIdx != INVALID_INDEX ) != !( posCellParentIdx != INVALID_INDEX ) ) // XOR
            {
                if( !isCellProperties( negCellProperties, CELL_PROPERTIES_GHOST ) && 
                    !isCellProperties( posCellProperties, CELL_PROPERTIES_GHOST ) )
                {
                    if (negCellParentIdx != INVALID_INDEX)
                    {
                        applyFluxToNegCell(dataBase, negCellParentIdx, flux, direction, parameters.dt);
                    }

                    if (posCellParentIdx != INVALID_INDEX)
                    {
                        applyFluxToPosCell(dataBase, posCellParentIdx, flux, direction, parameters.dt);
                    }
                }
            }

            applyFluxToNegCell(dataBase, negCellIdx, flux, direction, parameters.dt);
            applyFluxToPosCell(dataBase, posCellIdx, flux, direction, parameters.dt);
        }
    }
}
