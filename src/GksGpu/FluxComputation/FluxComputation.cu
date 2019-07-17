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
#include "FluxComputation/Smagorinsky.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void fluxKernel  ( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void fluxFunction( DataBaseStruct dataBase, Parameters parameters, char direction, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FluxComputation::run( SPtr<DataBase> dataBase, Parameters parameters, uint level, bool evaluateCommFaces )
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
    //{
    //    CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFaces, 64);

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
    //////////////////////////////////////////////////////////////////////////
    if( evaluateCommFaces )
    {
        CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfFaces - dataBase->perLevelCount[level].numberOfInnerFaces, 64, CudaUtility::communicationStream);

        if( grid.numberOfEntities <= 0 ) return;

        runKernel(fluxKernel,
                  fluxFunction,
                  dataBase->getDeviceType(), grid,
                  dataBase->toStruct(),
                  parameters,
                  'x',
                  dataBase->perLevelCount[level].startOfFacesX + dataBase->perLevelCount[level].numberOfInnerFaces);
    }
    else
    {
        CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfInnerFaces, 64, CudaUtility::computeStream);

        runKernel(fluxKernel,
                  fluxFunction,
                  dataBase->getDeviceType(), grid,
                  dataBase->toStruct(),
                  parameters,
                  'x',
                  dataBase->perLevelCount[level].startOfFacesX);
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

    real K = parameters.K;

    direction = dataBase.faceOrientation[ faceIndex ];

    parameters.D1 = parameters.D;
    parameters.D2 = parameters.D;

    //////////////////////////////////////////////////////////////////////////

    if( parameters.useSpongeLayer )
    {
        if( parameters.spongeLayerIdx == 0 )
        {
            real x = dataBase.faceCenter[VEC_X(faceIndex, dataBase.numberOfFaces)];
            real z = dataBase.faceCenter[VEC_Z(faceIndex, dataBase.numberOfFaces)];

            real muNew = parameters.mu;

            real zStart = real(0.35);

            if (fabsf(z) > zStart)
            {
                muNew += (fabs(z) - zStart) * c10o1 * c10o1 * c10o1 * parameters.mu;
            }

            parameters.mu = muNew;
        }
        if( parameters.spongeLayerIdx == 1 )
        {
            real x = dataBase.faceCenter[VEC_X(faceIndex, dataBase.numberOfFaces)];
            real z = dataBase.faceCenter[VEC_Z(faceIndex, dataBase.numberOfFaces)];

            real muNew = parameters.mu;

            real zStart = real(3.5);

            if (fabsf(z) > zStart)
            {
                muNew += (fabs(z) - zStart) * c10o1 * c10o1 * c10o1 * parameters.mu;
            }

            parameters.mu = muNew;
        }
        if( parameters.spongeLayerIdx == 2 )
        {
            real y = dataBase.faceCenter[VEC_Y(faceIndex, dataBase.numberOfFaces)];

            real muNew = parameters.mu;

            real yStart = real(3.0);

            if (fabsf(y) > yStart)
            {
                muNew += (fabs(y) - yStart) * c10o1 * c10o1 * parameters.mu;
            }

            parameters.mu = muNew;
        }
    }

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
        ax[i] = c0o1; 
        ay[i] = c0o1; 
        az[i] = c0o1; 
        at[i] = c0o1;
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
                                     facePrim,
                                     K);

        transformGlobalToLocal( gradN , direction );
        transformGlobalToLocal( gradT1, direction );
        transformGlobalToLocal( gradT2, direction );

        transformGlobalToLocal( facePrim, direction );

        computeExpansionCoefficients(facePrim, gradN , K, ax);
        computeExpansionCoefficients(facePrim, gradT1, K, ay);
        computeExpansionCoefficients(facePrim, gradT2, K, az);

        //////////////////////////////////////////////////////////////////////////

        if(parameters.useSmagorinsky)
        {
            real muTurb = getTurbulentViscositySmagorinsky( parameters, facePrim, gradN, gradT1, gradT2 );

            if( muTurb > parameters.mu )
            {
                real turbSc = real(0.3);
                real turbPr = real(0.5);

                parameters.mu = muTurb;

                parameters.D  = muTurb / turbSc;
                parameters.Pr = turbPr;
            }
        }

        //////////////////////////////////////////////////////////////////////////

        if(parameters.useTemperatureLimiter){
            real k = parameters.mu / parameters.Pr;

            real dUdx1 = ( gradN.rhoU  - facePrim.U * gradN.rho  );
            real dUdx2 = ( gradT1.rhoU - facePrim.U * gradT1.rho );
            real dUdx3 = ( gradT2.rhoU - facePrim.U * gradT2.rho );
    
            real dVdx1 = ( gradN.rhoV  - facePrim.V * gradN.rho  );
            real dVdx2 = ( gradT1.rhoV - facePrim.V * gradT1.rho );
            real dVdx3 = ( gradT2.rhoV - facePrim.V * gradT2.rho );
    
            real dWdx1 = ( gradN.rhoW  - facePrim.W * gradN.rho  );
            real dWdx2 = ( gradT1.rhoW - facePrim.W * gradT1.rho );
            real dWdx3 = ( gradT2.rhoW - facePrim.W * gradT2.rho );
    
            real dEdx1 = ( gradN.rhoE  - facePrim.W * gradN.rho  );
            real dEdx2 = ( gradT1.rhoE - facePrim.W * gradT1.rho );
            real dEdx3 = ( gradT2.rhoE - facePrim.W * gradT2.rho );

            real dTdx1 = dEdx1 - c2o1 * facePrim.U * dUdx1 - c2o1 * facePrim.V * dVdx1 - c2o1 * facePrim.W * dWdx1;
            real dTdx2 = dEdx2 - c2o1 * facePrim.U * dUdx2 - c2o1 * facePrim.V * dVdx2 - c2o1 * facePrim.W * dWdx2;
            real dTdx3 = dEdx3 - c2o1 * facePrim.U * dUdx3 - c2o1 * facePrim.V * dVdx3 - c2o1 * facePrim.W * dWdx3;

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // this one works for some time
            //real S = parameters.dx * parameters.dx * ( fabsf(dTdx1) + fabsf(dTdx2) + fabsf(dTdx3) );
            //k += real(0.00002) / real(0.015625) * S;

            //real T = getT(facePrim);
            //if( T > 20 )
                //k += parameters.temperatureLimiter * S;

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            //real S = parameters.dx * ( fabsf(dTdx1) + fabsf(dTdx2) + fabsf(dTdx3) );
            //k += real(0.00001) * real(0.0025) * S * S;
            
            //real kMax = real(0.01) * c1o2 * parameters.dx * parameters.dx / parameters.dt;
            real kMax = real(0.01);

            real S = parameters.dx * parameters.dx * ( dTdx1 * dTdx1 + dTdx2 * dTdx2 + dTdx3 * dTdx3 );
            k += fminf(kMax, parameters.temperatureLimiter * S);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // this one works for some time
            //real S = ( fabsf(dTdx1) + fabsf(dTdx2) + fabsf(dTdx3) );
            //k += real(0.00002) / real(0.015625) * S;
            //k += real(1.28e-4) * parameters.dx * parameters.dx * parameters.dx * S * S;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            parameters.Pr = parameters.mu / k;
        }
    }

    //////////////////////////////////////////////////////////////////////////

    parameters.D1 = parameters.D;
    parameters.D2 = parameters.D;

    if(parameters.usePassiveScalarLimiter){
    #ifdef USE_PASSIVE_SCALAR

        if( facePrim.S_1 < c0o1 ) parameters.D1 += - parameters.passiveScalarLimiter *   facePrim.S_1;
        if( facePrim.S_1 > c1o1  ) parameters.D1 +=   parameters.passiveScalarLimiter * ( facePrim.S_1 - c1o1 );
        
        parameters.D2 = parameters.D1;

        if( facePrim.S_2 < c0o1 ) parameters.D2 += - real(0.1)*parameters.passiveScalarLimiter *   facePrim.S_2;
        if( facePrim.S_2 > c1o1  ) parameters.D2 +=   real(0.1)*parameters.passiveScalarLimiter * ( facePrim.S_2 - c1o1 );
    #endif // USE_PASSIVE_SCALAR
    }

    //////////////////////////////////////////////////////////////////////////

    //{
    //#ifdef USE_PASSIVE_SCALAR
    //    if( facePrim.S_1 < zero )
    //    {
    //        parameters.D += - real(0.1) * facePrim.S_1;
    //    }
    //    if( facePrim.S_1 > one )
    //    {
    //        parameters.D +=   real(0.1) * ( facePrim.S_1 - one );
    //    }

    //#endif // USE_PASSIVE_SCALAR
    //}

    //////////////////////////////////////////////////////////////////////////

    {
        ConservedVariables flux;

        {
            real momentU [ NUMBER_OF_MOMENTS ]; 
            real momentV [ NUMBER_OF_MOMENTS ]; 
            real momentW [ NUMBER_OF_MOMENTS ]; 
            real momentXi[ NUMBER_OF_MOMENTS ];

            computeMoments( facePrim, K, momentU, momentV, momentW, momentXi );

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

                computeExpansionCoefficients( facePrim, timeGrad, K, at );
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

        #if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
            atomicAdd( &( dataBase.diffusivity[ negCellIdx ] ), parameters.D * parameters.dx * parameters.dx * parameters.dt );
            atomicAdd( &( dataBase.diffusivity[ posCellIdx ] ), parameters.D * parameters.dx * parameters.dx * parameters.dt );
        #endif

            CellProperties negCellProperties = dataBase.cellProperties[ negCellIdx ];
            CellProperties posCellProperties = dataBase.cellProperties[ posCellIdx ];

            //if( isCellProperties( negCellProperties, CELL_PROPERTIES_IS_FLUX_BC ) || 
            //    isCellProperties( posCellProperties, CELL_PROPERTIES_IS_FLUX_BC ) )
            //    return;

            if( isCellProperties( negCellProperties, CELL_PROPERTIES_WALL ) || 
                isCellProperties( posCellProperties, CELL_PROPERTIES_WALL ) )
            {
                flux.rho    = c0o1;
            #ifdef USE_PASSIVE_SCALAR
                flux.rhoS_1 = c0o1;
                flux.rhoS_2 = c0o1;
            #endif //USE_PASSIVE_SCALAR
            }

            if( isCellProperties( negCellProperties, CELL_PROPERTIES_IS_INSULATED ) || 
                isCellProperties( posCellProperties, CELL_PROPERTIES_IS_INSULATED ) )
            {
                flux.rhoE   = c0o1;
            }

            uint negCellParentIdx = dataBase.parentCell[ negCellIdx ];
            uint posCellParentIdx = dataBase.parentCell[ posCellIdx ];

            //if( !( negCellParentIdx != INVALID_INDEX ) != !( posCellParentIdx != INVALID_INDEX ) ) // XOR
            if( ( negCellParentIdx == INVALID_INDEX ) != ( posCellParentIdx == INVALID_INDEX ) ) // XOR
            {
                if( !isCellProperties( negCellProperties, CELL_PROPERTIES_GHOST ) && 
                    !isCellProperties( posCellProperties, CELL_PROPERTIES_GHOST ) )
                {
                    if (negCellParentIdx != INVALID_INDEX)
                    {
                        applyFluxToNegCell(dataBase, negCellParentIdx, flux, direction, parameters);
                    }

                    if (posCellParentIdx != INVALID_INDEX)
                    {
                        applyFluxToPosCell(dataBase, posCellParentIdx, flux, direction, parameters);
                    }
                }
            }

            applyFluxToNegCell(dataBase, negCellIdx, flux, direction, parameters);
            applyFluxToPosCell(dataBase, posCellIdx, flux, direction, parameters);
        }
    }
}
