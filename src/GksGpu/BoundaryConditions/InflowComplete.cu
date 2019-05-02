#include "InflowComplete.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "FluxComputation/Moments.cuh"
#include "FluxComputation/ApplyFlux.cuh"
#include "FluxComputation/Transformation.cuh"
#include "FluxComputation/AssembleFlux.cuh"
#include "FluxComputation/ExpansionCoefficients.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const InflowCompleteStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const InflowCompleteStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void InflowComplete::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
                                                const Parameters parameters, 
                                                const uint level)
{    
    CudaUtility::CudaGrid grid( this->numberOfCellsPerLevel[ level ], 32 );

    runKernel( boundaryConditionKernel,
               boundaryConditionFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               this->toStruct(),
               parameters,
               this->startOfCellsPerLevel[ level ] );

    cudaDeviceSynchronize();

    getLastCudaError("Inflow::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const InflowCompleteStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const InflowCompleteStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( true )
    {
        PrimitiveVariables ghostCellPrim;
        {
            PrimitiveVariables domainCellPrim;

            {
                ConservedVariables domainCellData;
                readCellData(domainCellIdx, dataBase, domainCellData);
                domainCellPrim = toPrimitiveVariables(domainCellData, parameters.K);
            }

            //    ghostCellPrim.rho    = two * boundaryCondition.prim.rho    - domainCellPrim.rho;
            //    ghostCellPrim.U      = two * boundaryCondition.prim.U      - domainCellPrim.U;
            //    ghostCellPrim.V      = two * boundaryCondition.prim.V      - domainCellPrim.V;
            //    ghostCellPrim.W      = two * boundaryCondition.prim.W      - domainCellPrim.W;
            ghostCellPrim.lambda = /*two * boundaryCondition.prim.lambda -*/ domainCellPrim.lambda;
            //#ifdef USE_PASSIVE_SCALAR
            //    ghostCellPrim.S_1    = two * boundaryCondition.prim.S_1    - domainCellPrim.S_1;
            //    ghostCellPrim.S_2    = two * boundaryCondition.prim.S_2    - domainCellPrim.S_2;
            //#endif // USE_PASSIVE_SCALAR

            ghostCellPrim.rho = boundaryCondition.prim.rho;
            ghostCellPrim.U = two * boundaryCondition.prim.U - domainCellPrim.U;
            ghostCellPrim.V = two * boundaryCondition.prim.V - domainCellPrim.V;
            //ghostCellPrim.W = two * boundaryCondition.prim.W - domainCellPrim.W;
            ghostCellPrim.W      = boundaryCondition.prim.W;
            //ghostCellPrim.lambda = boundaryCondition.prim.lambda;
#ifdef USE_PASSIVE_SCALAR
            ghostCellPrim.S_1 = boundaryCondition.prim.S_1;
            ghostCellPrim.S_2 = boundaryCondition.prim.S_2;
#endif // USE_PASSIVE_SCALAR

            real y = dataBase.cellCenter[VEC_Y(ghostCellIdx, dataBase.numberOfCells)];
            real x = dataBase.cellCenter[VEC_X(ghostCellIdx, dataBase.numberOfCells)];

            real r = sqrt(y*y + x*x);

            ghostCellPrim.W *= (one - four*r*r);
#ifdef USE_PASSIVE_SCALAR
            ghostCellPrim.S_1 *= (one - four*r*r);
            ghostCellPrim.S_2 = boundaryCondition.prim.S_2 - ghostCellPrim.S_1;
#endif // USE_PASSIVE_SCALAR
        }

        {
            ConservedVariables ghostCons = toConservedVariables(ghostCellPrim, parameters.K);

            writeCellData(ghostCellIdx, dataBase, ghostCons);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( false )
    {
        PrimitiveVariables domainCellPrim;
        {
            ConservedVariables domainCellData;
            readCellData(domainCellIdx, dataBase, domainCellData);
            domainCellPrim = toPrimitiveVariables(domainCellData, parameters.K);
        }

        real momentU [NUMBER_OF_MOMENTS];
        real momentV [NUMBER_OF_MOMENTS];
        real momentW [NUMBER_OF_MOMENTS];
        real momentXi[NUMBER_OF_MOMENTS];

        PrimitiveVariables facePrim = boundaryCondition.prim;

        //facePrim.lambda = domainCellPrim.lambda;

        transformGlobalToLocal( facePrim, 'z' );

        computeMoments(facePrim, parameters.K, momentU, momentV, momentW, momentXi);

        ConservedVariables flux;

        flux.rho  = momentU[0 + 1];
        //flux.rhoU = momentU[1 + 1];

        //flux.rhoE = c1o2 * ( momentU[2 + 1]
        //                   + momentU[0 + 1] * momentV [2]
        //                   + momentU[0 + 1] * momentW [2]
        //                   + momentU[0 + 1] * momentXi[2] );

        flux.rhoE = momentU[0 + 1] * c1o4 * ( parameters.K + five ) / boundaryCondition.prim.lambda;

        //////////////////////////////////////////////////////////////////////////

#ifdef USE_PASSIVE_SCALAR
        flux.rhoS_1 = flux.rho * boundaryCondition.prim.S_1;
        flux.rhoS_2 = flux.rho * boundaryCondition.prim.S_2;
#endif // USE_PASSIVE_SCALAR

        flux   = ( parameters.dt * parameters.dx * parameters.dx * facePrim.rho ) * flux;

        transformLocalToGlobal( flux, 'z' );

        applyFluxToPosCell(dataBase, domainCellIdx, flux, 'z', parameters.dt);

        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( false )
    {
        PrimitiveVariables domainCellPrim;
        {
            ConservedVariables domainCellData;
            readCellData(domainCellIdx, dataBase, domainCellData);
            domainCellPrim = toPrimitiveVariables(domainCellData, parameters.K);
        }    

        PrimitiveVariables facePrim = boundaryCondition.prim;

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
        
        {
            ConservedVariables gradN, gradT1, gradT2;

            transformGlobalToLocal( gradN , 'z' );
            transformGlobalToLocal( gradT1, 'z' );
            transformGlobalToLocal( gradT2, 'z' );

            transformGlobalToLocal( facePrim, 'z' );

            computeExpansionCoefficients(facePrim, gradN , parameters.K, ax);
            computeExpansionCoefficients(facePrim, gradT1, parameters.K, ay);
            computeExpansionCoefficients(facePrim, gradT2, parameters.K, az);
        }

        //////////////////////////////////////////////////////////////////////////

        {
            ConservedVariables flux;
            {
                real momentU [ NUMBER_OF_MOMENTS ]; 
                real momentV [ NUMBER_OF_MOMENTS ]; 
                real momentW [ NUMBER_OF_MOMENTS ]; 
                real momentXi[ NUMBER_OF_MOMENTS ];

                computeMoments( facePrim, parameters.K, momentU, momentV, momentW, momentXi );

                Vec3 force = parameters.force;

                transformGlobalToLocal(force, 'z');

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

                    transformLocalToGlobal( flux, 'z' );
                }
            }

            applyFluxToPosCell(dataBase, domainCellIdx, flux, 'z', parameters.dt);
            applyFluxToNegCell(dataBase, ghostCellIdx , flux, 'z', parameters.dt);
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

InflowComplete::InflowComplete(SPtr<DataBase> dataBase, PrimitiveVariables prim)
    : BoundaryCondition( dataBase )
{
    this->prim = prim;
}

bool InflowComplete::isWall()
{
    return false;
}

bool InflowComplete::isFluxBC()
{
    return false;
}

bool InflowComplete::secondCellsNeeded()
{
    return false;
}

