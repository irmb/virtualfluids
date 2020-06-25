#include "Inflow.h"

#define _USE_MATH_DEFINES
#include <math.h>

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

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu{

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const InflowStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const InflowStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Inflow::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
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
                                        const InflowStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const InflowStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];
    uint secondCellIdx = boundaryCondition.secondCells[ startIndex + index ];

    PrimitiveVariables ghostCellPrim;
    {
        PrimitiveVariables domainCellPrim;
        PrimitiveVariables secondCellPrim;

        {
            ConservedVariables domainCellData;
            readCellData( domainCellIdx, dataBase, domainCellData );
            domainCellPrim = toPrimitiveVariables( domainCellData, parameters.K );
        }

        real factor = c1o1;
        //if( fabs(boundaryCondition.a1) > real(1.0e-6) )
        {
            real y = dataBase.cellCenter[ VEC_Y(ghostCellIdx, dataBase.numberOfCells) ];
            real z = dataBase.cellCenter[ VEC_Z(ghostCellIdx, dataBase.numberOfCells) ];

            real r = sqrt( y*y + z*z );

            //factor =  ( boundaryCondition.a0 
            //          + boundaryCondition.a1*y 
            //          + boundaryCondition.a2*y*y  ) * ( four / boundaryCondition.a1 / boundaryCondition.a1 );

            factor = ( boundaryCondition.a0 
                      + boundaryCondition.a1*r 
                      + boundaryCondition.a2*r*r  );

            //factor = one;
        }

        //ghostCellPrim.rho    = two *          boundaryCondition.rho        - domainCellPrim.rho;
        ghostCellPrim.U      = c2o1 * factor * boundaryCondition.velocity.x - domainCellPrim.U;
        ghostCellPrim.V      = c2o1 * factor * boundaryCondition.velocity.y - domainCellPrim.V;
        ghostCellPrim.W      = c2o1 * factor * boundaryCondition.velocity.z - domainCellPrim.W;
        ghostCellPrim.lambda = c2o1 *          boundaryCondition.lambda     - domainCellPrim.lambda;
    #ifdef USE_PASSIVE_SCALAR
        ghostCellPrim.S_1    = c2o1 *          boundaryCondition.S_1        - domainCellPrim.S_1;
        ghostCellPrim.S_2    = c2o1 *          boundaryCondition.S_2        - domainCellPrim.S_2;
    #endif // USE_PASSIVE_SCALAR
        
        real p = c1o2 * domainCellPrim.rho / domainCellPrim.lambda;
        ghostCellPrim.rho = c2o1 * p * ghostCellPrim.lambda;
    }

    {
        ConservedVariables ghostCons = toConservedVariables( ghostCellPrim, parameters.K );

        writeCellData( ghostCellIdx, dataBase, ghostCons );
    }
}

Inflow::Inflow(SPtr<DataBase> dataBase, Vec3 velocity, real lambda, real rho, real a0, real a1, real a2, real S_1, real S_2)
    : BoundaryCondition( dataBase )
{
    this->velocity       = velocity;
    this->lambda         = lambda;
    this->rho            = rho;
    this->S_1            = S_1;
    this->S_2            = S_2;

    this->a0             = a0;
    this->a1             = a1;
    this->a2             = a2;
}

bool Inflow::isWall()
{
    return false;
}

bool Inflow::secondCellsNeeded()
{
    return false;
}

} // namespace GksGpu

