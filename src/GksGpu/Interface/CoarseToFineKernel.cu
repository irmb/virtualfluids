#include "Interface.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"

#include "DataBase/DataBaseStruct.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"

#include "FlowStateData/AccessDeviceData.cuh"

#include "Definitions/PassiveScalar.h"
#include "Definitions/MemoryAccessPattern.h"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void coarseToFineKernel  ( DataBaseStruct dataBase, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void coarseToFineFunction                      ( DataBaseStruct dataBase, uint startIndex, uint index );
__host__ __device__ inline void coarseToFineFunctionPrimitiveInterpolation( DataBaseStruct dataBase, uint startIndex, uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void Interface::runCoarseToFine( SPtr<DataBase> dataBase, uint level )
{
    CudaUtility::CudaGrid grid(dataBase->perLevelCount[level].numberOfCoarseToFine, 128);

    runKernel(coarseToFineKernel,
              coarseToFineFunction,
              dataBase->getDeviceType(), grid,
              dataBase->toStruct(),
              dataBase->perLevelCount[level].startOfCoarseToFine);

    cudaDeviceSynchronize();

    getLastCudaError("void Interface::runCoarseToFine( SPtr<DataBase> dataBase, uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void coarseToFineKernel( DataBaseStruct dataBase, uint startIndex, uint numberOfEntities )
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    coarseToFineFunction( dataBase, startIndex, index );
    //coarseToFineFunctionPrimitiveInterpolation( dataBase, startIndex, index );
}

//__host__ __device__ inline void coarseToFineFunction( DataBaseStruct dataBase, uint startIndex, uint index )
//{
//    index += startIndex;
//
//    uint cellIndex = dataBase.coarseToFine[ COARSE_TO_FINE( index, 0, dataBase.numberOfFineGhostCells ) ];
//
//    uint cellToCell [6];
//
//    cellToCell[0] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 0, dataBase.numberOfCells ) ];
//    cellToCell[1] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 1, dataBase.numberOfCells ) ];
//    cellToCell[2] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 2, dataBase.numberOfCells ) ];
//    cellToCell[3] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 3, dataBase.numberOfCells ) ];
//    cellToCell[4] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 4, dataBase.numberOfCells ) ];
//    cellToCell[5] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 5, dataBase.numberOfCells ) ];
//
//    ConservedVariables childCons [8];
//
//    {
//        real data [7];
//
//        data[0] = dataBase.data[ RHO__(cellToCell[0], dataBase.numberOfCells) ];
//        data[1] = dataBase.data[ RHO__(cellToCell[1], dataBase.numberOfCells) ];
//        data[2] = dataBase.data[ RHO__(cellToCell[2], dataBase.numberOfCells) ];
//        data[3] = dataBase.data[ RHO__(cellToCell[3], dataBase.numberOfCells) ];
//        data[4] = dataBase.data[ RHO__(cellToCell[4], dataBase.numberOfCells) ];
//        data[5] = dataBase.data[ RHO__(cellToCell[5], dataBase.numberOfCells) ];
//        data[6] = dataBase.data[ RHO__(cellIndex    , dataBase.numberOfCells) ];
//
//        //                                      PX        PY        PZ            MX        MY        MZ
//        childCons[0].rho  = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
//        childCons[1].rho  = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
//        childCons[2].rho  = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
//        childCons[3].rho  = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
//        childCons[4].rho  = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
//        childCons[5].rho  = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
//        childCons[6].rho  = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
//        childCons[7].rho  = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
//    }
//
//    {
//        real data [7];
//
//        data[0] = dataBase.data[ RHO_U(cellToCell[0], dataBase.numberOfCells) ];
//        data[1] = dataBase.data[ RHO_U(cellToCell[1], dataBase.numberOfCells) ];
//        data[2] = dataBase.data[ RHO_U(cellToCell[2], dataBase.numberOfCells) ];
//        data[3] = dataBase.data[ RHO_U(cellToCell[3], dataBase.numberOfCells) ];
//        data[4] = dataBase.data[ RHO_U(cellToCell[4], dataBase.numberOfCells) ];
//        data[5] = dataBase.data[ RHO_U(cellToCell[5], dataBase.numberOfCells) ];
//        data[6] = dataBase.data[ RHO_U(cellIndex    , dataBase.numberOfCells) ];
//
//        //                                      PX        PY        PZ            MX        MY        MZ
//        childCons[0].rhoU = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
//        childCons[1].rhoU = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
//        childCons[2].rhoU = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
//        childCons[3].rhoU = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
//        childCons[4].rhoU = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
//        childCons[5].rhoU = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
//        childCons[6].rhoU = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
//        childCons[7].rhoU = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
//    }
//
//    {
//        real data [7];
//
//        data[0] = dataBase.data[ RHO_V(cellToCell[0], dataBase.numberOfCells) ];
//        data[1] = dataBase.data[ RHO_V(cellToCell[1], dataBase.numberOfCells) ];
//        data[2] = dataBase.data[ RHO_V(cellToCell[2], dataBase.numberOfCells) ];
//        data[3] = dataBase.data[ RHO_V(cellToCell[3], dataBase.numberOfCells) ];
//        data[4] = dataBase.data[ RHO_V(cellToCell[4], dataBase.numberOfCells) ];
//        data[5] = dataBase.data[ RHO_V(cellToCell[5], dataBase.numberOfCells) ];
//        data[6] = dataBase.data[ RHO_V(cellIndex    , dataBase.numberOfCells) ];
//
//        //                                      PX        PY        PZ            MX        MY        MZ
//        childCons[0].rhoV = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
//        childCons[1].rhoV = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
//        childCons[2].rhoV = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
//        childCons[3].rhoV = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
//        childCons[4].rhoV = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
//        childCons[5].rhoV = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
//        childCons[6].rhoV = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
//        childCons[7].rhoV = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
//    }
//
//    {
//        real data [7];
//
//        data[0] = dataBase.data[ RHO_W(cellToCell[0], dataBase.numberOfCells) ];
//        data[1] = dataBase.data[ RHO_W(cellToCell[1], dataBase.numberOfCells) ];
//        data[2] = dataBase.data[ RHO_W(cellToCell[2], dataBase.numberOfCells) ];
//        data[3] = dataBase.data[ RHO_W(cellToCell[3], dataBase.numberOfCells) ];
//        data[4] = dataBase.data[ RHO_W(cellToCell[4], dataBase.numberOfCells) ];
//        data[5] = dataBase.data[ RHO_W(cellToCell[5], dataBase.numberOfCells) ];
//        data[6] = dataBase.data[ RHO_W(cellIndex    , dataBase.numberOfCells) ];
//
//        //                                      PX        PY        PZ            MX        MY        MZ
//        childCons[0].rhoW = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
//        childCons[1].rhoW = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
//        childCons[2].rhoW = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
//        childCons[3].rhoW = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
//        childCons[4].rhoW = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
//        childCons[5].rhoW = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
//        childCons[6].rhoW = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
//        childCons[7].rhoW = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
//    }
//
//    {
//        real data [7];
//
//        data[0] = dataBase.data[ RHO_E(cellToCell[0], dataBase.numberOfCells) ];
//        data[1] = dataBase.data[ RHO_E(cellToCell[1], dataBase.numberOfCells) ];
//        data[2] = dataBase.data[ RHO_E(cellToCell[2], dataBase.numberOfCells) ];
//        data[3] = dataBase.data[ RHO_E(cellToCell[3], dataBase.numberOfCells) ];
//        data[4] = dataBase.data[ RHO_E(cellToCell[4], dataBase.numberOfCells) ];
//        data[5] = dataBase.data[ RHO_E(cellToCell[5], dataBase.numberOfCells) ];
//        data[6] = dataBase.data[ RHO_E(cellIndex    , dataBase.numberOfCells) ];
//
//        //                                      PX        PY        PZ            MX        MY        MZ
//        childCons[0].rhoE = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
//        childCons[1].rhoE = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
//        childCons[2].rhoE = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
//        childCons[3].rhoE = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
//        childCons[4].rhoE = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
//        childCons[5].rhoE = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
//        childCons[6].rhoE = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
//        childCons[7].rhoE = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
//    }
//
//    #ifdef USE_PASSIVE_SCALAR
//    {
//        {
//            real data[7];
//
//            data[0] = dataBase.data[RHO_S_1(cellToCell[0], dataBase.numberOfCells)];
//            data[1] = dataBase.data[RHO_S_1(cellToCell[1], dataBase.numberOfCells)];
//            data[2] = dataBase.data[RHO_S_1(cellToCell[2], dataBase.numberOfCells)];
//            data[3] = dataBase.data[RHO_S_1(cellToCell[3], dataBase.numberOfCells)];
//            data[4] = dataBase.data[RHO_S_1(cellToCell[4], dataBase.numberOfCells)];
//            data[5] = dataBase.data[RHO_S_1(cellToCell[5], dataBase.numberOfCells)];
//            data[6] = dataBase.data[RHO_S_1(cellIndex, dataBase.numberOfCells)];
//
//            //                                      PX        PY        PZ            MX        MY        MZ
//            childCons[0].rhoS_1 = data[6] + c1o8 * (+data[0] + data[2] + data[4] - data[1] - data[3] - data[5]); // PX PY PZ
//            childCons[1].rhoS_1 = data[6] + c1o8 * (+data[0] + data[2] - data[4] - data[1] - data[3] + data[5]); // PX PY MZ
//            childCons[2].rhoS_1 = data[6] + c1o8 * (+data[0] - data[2] + data[4] - data[1] + data[3] - data[5]); // PX MY PZ
//            childCons[3].rhoS_1 = data[6] + c1o8 * (+data[0] - data[2] - data[4] - data[1] + data[3] + data[5]); // PX MY MZ
//            childCons[4].rhoS_1 = data[6] + c1o8 * (-data[0] + data[2] + data[4] + data[1] - data[3] - data[5]); // MX PY PZ
//            childCons[5].rhoS_1 = data[6] + c1o8 * (-data[0] + data[2] - data[4] + data[1] - data[3] + data[5]); // MX PY MZ
//            childCons[6].rhoS_1 = data[6] + c1o8 * (-data[0] - data[2] + data[4] + data[1] + data[3] - data[5]); // MX MY PZ
//            childCons[7].rhoS_1 = data[6] + c1o8 * (-data[0] - data[2] - data[4] + data[1] + data[3] + data[5]); // MX MY MZ
//        }
//
//        {
//            real data[7];
//
//            data[0] = dataBase.data[RHO_S_2(cellToCell[0], dataBase.numberOfCells)];
//            data[1] = dataBase.data[RHO_S_2(cellToCell[1], dataBase.numberOfCells)];
//            data[2] = dataBase.data[RHO_S_2(cellToCell[2], dataBase.numberOfCells)];
//            data[3] = dataBase.data[RHO_S_2(cellToCell[3], dataBase.numberOfCells)];
//            data[4] = dataBase.data[RHO_S_2(cellToCell[4], dataBase.numberOfCells)];
//            data[5] = dataBase.data[RHO_S_2(cellToCell[5], dataBase.numberOfCells)];
//            data[6] = dataBase.data[RHO_S_2(cellIndex, dataBase.numberOfCells)];
//
//            //                                      PX        PY        PZ            MX        MY        MZ
//            childCons[0].rhoS_2 = data[6] + c1o8 * (+data[0] + data[2] + data[4] - data[1] - data[3] - data[5]); // PX PY PZ
//            childCons[1].rhoS_2 = data[6] + c1o8 * (+data[0] + data[2] - data[4] - data[1] - data[3] + data[5]); // PX PY MZ
//            childCons[2].rhoS_2 = data[6] + c1o8 * (+data[0] - data[2] + data[4] - data[1] + data[3] - data[5]); // PX MY PZ
//            childCons[3].rhoS_2 = data[6] + c1o8 * (+data[0] - data[2] - data[4] - data[1] + data[3] + data[5]); // PX MY MZ
//            childCons[4].rhoS_2 = data[6] + c1o8 * (-data[0] + data[2] + data[4] + data[1] - data[3] - data[5]); // MX PY PZ
//            childCons[5].rhoS_2 = data[6] + c1o8 * (-data[0] + data[2] - data[4] + data[1] - data[3] + data[5]); // MX PY MZ
//            childCons[6].rhoS_2 = data[6] + c1o8 * (-data[0] - data[2] + data[4] + data[1] + data[3] - data[5]); // MX MY PZ
//            childCons[7].rhoS_2 = data[6] + c1o8 * (-data[0] - data[2] - data[4] + data[1] + data[3] + data[5]); // MX MY MZ
//        }
//    }
//    #endif // USE_PASSIVE_SCALAR
//
//#pragma unroll
//    for( uint childIndex = 0; childIndex < 8; childIndex++ ){
//
//        uint childCellIndex = dataBase.coarseToFine[ COARSE_TO_FINE( index, ( 1 + childIndex ), dataBase.numberOfFineGhostCells ) ];
//
//        dataBase.data[ RHO__(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rho ;
//        dataBase.data[ RHO_U(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoU;
//        dataBase.data[ RHO_V(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoV;
//        dataBase.data[ RHO_W(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoW;
//        dataBase.data[ RHO_E(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoE;
//    #ifdef USE_PASSIVE_SCALAR
//	    dataBase.data[ RHO_S_1(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoS_1;
//	    dataBase.data[ RHO_S_2(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoS_2;
//    #endif // USE_PASSIVE_SCALAR
//    }
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void coarseToFineFunction( DataBaseStruct dataBase, uint startIndex, uint index )
{
    index += startIndex;

    uint cellIndex = dataBase.coarseToFine[ COARSE_TO_FINE( index, 0, dataBase.numberOfFineGhostCells ) ];

    uint cellToCell [6];

    cellToCell[0] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 0, dataBase.numberOfCells ) ];
    cellToCell[1] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 1, dataBase.numberOfCells ) ];
    cellToCell[2] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 2, dataBase.numberOfCells ) ];
    cellToCell[3] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 3, dataBase.numberOfCells ) ];
    cellToCell[4] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 4, dataBase.numberOfCells ) ];
    cellToCell[5] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 5, dataBase.numberOfCells ) ];

    ConservedVariables cons[7];

    readCellData(cellToCell[0], dataBase, cons[0]);
    readCellData(cellToCell[1], dataBase, cons[1]);
    readCellData(cellToCell[2], dataBase, cons[2]);
    readCellData(cellToCell[3], dataBase, cons[3]);
    readCellData(cellToCell[4], dataBase, cons[4]);
    readCellData(cellToCell[5], dataBase, cons[5]);
    readCellData(cellIndex, dataBase, cons[6]);

    ConservedVariables childCons [8];
    ConservedVariables zeroCons;

    //                                                 PX           PY           PZ               MX           MY           MZ
    childCons[0]    = cons[6]    + c1o8 * ( zeroCons + cons[0]    + cons[2]    + cons[4]        - cons[1]    - cons[3]    - cons[5]    ); // PX PY PZ
    childCons[1]    = cons[6]    + c1o8 * ( zeroCons + cons[0]    + cons[2]    - cons[4]        - cons[1]    - cons[3]    + cons[5]    ); // PX PY MZ
    childCons[2]    = cons[6]    + c1o8 * ( zeroCons + cons[0]    - cons[2]    + cons[4]        - cons[1]    + cons[3]    - cons[5]    ); // PX MY PZ
    childCons[3]    = cons[6]    + c1o8 * ( zeroCons + cons[0]    - cons[2]    - cons[4]        - cons[1]    + cons[3]    + cons[5]    ); // PX MY MZ
    childCons[4]    = cons[6]    + c1o8 * ( zeroCons - cons[0]    + cons[2]    + cons[4]        + cons[1]    - cons[3]    - cons[5]    ); // MX PY PZ
    childCons[5]    = cons[6]    + c1o8 * ( zeroCons - cons[0]    + cons[2]    - cons[4]        + cons[1]    - cons[3]    + cons[5]    ); // MX PY MZ
    childCons[6]    = cons[6]    + c1o8 * ( zeroCons - cons[0]    - cons[2]    + cons[4]        + cons[1]    + cons[3]    - cons[5]    ); // MX MY PZ
    childCons[7]    = cons[6]    + c1o8 * ( zeroCons - cons[0]    - cons[2]    - cons[4]        + cons[1]    + cons[3]    + cons[5]    ); // MX MY MZ
    
#ifdef USE_PASSIVE_SCALAR
    ConservedVariables min(  1.0e99,  1.0e99,  1.0e99,  1.0e99,  1.0e99,  1.0e99,  1.0e99 );
    ConservedVariables max( -1.0e99, -1.0e99, -1.0e99, -1.0e99, -1.0e99, -1.0e99, -1.0e99 );
#else
    ConservedVariables min(  1.0e99,  1.0e99,  1.0e99,  1.0e99,  1.0e99 );
    ConservedVariables max( -1.0e99, -1.0e99, -1.0e99, -1.0e99, -1.0e99 );
#endif

    for( uint index = 0; index < 7; index++ )
    {
        if( cons[ index ].rho    < min.rho    ) min.rho    = cons[ index ].rho   ;
        if( cons[ index ].rhoU   < min.rhoU   ) min.rhoU   = cons[ index ].rhoU  ;
        if( cons[ index ].rhoV   < min.rhoV   ) min.rhoV   = cons[ index ].rhoV  ;
        if( cons[ index ].rhoW   < min.rhoW   ) min.rhoW   = cons[ index ].rhoW  ;
        if( cons[ index ].rhoE   < min.rhoE   ) min.rhoE   = cons[ index ].rhoE  ;
    #ifdef USE_PASSIVE_SCALAR
        if( cons[ index ].rhoS_1 < min.rhoS_1 ) min.rhoS_1 = cons[ index ].rhoS_1;
        if( cons[ index ].rhoS_2 < min.rhoS_2 ) min.rhoS_2 = cons[ index ].rhoS_2;
    #endif

        if( cons[ index ].rho    > max.rho    ) max.rho    = cons[ index ].rho   ;
        if( cons[ index ].rhoU   > max.rhoU   ) max.rhoU   = cons[ index ].rhoU  ;
        if( cons[ index ].rhoV   > max.rhoV   ) max.rhoV   = cons[ index ].rhoV  ;
        if( cons[ index ].rhoW   > max.rhoW   ) max.rhoW   = cons[ index ].rhoW  ;
        if( cons[ index ].rhoE   > max.rhoE   ) max.rhoE   = cons[ index ].rhoE  ;
    #ifdef USE_PASSIVE_SCALAR
        if( cons[ index ].rhoS_1 > max.rhoS_1 ) max.rhoS_1 = cons[ index ].rhoS_1;
        if( cons[ index ].rhoS_2 > max.rhoS_2 ) max.rhoS_2 = cons[ index ].rhoS_2;
    #endif
    }

#pragma unroll
    for( uint index = 0; index < 8; index++ )
    {
        if( childCons[ index ].rho    < min.rho    ) childCons[ index ].rho    = min.rho    ;
        if( childCons[ index ].rhoU   < min.rhoU   ) childCons[ index ].rhoU   = min.rhoU   ;
        if( childCons[ index ].rhoV   < min.rhoV   ) childCons[ index ].rhoV   = min.rhoV   ;
        if( childCons[ index ].rhoW   < min.rhoW   ) childCons[ index ].rhoW   = min.rhoW   ;
        if( childCons[ index ].rhoE   < min.rhoE   ) childCons[ index ].rhoE   = min.rhoE   ;
    #ifdef USE_PASSIVE_SCALAR
        if( childCons[ index ].rhoS_1 < min.rhoS_1 ) childCons[ index ].rhoS_1 = min.rhoS_1 ;
        if( childCons[ index ].rhoS_2 < min.rhoS_2 ) childCons[ index ].rhoS_2 = min.rhoS_2 ;
    #endif
        
        if( childCons[ index ].rho    > max.rho    ) childCons[ index ].rho    = max.rho    ;
        if( childCons[ index ].rhoU   > max.rhoU   ) childCons[ index ].rhoU   = max.rhoU   ;
        if( childCons[ index ].rhoV   > max.rhoV   ) childCons[ index ].rhoV   = max.rhoV   ;
        if( childCons[ index ].rhoW   > max.rhoW   ) childCons[ index ].rhoW   = max.rhoW   ;
        if( childCons[ index ].rhoE   > max.rhoE   ) childCons[ index ].rhoE   = max.rhoE   ;
    #ifdef USE_PASSIVE_SCALAR
        if( childCons[ index ].rhoS_1 > max.rhoS_1 ) childCons[ index ].rhoS_1 = max.rhoS_1 ;
        if( childCons[ index ].rhoS_2 > max.rhoS_2 ) childCons[ index ].rhoS_2 = max.rhoS_2 ;
    #endif
    }

#pragma unroll
    for( uint childIndex = 0; childIndex < 8; childIndex++ ){

        uint childCellIndex = dataBase.coarseToFine[ COARSE_TO_FINE( index, ( 1 + childIndex ), dataBase.numberOfFineGhostCells ) ];

        writeCellData(childCellIndex, dataBase, childCons[childIndex]);
    }
}

__host__ __device__ inline void coarseToFineFunctionPrimitiveInterpolation( DataBaseStruct dataBase, uint startIndex, uint index )
{
    index += startIndex;

    uint cellIndex = dataBase.coarseToFine[ COARSE_TO_FINE( index, 0, dataBase.numberOfFineGhostCells ) ];

    uint cellToCell [6];

    cellToCell[0] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 0, dataBase.numberOfCells ) ];
    cellToCell[1] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 1, dataBase.numberOfCells ) ];
    cellToCell[2] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 2, dataBase.numberOfCells ) ];
    cellToCell[3] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 3, dataBase.numberOfCells ) ];
    cellToCell[4] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 4, dataBase.numberOfCells ) ];
    cellToCell[5] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 5, dataBase.numberOfCells ) ];

    PrimitiveVariables prim [7];
    ConservedVariables cons[7];

    readCellData(cellToCell[0], dataBase, cons[0]);
    readCellData(cellToCell[1], dataBase, cons[1]);
    readCellData(cellToCell[2], dataBase, cons[2]);
    readCellData(cellToCell[3], dataBase, cons[3]);
    readCellData(cellToCell[4], dataBase, cons[4]);
    readCellData(cellToCell[5], dataBase, cons[5]);
    readCellData(cellIndex, dataBase, cons[6]);

    prim[0] = toPrimitiveVariables(cons[0], c2o1);
    prim[1] = toPrimitiveVariables(cons[1], c2o1);
    prim[2] = toPrimitiveVariables(cons[2], c2o1);
    prim[3] = toPrimitiveVariables(cons[3], c2o1);
    prim[4] = toPrimitiveVariables(cons[4], c2o1);
    prim[5] = toPrimitiveVariables(cons[5], c2o1);
    prim[6] = toPrimitiveVariables(cons[6], c2o1);

    PrimitiveVariables childPrim [8];
    PrimitiveVariables zeroPrim;

    //                                                     PX           PY           PZ               MX           MY           MZ
        childPrim[0]    = prim[6]    + c1o8 * ( zeroPrim + prim[0]    + prim[2]    + prim[4]        - prim[1]    - prim[3]    - prim[5]    ); // PX PY PZ
        childPrim[1]    = prim[6]    + c1o8 * ( zeroPrim + prim[0]    + prim[2]    - prim[4]        - prim[1]    - prim[3]    + prim[5]    ); // PX PY MZ
        childPrim[2]    = prim[6]    + c1o8 * ( zeroPrim + prim[0]    - prim[2]    + prim[4]        - prim[1]    + prim[3]    - prim[5]    ); // PX MY PZ
        childPrim[3]    = prim[6]    + c1o8 * ( zeroPrim + prim[0]    - prim[2]    - prim[4]        - prim[1]    + prim[3]    + prim[5]    ); // PX MY MZ
        childPrim[4]    = prim[6]    + c1o8 * ( zeroPrim - prim[0]    + prim[2]    + prim[4]        + prim[1]    - prim[3]    - prim[5]    ); // MX PY PZ
        childPrim[5]    = prim[6]    + c1o8 * ( zeroPrim - prim[0]    + prim[2]    - prim[4]        + prim[1]    - prim[3]    + prim[5]    ); // MX PY MZ
        childPrim[6]    = prim[6]    + c1o8 * ( zeroPrim - prim[0]    - prim[2]    + prim[4]        + prim[1]    + prim[3]    - prim[5]    ); // MX MY PZ
        childPrim[7]    = prim[6]    + c1o8 * ( zeroPrim - prim[0]    - prim[2]    - prim[4]        + prim[1]    + prim[3]    + prim[5]    ); // MX MY MZ

#pragma unroll
    for( uint childIndex = 0; childIndex < 8; childIndex++ ){

        uint childCellIndex = dataBase.coarseToFine[ COARSE_TO_FINE( index, ( 1 + childIndex ), dataBase.numberOfFineGhostCells ) ];

        ConservedVariables childCons = toConservedVariables(childPrim[childIndex], c2o1);

        writeCellData(childCellIndex, dataBase, childCons);
    }
}
