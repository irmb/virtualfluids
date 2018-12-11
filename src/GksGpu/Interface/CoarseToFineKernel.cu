#include "Interface.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"

#include "DataBase/DataBaseStruct.h"

#include "FlowStateData/FlowStateData.cuh"

#include "Definitions/PassiveScalar.h"
#include "Definitions/MemoryAccessPattern.h"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void coarseToFineKernel  ( DataBaseStruct dataBase, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void coarseToFineFunction( DataBaseStruct dataBase, uint startIndex, uint index );

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
}

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

    ConservedVariables childCons [8];

    {
        real data [7];

        data[0] = dataBase.data[ RHO__(cellToCell[0], dataBase.numberOfCells) ];
        data[1] = dataBase.data[ RHO__(cellToCell[1], dataBase.numberOfCells) ];
        data[2] = dataBase.data[ RHO__(cellToCell[2], dataBase.numberOfCells) ];
        data[3] = dataBase.data[ RHO__(cellToCell[3], dataBase.numberOfCells) ];
        data[4] = dataBase.data[ RHO__(cellToCell[4], dataBase.numberOfCells) ];
        data[5] = dataBase.data[ RHO__(cellToCell[5], dataBase.numberOfCells) ];
        data[6] = dataBase.data[ RHO__(cellIndex    , dataBase.numberOfCells) ];

        //                                      PX        PY        PZ            MX        MY        MZ
        childCons[0].rho  = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
        childCons[1].rho  = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
        childCons[2].rho  = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
        childCons[3].rho  = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
        childCons[4].rho  = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
        childCons[5].rho  = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
        childCons[6].rho  = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
        childCons[7].rho  = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
    }

    {
        real data [7];

        data[0] = dataBase.data[ RHO_U(cellToCell[0], dataBase.numberOfCells) ];
        data[1] = dataBase.data[ RHO_U(cellToCell[1], dataBase.numberOfCells) ];
        data[2] = dataBase.data[ RHO_U(cellToCell[2], dataBase.numberOfCells) ];
        data[3] = dataBase.data[ RHO_U(cellToCell[3], dataBase.numberOfCells) ];
        data[4] = dataBase.data[ RHO_U(cellToCell[4], dataBase.numberOfCells) ];
        data[5] = dataBase.data[ RHO_U(cellToCell[5], dataBase.numberOfCells) ];
        data[6] = dataBase.data[ RHO_U(cellIndex    , dataBase.numberOfCells) ];

        //                                      PX        PY        PZ            MX        MY        MZ
        childCons[0].rhoU = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
        childCons[1].rhoU = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
        childCons[2].rhoU = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
        childCons[3].rhoU = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
        childCons[4].rhoU = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
        childCons[5].rhoU = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
        childCons[6].rhoU = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
        childCons[7].rhoU = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
    }

    {
        real data [7];

        data[0] = dataBase.data[ RHO_V(cellToCell[0], dataBase.numberOfCells) ];
        data[1] = dataBase.data[ RHO_V(cellToCell[1], dataBase.numberOfCells) ];
        data[2] = dataBase.data[ RHO_V(cellToCell[2], dataBase.numberOfCells) ];
        data[3] = dataBase.data[ RHO_V(cellToCell[3], dataBase.numberOfCells) ];
        data[4] = dataBase.data[ RHO_V(cellToCell[4], dataBase.numberOfCells) ];
        data[5] = dataBase.data[ RHO_V(cellToCell[5], dataBase.numberOfCells) ];
        data[6] = dataBase.data[ RHO_V(cellIndex    , dataBase.numberOfCells) ];

        //                                      PX        PY        PZ            MX        MY        MZ
        childCons[0].rhoV = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
        childCons[1].rhoV = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
        childCons[2].rhoV = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
        childCons[3].rhoV = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
        childCons[4].rhoV = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
        childCons[5].rhoV = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
        childCons[6].rhoV = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
        childCons[7].rhoV = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
    }

    {
        real data [7];

        data[0] = dataBase.data[ RHO_W(cellToCell[0], dataBase.numberOfCells) ];
        data[1] = dataBase.data[ RHO_W(cellToCell[1], dataBase.numberOfCells) ];
        data[2] = dataBase.data[ RHO_W(cellToCell[2], dataBase.numberOfCells) ];
        data[3] = dataBase.data[ RHO_W(cellToCell[3], dataBase.numberOfCells) ];
        data[4] = dataBase.data[ RHO_W(cellToCell[4], dataBase.numberOfCells) ];
        data[5] = dataBase.data[ RHO_W(cellToCell[5], dataBase.numberOfCells) ];
        data[6] = dataBase.data[ RHO_W(cellIndex    , dataBase.numberOfCells) ];

        //                                      PX        PY        PZ            MX        MY        MZ
        childCons[0].rhoW = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
        childCons[1].rhoW = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
        childCons[2].rhoW = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
        childCons[3].rhoW = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
        childCons[4].rhoW = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
        childCons[5].rhoW = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
        childCons[6].rhoW = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
        childCons[7].rhoW = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
    }

    {
        real data [7];

        data[0] = dataBase.data[ RHO_E(cellToCell[0], dataBase.numberOfCells) ];
        data[1] = dataBase.data[ RHO_E(cellToCell[1], dataBase.numberOfCells) ];
        data[2] = dataBase.data[ RHO_E(cellToCell[2], dataBase.numberOfCells) ];
        data[3] = dataBase.data[ RHO_E(cellToCell[3], dataBase.numberOfCells) ];
        data[4] = dataBase.data[ RHO_E(cellToCell[4], dataBase.numberOfCells) ];
        data[5] = dataBase.data[ RHO_E(cellToCell[5], dataBase.numberOfCells) ];
        data[6] = dataBase.data[ RHO_E(cellIndex    , dataBase.numberOfCells) ];

        //                                      PX        PY        PZ            MX        MY        MZ
        childCons[0].rhoE = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
        childCons[1].rhoE = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
        childCons[2].rhoE = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
        childCons[3].rhoE = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
        childCons[4].rhoE = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
        childCons[5].rhoE = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
        childCons[6].rhoE = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
        childCons[7].rhoE = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
    }

    #ifdef USE_PASSIVE_SCALAR
    {
        real data [7];

        data[0] = dataBase.data[ RHO_S(cellToCell[0], dataBase.numberOfCells) ];
        data[1] = dataBase.data[ RHO_S(cellToCell[1], dataBase.numberOfCells) ];
        data[2] = dataBase.data[ RHO_S(cellToCell[2], dataBase.numberOfCells) ];
        data[3] = dataBase.data[ RHO_S(cellToCell[3], dataBase.numberOfCells) ];
        data[4] = dataBase.data[ RHO_S(cellToCell[4], dataBase.numberOfCells) ];
        data[5] = dataBase.data[ RHO_S(cellToCell[5], dataBase.numberOfCells) ];
        data[6] = dataBase.data[ RHO_S(cellIndex    , dataBase.numberOfCells) ];

        //                                      PX        PY        PZ            MX        MY        MZ
        childCons[0].rhoS = data[6] + c1o8 * ( + data[0] + data[2] + data[4]     - data[1] - data[3] - data[5] ); // PX PY PZ
        childCons[1].rhoS = data[6] + c1o8 * ( + data[0] + data[2] - data[4]     - data[1] - data[3] + data[5] ); // PX PY MZ
        childCons[2].rhoS = data[6] + c1o8 * ( + data[0] - data[2] + data[4]     - data[1] + data[3] - data[5] ); // PX MY PZ
        childCons[3].rhoS = data[6] + c1o8 * ( + data[0] - data[2] - data[4]     - data[1] + data[3] + data[5] ); // PX MY MZ
        childCons[4].rhoS = data[6] + c1o8 * ( - data[0] + data[2] + data[4]     + data[1] - data[3] - data[5] ); // MX PY PZ
        childCons[5].rhoS = data[6] + c1o8 * ( - data[0] + data[2] - data[4]     + data[1] - data[3] + data[5] ); // MX PY MZ
        childCons[6].rhoS = data[6] + c1o8 * ( - data[0] - data[2] + data[4]     + data[1] + data[3] - data[5] ); // MX MY PZ
        childCons[7].rhoS = data[6] + c1o8 * ( - data[0] - data[2] - data[4]     + data[1] + data[3] + data[5] ); // MX MY MZ
    }
    #endif // USE_PASSIVE_SCALAR

#pragma unroll
    for( uint childIndex = 0; childIndex < 8; childIndex++ ){

        uint childCellIndex = dataBase.coarseToFine[ COARSE_TO_FINE( index, ( 1 + childIndex ), dataBase.numberOfFineGhostCells ) ];

        dataBase.data[ RHO__(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rho ;
        dataBase.data[ RHO_U(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoU;
        dataBase.data[ RHO_V(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoV;
        dataBase.data[ RHO_W(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoW;
        dataBase.data[ RHO_E(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoE;
    #ifdef USE_PASSIVE_SCALAR
	    dataBase.data[ RHO_S(childCellIndex, dataBase.numberOfCells) ] = childCons[childIndex].rhoS;
    #endif // USE_PASSIVE_SCALAR
    }
}
