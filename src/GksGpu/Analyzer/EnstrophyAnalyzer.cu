#include "EnstrophyAnalyzer.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cmath>
#include <sstream>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

#include <iomanip>

#include "Core/Logger/Logger.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu {

__global__                 void enstrophyKernel  ( DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint nx, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void enstrophyFunction( DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint nx, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool EnstrophyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

    thrust::device_vector<real> enstrophy( this->dataBase->perLevelCount[ 0 ].numberOfBulkCells );

    CudaUtility::CudaGrid grid( dataBase->perLevelCount[ 0 ].numberOfBulkCells, 32 );

    uint nx;
    if     ( dataBase->perLevelCount[ 0 ].numberOfBulkCells ==  64* 64* 64 ) nx =  64;
    else if( dataBase->perLevelCount[ 0 ].numberOfBulkCells == 128*128*128 ) nx = 128;
    else if( dataBase->perLevelCount[ 0 ].numberOfBulkCells == 256*256*256 ) nx = 256;

    runKernel( enstrophyKernel,
               enstrophyFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               parameters,
               enstrophy.data().get(),
               nx,
               dataBase->perLevelCount[ 0 ].startOfCells );

    getLastCudaError("KineticEnergyAnalyzer::run(uint iter)");

    real EnstrophyTmp = thrust::reduce( enstrophy.begin(), enstrophy.end(), c0o1, thrust::plus<real>() )
                      / real(dataBase->perLevelCount[ 0 ].numberOfBulkCells);

    this->enstrophyTimeSeries.push_back( EnstrophyTmp );

    //*logging::out << logging::Logger::INFO_HIGH << "EKin = " << EKin << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void enstrophyKernel(DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint nx, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    enstrophyFunction( dataBase, parameters, enstrophy, nx, startIndex, index );
}

__host__ __device__ void enstrophyFunction(DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint nx, uint startIndex, uint index)
{
    uint cellIndex = startIndex + index;

    //////////////////////////////////////////////////////////////////////////

    uint xIndex = ( cellIndex % ( nx*nx ) ) % nx;
    uint yIndex = ( cellIndex % ( nx*nx ) ) / nx;
    uint zIndex = ( cellIndex / ( nx*nx ) );

    uint xP1 = (( xIndex + 1 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint xP2 = (( xIndex + 2 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint xP3 = (( xIndex + 3 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint xP4 = (( xIndex + 4 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint xM1 = (( xIndex - 1 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint xM2 = (( xIndex - 2 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint xM3 = (( xIndex - 3 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint xM4 = (( xIndex - 4 )%nx) + (( yIndex     )%nx)*nx + (( zIndex     )%nx)*nx*nx;

    uint yP1 = (( xIndex     )%nx) + (( yIndex + 1 )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint yP2 = (( xIndex     )%nx) + (( yIndex + 2 )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint yP3 = (( xIndex     )%nx) + (( yIndex + 3 )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint yP4 = (( xIndex     )%nx) + (( yIndex + 4 )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint yM1 = (( xIndex     )%nx) + (( yIndex - 1 )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint yM2 = (( xIndex     )%nx) + (( yIndex - 2 )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint yM3 = (( xIndex     )%nx) + (( yIndex - 3 )%nx)*nx + (( zIndex     )%nx)*nx*nx;
    uint yM4 = (( xIndex     )%nx) + (( yIndex - 4 )%nx)*nx + (( zIndex     )%nx)*nx*nx;

    uint zP1 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex + 1 )%nx)*nx*nx;
    uint zP2 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex + 2 )%nx)*nx*nx;
    uint zP3 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex + 3 )%nx)*nx*nx;
    uint zP4 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex + 4 )%nx)*nx*nx;
    uint zM1 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex - 1 )%nx)*nx*nx;
    uint zM2 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex - 2 )%nx)*nx*nx;
    uint zM3 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex - 3 )%nx)*nx*nx;
    uint zM4 = (( xIndex     )%nx) + (( yIndex     )%nx)*nx + (( zIndex - 4 )%nx)*nx*nx;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real rho_xP1 = dataBase.data[ RHO__( xP1, dataBase.numberOfCells ) ];
    real rho_xP2 = dataBase.data[ RHO__( xP2, dataBase.numberOfCells ) ];
    real rho_xP3 = dataBase.data[ RHO__( xP3, dataBase.numberOfCells ) ];
    real rho_xP4 = dataBase.data[ RHO__( xP4, dataBase.numberOfCells ) ];
    real rho_xM1 = dataBase.data[ RHO__( xM1, dataBase.numberOfCells ) ];
    real rho_xM2 = dataBase.data[ RHO__( xM2, dataBase.numberOfCells ) ];
    real rho_xM3 = dataBase.data[ RHO__( xM3, dataBase.numberOfCells ) ];
    real rho_xM4 = dataBase.data[ RHO__( xM4, dataBase.numberOfCells ) ];

    real rho_yP1 = dataBase.data[ RHO__( yP1, dataBase.numberOfCells ) ];
    real rho_yP2 = dataBase.data[ RHO__( yP2, dataBase.numberOfCells ) ];
    real rho_yP3 = dataBase.data[ RHO__( yP3, dataBase.numberOfCells ) ];
    real rho_yP4 = dataBase.data[ RHO__( yP4, dataBase.numberOfCells ) ];
    real rho_yM1 = dataBase.data[ RHO__( yM1, dataBase.numberOfCells ) ];
    real rho_yM2 = dataBase.data[ RHO__( yM2, dataBase.numberOfCells ) ];
    real rho_yM3 = dataBase.data[ RHO__( yM3, dataBase.numberOfCells ) ];
    real rho_yM4 = dataBase.data[ RHO__( yM4, dataBase.numberOfCells ) ];

    real rho_zP1 = dataBase.data[ RHO__( zP1, dataBase.numberOfCells ) ];
    real rho_zP2 = dataBase.data[ RHO__( zP2, dataBase.numberOfCells ) ];
    real rho_zP3 = dataBase.data[ RHO__( zP3, dataBase.numberOfCells ) ];
    real rho_zP4 = dataBase.data[ RHO__( zP4, dataBase.numberOfCells ) ];
    real rho_zM1 = dataBase.data[ RHO__( zM1, dataBase.numberOfCells ) ];
    real rho_zM2 = dataBase.data[ RHO__( zM2, dataBase.numberOfCells ) ];
    real rho_zM3 = dataBase.data[ RHO__( zM3, dataBase.numberOfCells ) ];
    real rho_zM4 = dataBase.data[ RHO__( zM4, dataBase.numberOfCells ) ];

    //////////////////////////////////////////////////////////////////////////

    real U_xP1   = dataBase.data[ RHO_U( xP1, dataBase.numberOfCells ) ] / rho_xP1;
    real U_xP2   = dataBase.data[ RHO_U( xP2, dataBase.numberOfCells ) ] / rho_xP2;
    real U_xP3   = dataBase.data[ RHO_U( xP3, dataBase.numberOfCells ) ] / rho_xP3;
    real U_xP4   = dataBase.data[ RHO_U( xP4, dataBase.numberOfCells ) ] / rho_xP4;
    real U_xM1   = dataBase.data[ RHO_U( xM1, dataBase.numberOfCells ) ] / rho_xM1;
    real U_xM2   = dataBase.data[ RHO_U( xM2, dataBase.numberOfCells ) ] / rho_xM2;
    real U_xM3   = dataBase.data[ RHO_U( xM3, dataBase.numberOfCells ) ] / rho_xM3;
    real U_xM4   = dataBase.data[ RHO_U( xM4, dataBase.numberOfCells ) ] / rho_xM4;

    real U_yP1   = dataBase.data[ RHO_U( yP1, dataBase.numberOfCells ) ] / rho_yP1;
    real U_yP2   = dataBase.data[ RHO_U( yP2, dataBase.numberOfCells ) ] / rho_yP2;
    real U_yP3   = dataBase.data[ RHO_U( yP3, dataBase.numberOfCells ) ] / rho_yP3;
    real U_yP4   = dataBase.data[ RHO_U( yP4, dataBase.numberOfCells ) ] / rho_yP4;
    real U_yM1   = dataBase.data[ RHO_U( yM1, dataBase.numberOfCells ) ] / rho_yM1;
    real U_yM2   = dataBase.data[ RHO_U( yM2, dataBase.numberOfCells ) ] / rho_yM2;
    real U_yM3   = dataBase.data[ RHO_U( yM3, dataBase.numberOfCells ) ] / rho_yM3;
    real U_yM4   = dataBase.data[ RHO_U( yM4, dataBase.numberOfCells ) ] / rho_yM4;

    real U_zP1   = dataBase.data[ RHO_U( zP1, dataBase.numberOfCells ) ] / rho_zP1;
    real U_zP2   = dataBase.data[ RHO_U( zP2, dataBase.numberOfCells ) ] / rho_zP2;
    real U_zP3   = dataBase.data[ RHO_U( zP3, dataBase.numberOfCells ) ] / rho_zP3;
    real U_zP4   = dataBase.data[ RHO_U( zP4, dataBase.numberOfCells ) ] / rho_zP4;
    real U_zM1   = dataBase.data[ RHO_U( zM1, dataBase.numberOfCells ) ] / rho_zM1;
    real U_zM2   = dataBase.data[ RHO_U( zM2, dataBase.numberOfCells ) ] / rho_zM2;
    real U_zM3   = dataBase.data[ RHO_U( zM3, dataBase.numberOfCells ) ] / rho_zM3;
    real U_zM4   = dataBase.data[ RHO_U( zM4, dataBase.numberOfCells ) ] / rho_zM4;

    //////////////////////////////////////////////////////////////////////////

    real V_xP1   = dataBase.data[ RHO_V( xP1, dataBase.numberOfCells ) ] / rho_xP1;
    real V_xP2   = dataBase.data[ RHO_V( xP2, dataBase.numberOfCells ) ] / rho_xP2;
    real V_xP3   = dataBase.data[ RHO_V( xP3, dataBase.numberOfCells ) ] / rho_xP3;
    real V_xP4   = dataBase.data[ RHO_V( xP4, dataBase.numberOfCells ) ] / rho_xP4;
    real V_xM1   = dataBase.data[ RHO_V( xM1, dataBase.numberOfCells ) ] / rho_xM1;
    real V_xM2   = dataBase.data[ RHO_V( xM2, dataBase.numberOfCells ) ] / rho_xM2;
    real V_xM3   = dataBase.data[ RHO_V( xM3, dataBase.numberOfCells ) ] / rho_xM3;
    real V_xM4   = dataBase.data[ RHO_V( xM4, dataBase.numberOfCells ) ] / rho_xM4;

    real V_yP1   = dataBase.data[ RHO_V( yP1, dataBase.numberOfCells ) ] / rho_yP1;
    real V_yP2   = dataBase.data[ RHO_V( yP2, dataBase.numberOfCells ) ] / rho_yP2;
    real V_yP3   = dataBase.data[ RHO_V( yP3, dataBase.numberOfCells ) ] / rho_yP3;
    real V_yP4   = dataBase.data[ RHO_V( yP4, dataBase.numberOfCells ) ] / rho_yP4;
    real V_yM1   = dataBase.data[ RHO_V( yM1, dataBase.numberOfCells ) ] / rho_yM1;
    real V_yM2   = dataBase.data[ RHO_V( yM2, dataBase.numberOfCells ) ] / rho_yM2;
    real V_yM3   = dataBase.data[ RHO_V( yM3, dataBase.numberOfCells ) ] / rho_yM3;
    real V_yM4   = dataBase.data[ RHO_V( yM4, dataBase.numberOfCells ) ] / rho_yM4;

    real V_zP1   = dataBase.data[ RHO_V( zP1, dataBase.numberOfCells ) ] / rho_zP1;
    real V_zP2   = dataBase.data[ RHO_V( zP2, dataBase.numberOfCells ) ] / rho_zP2;
    real V_zP3   = dataBase.data[ RHO_V( zP3, dataBase.numberOfCells ) ] / rho_zP3;
    real V_zP4   = dataBase.data[ RHO_V( zP4, dataBase.numberOfCells ) ] / rho_zP4;
    real V_zM1   = dataBase.data[ RHO_V( zM1, dataBase.numberOfCells ) ] / rho_zM1;
    real V_zM2   = dataBase.data[ RHO_V( zM2, dataBase.numberOfCells ) ] / rho_zM2;
    real V_zM3   = dataBase.data[ RHO_V( zM3, dataBase.numberOfCells ) ] / rho_zM3;
    real V_zM4   = dataBase.data[ RHO_V( zM4, dataBase.numberOfCells ) ] / rho_zM4;

    //////////////////////////////////////////////////////////////////////////

    real W_xP1   = dataBase.data[ RHO_W( xP1, dataBase.numberOfCells ) ] / rho_xP1;
    real W_xP2   = dataBase.data[ RHO_W( xP2, dataBase.numberOfCells ) ] / rho_xP2;
    real W_xP3   = dataBase.data[ RHO_W( xP3, dataBase.numberOfCells ) ] / rho_xP3;
    real W_xP4   = dataBase.data[ RHO_W( xP4, dataBase.numberOfCells ) ] / rho_xP4;
    real W_xM1   = dataBase.data[ RHO_W( xM1, dataBase.numberOfCells ) ] / rho_xM1;
    real W_xM2   = dataBase.data[ RHO_W( xM2, dataBase.numberOfCells ) ] / rho_xM2;
    real W_xM3   = dataBase.data[ RHO_W( xM3, dataBase.numberOfCells ) ] / rho_xM3;
    real W_xM4   = dataBase.data[ RHO_W( xM4, dataBase.numberOfCells ) ] / rho_xM4;

    real W_yP1   = dataBase.data[ RHO_W( yP1, dataBase.numberOfCells ) ] / rho_yP1;
    real W_yP2   = dataBase.data[ RHO_W( yP2, dataBase.numberOfCells ) ] / rho_yP2;
    real W_yP3   = dataBase.data[ RHO_W( yP3, dataBase.numberOfCells ) ] / rho_yP3;
    real W_yP4   = dataBase.data[ RHO_W( yP4, dataBase.numberOfCells ) ] / rho_yP4;
    real W_yM1   = dataBase.data[ RHO_W( yM1, dataBase.numberOfCells ) ] / rho_yM1;
    real W_yM2   = dataBase.data[ RHO_W( yM2, dataBase.numberOfCells ) ] / rho_yM2;
    real W_yM3   = dataBase.data[ RHO_W( yM3, dataBase.numberOfCells ) ] / rho_yM3;
    real W_yM4   = dataBase.data[ RHO_W( yM4, dataBase.numberOfCells ) ] / rho_yM4;

    real W_zP1   = dataBase.data[ RHO_W( zP1, dataBase.numberOfCells ) ] / rho_zP1;
    real W_zP2   = dataBase.data[ RHO_W( zP2, dataBase.numberOfCells ) ] / rho_zP2;
    real W_zP3   = dataBase.data[ RHO_W( zP3, dataBase.numberOfCells ) ] / rho_zP3;
    real W_zP4   = dataBase.data[ RHO_W( zP4, dataBase.numberOfCells ) ] / rho_zP4;
    real W_zM1   = dataBase.data[ RHO_W( zM1, dataBase.numberOfCells ) ] / rho_zM1;
    real W_zM2   = dataBase.data[ RHO_W( zM2, dataBase.numberOfCells ) ] / rho_zM2;
    real W_zM3   = dataBase.data[ RHO_W( zM3, dataBase.numberOfCells ) ] / rho_zM3;
    real W_zM4   = dataBase.data[ RHO_W( zM4, dataBase.numberOfCells ) ] / rho_zM4;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    real dVdx = ( (c28o1 * c8o1) * ( V_xP1 - V_xM1 ) - (c7o1 * c8o1) * ( V_xP2 - V_xM2 ) + (c8o1 * c4o1 * c1o3) * ( V_xP3 - V_xM3 ) - ( V_xP4 - V_xM4 ) ) / (c7o1 * c10o1 * c4o1 * parameters.dx);
    real dWdx = ( (c28o1 * c8o1) * ( W_xP1 - W_xM1 ) - (c7o1 * c8o1) * ( W_xP2 - W_xM2 ) + (c8o1 * c4o1 * c1o3) * ( W_xP3 - W_xM3 ) - ( W_xP4 - W_xM4 ) ) / (c7o1 * c10o1 * c4o1 * parameters.dx);
    real dUdy = ( (c28o1 * c8o1) * ( U_yP1 - U_yM1 ) - (c7o1 * c8o1) * ( U_yP2 - U_yM2 ) + (c8o1 * c4o1 * c1o3) * ( U_yP3 - U_yM3 ) - ( U_yP4 - U_yM4 ) ) / (c7o1 * c10o1 * c4o1 * parameters.dx);
    real dWdy = ( (c28o1 * c8o1) * ( W_yP1 - W_yM1 ) - (c7o1 * c8o1) * ( W_yP2 - W_yM2 ) + (c8o1 * c4o1 * c1o3) * ( W_yP3 - W_yM3 ) - ( W_yP4 - W_yM4 ) ) / (c7o1 * c10o1 * c4o1 * parameters.dx);
    real dUdz = ( (c28o1 * c8o1) * ( U_zP1 - U_zM1 ) - (c7o1 * c8o1) * ( U_zP2 - U_zM2 ) + (c8o1 * c4o1 * c1o3) * ( U_zP3 - U_zM3 ) - ( U_zP4 - U_zM4 ) ) / (c7o1 * c10o1 * c4o1 * parameters.dx);
    real dVdz = ( (c28o1 * c8o1) * ( V_zP1 - V_zM1 ) - (c7o1 * c8o1) * ( V_zP2 - V_zM2 ) + (c8o1 * c4o1 * c1o3) * ( V_zP3 - V_zM3 ) - ( V_zP4 - V_zM4 ) ) / (c7o1 * c10o1 * c4o1 * parameters.dx);

    real tmpX = dWdy - dVdz;
    real tmpY = dUdz - dWdx;
    real tmpZ = dVdx - dUdy;

    //////////////////////////////////////////////////////////////////////////

    real rho = dataBase.data[ RHO__( cellIndex, dataBase.numberOfCells ) ];

    enstrophy[ cellIndex ] = c1o2 * rho * ( tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ );
}

//__host__ __device__ void enstrophyFunction(DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint startIndex, uint index)
//{
//    uint cellIndex = startIndex + index;
//
//    //////////////////////////////////////////////////////////////////////////
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
//    real rho [7];
//    real U   [6];
//    real V   [6];
//    real W   [6];
//
//    rho[0] = dataBase.data[ RHO__( cellToCell[0], dataBase.numberOfCells ) ];
//    rho[1] = dataBase.data[ RHO__( cellToCell[1], dataBase.numberOfCells ) ];
//    rho[2] = dataBase.data[ RHO__( cellToCell[2], dataBase.numberOfCells ) ];
//    rho[3] = dataBase.data[ RHO__( cellToCell[3], dataBase.numberOfCells ) ];
//    rho[4] = dataBase.data[ RHO__( cellToCell[4], dataBase.numberOfCells ) ];
//    rho[5] = dataBase.data[ RHO__( cellToCell[5], dataBase.numberOfCells ) ];
//    rho[6] = dataBase.data[ RHO__( cellIndex    , dataBase.numberOfCells ) ];
//
//    U  [0] = dataBase.data[ RHO_U( cellToCell[0], dataBase.numberOfCells ) ] / rho[0];
//    U  [1] = dataBase.data[ RHO_U( cellToCell[1], dataBase.numberOfCells ) ] / rho[1];
//    U  [2] = dataBase.data[ RHO_U( cellToCell[2], dataBase.numberOfCells ) ] / rho[2];
//    U  [3] = dataBase.data[ RHO_U( cellToCell[3], dataBase.numberOfCells ) ] / rho[3];
//    U  [4] = dataBase.data[ RHO_U( cellToCell[4], dataBase.numberOfCells ) ] / rho[4];
//    U  [5] = dataBase.data[ RHO_U( cellToCell[5], dataBase.numberOfCells ) ] / rho[5];
//
//    V  [0] = dataBase.data[ RHO_V( cellToCell[0], dataBase.numberOfCells ) ] / rho[0];
//    V  [1] = dataBase.data[ RHO_V( cellToCell[1], dataBase.numberOfCells ) ] / rho[1];
//    V  [2] = dataBase.data[ RHO_V( cellToCell[2], dataBase.numberOfCells ) ] / rho[2];
//    V  [3] = dataBase.data[ RHO_V( cellToCell[3], dataBase.numberOfCells ) ] / rho[3];
//    V  [4] = dataBase.data[ RHO_V( cellToCell[4], dataBase.numberOfCells ) ] / rho[4];
//    V  [5] = dataBase.data[ RHO_V( cellToCell[5], dataBase.numberOfCells ) ] / rho[5];
//
//    W  [0] = dataBase.data[ RHO_W( cellToCell[0], dataBase.numberOfCells ) ] / rho[0];
//    W  [1] = dataBase.data[ RHO_W( cellToCell[1], dataBase.numberOfCells ) ] / rho[1];
//    W  [2] = dataBase.data[ RHO_W( cellToCell[2], dataBase.numberOfCells ) ] / rho[2];
//    W  [3] = dataBase.data[ RHO_W( cellToCell[3], dataBase.numberOfCells ) ] / rho[3];
//    W  [4] = dataBase.data[ RHO_W( cellToCell[4], dataBase.numberOfCells ) ] / rho[4];
//    W  [5] = dataBase.data[ RHO_W( cellToCell[5], dataBase.numberOfCells ) ] / rho[5];
//
//    real dVdx = c1o2 * ( V[1] - V[0] ) / parameters.dx;
//    real dWdx = c1o2 * ( W[1] - W[0] ) / parameters.dx;
//
//    real dUdy = c1o2 * ( U[3] - U[2] ) / parameters.dx;
//    real dWdy = c1o2 * ( W[3] - W[2] ) / parameters.dx;
//
//    real dUdz = c1o2 * ( U[5] - U[4] ) / parameters.dx;
//    real dVdz = c1o2 * ( V[5] - V[4] ) / parameters.dx;
//
//    real tmpX = dWdy - dVdz;
//    real tmpY = dUdz - dWdx;
//    real tmpZ = dVdx - dUdy;
//
//    //////////////////////////////////////////////////////////////////////////
//
//    enstrophy[ cellIndex ] = c1o2 * rho[6] * ( tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ );
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EnstrophyAnalyzer::EnstrophyAnalyzer(SPtr<DataBase> dataBase, Parameters parameters, uint analyzeIter, uint outputIter)
{
    this->dataBase   = dataBase;
    this->parameters = parameters;

    this->analyzeIter = analyzeIter;
    this->outputIter  = outputIter;
}

void EnstrophyAnalyzer::writeToFile(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "EnstrophyAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + ".dat" );

    for( auto& EKin : this->enstrophyTimeSeries )
        file << std::setprecision(15) << EKin << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

} // namespace GksGpu


