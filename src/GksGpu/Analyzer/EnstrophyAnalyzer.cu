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

__global__                 void enstrophyKernel  ( DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void enstrophyFunction( DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool EnstrophyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

    thrust::device_vector<real> enstrophy( this->dataBase->perLevelCount[ 0 ].numberOfBulkCells );

    CudaUtility::CudaGrid grid( dataBase->perLevelCount[ 0 ].numberOfBulkCells, 32 );

    runKernel( enstrophyKernel,
               enstrophyFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               parameters,
               enstrophy.data().get(),
               dataBase->perLevelCount[ 0 ].startOfCells );

    getLastCudaError("KineticEnergyAnalyzer::run(uint iter)");

    real EnstrophyTmp = thrust::reduce( enstrophy.begin(), enstrophy.end(), zero, thrust::plus<real>() )
                      / real(dataBase->perLevelCount[ 0 ].numberOfBulkCells);

    this->enstrophyTimeSeries.push_back( EnstrophyTmp );

    //*logging::out << logging::Logger::INFO_HIGH << "EKin = " << EKin << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void enstrophyKernel(DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    enstrophyFunction( dataBase, parameters, enstrophy, startIndex, index );
}

__host__ __device__ void enstrophyFunction(DataBaseStruct dataBase, Parameters parameters, real* enstrophy, uint startIndex, uint index)
{
    uint cellIndex = startIndex + index;

    //////////////////////////////////////////////////////////////////////////

    uint cellToCell [6];

    cellToCell[0] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 0, dataBase.numberOfCells ) ];
    cellToCell[1] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 1, dataBase.numberOfCells ) ];
    cellToCell[2] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 2, dataBase.numberOfCells ) ];
    cellToCell[3] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 3, dataBase.numberOfCells ) ];
    cellToCell[4] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 4, dataBase.numberOfCells ) ];
    cellToCell[5] = dataBase.cellToCell[ CELL_TO_CELL( cellIndex, 5, dataBase.numberOfCells ) ];

    real rho [7];
    real U   [6];
    real V   [6];
    real W   [6];

    rho[0] = dataBase.data[ RHO__( cellToCell[0], dataBase.numberOfCells ) ];
    rho[1] = dataBase.data[ RHO__( cellToCell[1], dataBase.numberOfCells ) ];
    rho[2] = dataBase.data[ RHO__( cellToCell[2], dataBase.numberOfCells ) ];
    rho[3] = dataBase.data[ RHO__( cellToCell[3], dataBase.numberOfCells ) ];
    rho[4] = dataBase.data[ RHO__( cellToCell[4], dataBase.numberOfCells ) ];
    rho[5] = dataBase.data[ RHO__( cellToCell[5], dataBase.numberOfCells ) ];
    rho[6] = dataBase.data[ RHO__( cellIndex    , dataBase.numberOfCells ) ];

    U  [0] = dataBase.data[ RHO_U( cellToCell[0], dataBase.numberOfCells ) ] / rho[0];
    U  [1] = dataBase.data[ RHO_U( cellToCell[1], dataBase.numberOfCells ) ] / rho[1];
    U  [2] = dataBase.data[ RHO_U( cellToCell[2], dataBase.numberOfCells ) ] / rho[2];
    U  [3] = dataBase.data[ RHO_U( cellToCell[3], dataBase.numberOfCells ) ] / rho[3];
    U  [4] = dataBase.data[ RHO_U( cellToCell[4], dataBase.numberOfCells ) ] / rho[4];
    U  [5] = dataBase.data[ RHO_U( cellToCell[5], dataBase.numberOfCells ) ] / rho[5];

    V  [0] = dataBase.data[ RHO_V( cellToCell[0], dataBase.numberOfCells ) ] / rho[0];
    V  [1] = dataBase.data[ RHO_V( cellToCell[1], dataBase.numberOfCells ) ] / rho[1];
    V  [2] = dataBase.data[ RHO_V( cellToCell[2], dataBase.numberOfCells ) ] / rho[2];
    V  [3] = dataBase.data[ RHO_V( cellToCell[3], dataBase.numberOfCells ) ] / rho[3];
    V  [4] = dataBase.data[ RHO_V( cellToCell[4], dataBase.numberOfCells ) ] / rho[4];
    V  [5] = dataBase.data[ RHO_V( cellToCell[5], dataBase.numberOfCells ) ] / rho[5];

    W  [0] = dataBase.data[ RHO_W( cellToCell[0], dataBase.numberOfCells ) ] / rho[0];
    W  [1] = dataBase.data[ RHO_W( cellToCell[1], dataBase.numberOfCells ) ] / rho[1];
    W  [2] = dataBase.data[ RHO_W( cellToCell[2], dataBase.numberOfCells ) ] / rho[2];
    W  [3] = dataBase.data[ RHO_W( cellToCell[3], dataBase.numberOfCells ) ] / rho[3];
    W  [4] = dataBase.data[ RHO_W( cellToCell[4], dataBase.numberOfCells ) ] / rho[4];
    W  [5] = dataBase.data[ RHO_W( cellToCell[5], dataBase.numberOfCells ) ] / rho[5];

    real dVdx = c1o2 * ( V[1] - V[0] ) / parameters.dx;
    real dWdx = c1o2 * ( W[1] - W[0] ) / parameters.dx;

    real dUdy = c1o2 * ( U[3] - U[2] ) / parameters.dx;
    real dWdy = c1o2 * ( W[3] - W[2] ) / parameters.dx;

    real dUdz = c1o2 * ( U[5] - U[4] ) / parameters.dx;
    real dVdz = c1o2 * ( V[5] - V[4] ) / parameters.dx;

    real tmpX = dWdy - dVdz;
    real tmpY = dUdz - dWdx;
    real tmpZ = dVdx - dUdy;

    //////////////////////////////////////////////////////////////////////////

    enstrophy[ cellIndex ] = c1o2 * rho[6] * ( tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ );
}

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
        file << EKin << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


