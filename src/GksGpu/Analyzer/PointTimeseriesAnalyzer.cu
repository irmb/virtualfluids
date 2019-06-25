#include "PointTimeSeriesAnalyzer.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cmath>
#include <sstream>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

#include <iomanip>

#include "Core/Logger/Logger.h"

#include "DataBase/DataBase.h"

#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void pointTimeSeriesKernel  ( DataBaseStruct dataBase, PointTimeSeriesAnalyzerStruct pointTimeSeriesAnalyzer, Parameters parameters, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void pointTimeSeriesFunction( DataBaseStruct dataBase, PointTimeSeriesAnalyzerStruct pointTimeSeriesAnalyzer, Parameters parameters, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool PointTimeSeriesAnalyzer::run(uint iter, Parameters parameters)
{

    CudaUtility::CudaGrid grid( 1, 1 );

    runKernel( pointTimeSeriesKernel,
               pointTimeSeriesFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               this->toStruct(),
               parameters,
               0 );

    getLastCudaError("PointTimeSeriesAnalyzer::run(uint iter, Parameters parameters)");

    this->counter++;

    if( this->counter == this->outputIter )
    {
        this->download();
        this->counter = 0;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void pointTimeSeriesKernel(DataBaseStruct dataBase, PointTimeSeriesAnalyzerStruct pointTimeSeriesAnalyzer, Parameters parameters, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    pointTimeSeriesFunction( dataBase, pointTimeSeriesAnalyzer, parameters, startIndex, index );
}

__host__ __device__ void pointTimeSeriesFunction(DataBaseStruct dataBase, PointTimeSeriesAnalyzerStruct pointTimeSeriesAnalyzer, Parameters parameters, uint startIndex, uint index)
{
    //////////////////////////////////////////////////////////////////////////

    ConservedVariables cons;

    readCellData(pointTimeSeriesAnalyzer.cellIndex, dataBase, cons);

    PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

    //////////////////////////////////////////////////////////////////////////

    if( pointTimeSeriesAnalyzer.quantity == 'U' ) pointTimeSeriesAnalyzer.deviceSeries [ pointTimeSeriesAnalyzer.counter ] = prim.U;
    if( pointTimeSeriesAnalyzer.quantity == 'V' ) pointTimeSeriesAnalyzer.deviceSeries [ pointTimeSeriesAnalyzer.counter ] = prim.V;
    if( pointTimeSeriesAnalyzer.quantity == 'W' ) pointTimeSeriesAnalyzer.deviceSeries [ pointTimeSeriesAnalyzer.counter ] = prim.W;

#ifdef USE_PASSIVE_SCALAR
    if( pointTimeSeriesAnalyzer.quantity == 'T' ) pointTimeSeriesAnalyzer.deviceSeries [ pointTimeSeriesAnalyzer.counter ] = getT(prim);
#else
    if( pointTimeSeriesAnalyzer.quantity == 'T' ) pointTimeSeriesAnalyzer.deviceSeries [ pointTimeSeriesAnalyzer.counter ] = one / prim.lambda;
#endif

    //////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PointTimeSeriesAnalyzer::~PointTimeSeriesAnalyzer()
{
    this->free();
}

PointTimeSeriesAnalyzer::PointTimeSeriesAnalyzer(SPtr<DataBase> dataBase, GksMeshAdapter & adapter, Vec3 coordinate, char quantity, uint outputIter)
    : dataBase(dataBase),
      deviceSeries(nullptr),
      counter(0),
      outputIter(outputIter),
      quantity(quantity)
{
    this->allocate();

    this->findCellIndex( adapter, coordinate );
}

void PointTimeSeriesAnalyzer::free()
{
    checkCudaErrors( cudaFree ( this->deviceSeries  ) );
}

void PointTimeSeriesAnalyzer::allocate()
{
    this->free();

    checkCudaErrors( cudaMalloc ( &this->deviceSeries , sizeof(real) * this->outputIter ) );
}

void PointTimeSeriesAnalyzer::findCellIndex(GksMeshAdapter & adapter, Vec3 coordinate)
{
    real minDistance = 1.0e99;

    for( uint level = 0; level < adapter.numberOfLevels; level++ )
    {
        uint startIdx = adapter.startOfCellsPerLevel[level] 
                      + adapter.numberOfBulkCellsPerLevel[level];

        uint endIdx   = adapter.startOfCellsPerLevel[level] 
                      + adapter.numberOfCellsPerLevel[level];

        for( uint cellIdx = startIdx ; cellIdx < endIdx; cellIdx++ )
        {
            MeshCell& cell = adapter.cells[ cellIdx ];

            Vec3 vec = cell.cellCenter - coordinate;

            real distance = sqrt( vec.x*vec.x + vec.z*vec.z + vec.z*vec.z );

            if( distance < minDistance )
            {
                this->cellIndex = cellIndex;
                minDistance = distance;
            }
        }
    }
}

void PointTimeSeriesAnalyzer::writeToFile(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "PointTimeSeriesAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + ".dat" );

    for( auto& value : this->hostSeries )
        file << std::setprecision(15) << value << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

PointTimeSeriesAnalyzerStruct PointTimeSeriesAnalyzer::toStruct()
{
    PointTimeSeriesAnalyzerStruct pointTimeSeriesAnalyzer;

    pointTimeSeriesAnalyzer.deviceSeries = this->deviceSeries;

    pointTimeSeriesAnalyzer.quantity     = this->quantity;

    pointTimeSeriesAnalyzer.counter      = this->counter;

    pointTimeSeriesAnalyzer.cellIndex    = this->cellIndex;

    return pointTimeSeriesAnalyzer;
}

void PointTimeSeriesAnalyzer::download()
{
    uint oldSize = hostSeries.size();

    this->hostSeries.resize( oldSize + this->outputIter, zero );

    checkCudaErrors( cudaMemcpy( this->hostSeries.data() + oldSize, this->deviceSeries , sizeof(real) * outputIter, cudaMemcpyDeviceToHost ) );
}


