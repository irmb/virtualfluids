#include "TurbulenceAnalyzer.h"

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

#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void turbulenceKernel  ( DataBaseStruct dataBase, TurbulenceAnalyzerStruct turbulenceAnalyzer, Parameters parameters, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void turbulenceFunction( DataBaseStruct dataBase, TurbulenceAnalyzerStruct turbulenceAnalyzer, Parameters parameters, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool TurbulenceAnalyzer::run(uint iter, Parameters parameters)
{
    if( iter < this->analyzeStartIter ) return false;

    thrust::device_vector<real> kineticEnergy( this->dataBase->perLevelCount[ 0 ].numberOfBulkCells );

    CudaUtility::CudaGrid grid( dataBase->numberOfCells, 32 );

    runKernel( turbulenceKernel,
               turbulenceFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               this->toStruct(),
               parameters,
               0 );

    getLastCudaError("TurbulenceAnalyzer::run(uint iter, Parameters parameters)");

    this->counter++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void turbulenceKernel(DataBaseStruct dataBase, TurbulenceAnalyzerStruct turbulenceAnalyzer, Parameters parameters, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    turbulenceFunction( dataBase, turbulenceAnalyzer, parameters, startIndex, index );
}

__host__ __device__ void turbulenceFunction(DataBaseStruct dataBase, TurbulenceAnalyzerStruct turbulenceAnalyzer, Parameters parameters, uint startIndex, uint index)
{
    uint cellIndex = startIndex + index;

    //////////////////////////////////////////////////////////////////////////

    ConservedVariables cons;

    cons.rho  = dataBase.data[ RHO__(cellIndex, dataBase.numberOfCells) ];
    cons.rhoU = dataBase.data[ RHO_U(cellIndex, dataBase.numberOfCells) ];
    cons.rhoV = dataBase.data[ RHO_V(cellIndex, dataBase.numberOfCells) ];
    cons.rhoW = dataBase.data[ RHO_W(cellIndex, dataBase.numberOfCells) ];
    cons.rhoE = dataBase.data[ RHO_E(cellIndex, dataBase.numberOfCells) ];
#ifdef USE_PASSIVE_SCALAR
    cons.rhoS = dataBase.data[ RHO_S(cellIndex, dataBase.numberOfCells) ];
#endif // USE_PASSIVE_SCALAR

    PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

    //////////////////////////////////////////////////////////////////////////

    turbulenceAnalyzer.U [ cellIndex ] += prim.U;
    turbulenceAnalyzer.V [ cellIndex ] += prim.V;
    turbulenceAnalyzer.W [ cellIndex ] += prim.W;

    turbulenceAnalyzer.UU[ cellIndex ] += prim.U * prim.U;
    turbulenceAnalyzer.VV[ cellIndex ] += prim.V * prim.V;
    turbulenceAnalyzer.WW[ cellIndex ] += prim.W * prim.W;

    turbulenceAnalyzer.UV[ cellIndex ] += prim.U * prim.V;
    turbulenceAnalyzer.UW[ cellIndex ] += prim.U * prim.W;
    turbulenceAnalyzer.VW[ cellIndex ] += prim.V * prim.W;

    turbulenceAnalyzer.T [ cellIndex ] += one / prim.lambda;
    turbulenceAnalyzer.T [ cellIndex ] += c1o2 * prim.rho / prim.lambda;

    //////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TurbulenceAnalyzer::~TurbulenceAnalyzer()
{
    checkCudaErrors( cudaFree ( this->U  ) );
    checkCudaErrors( cudaFree ( this->V  ) );
    checkCudaErrors( cudaFree ( this->W  ) );
    checkCudaErrors( cudaFree ( this->UU ) );
    checkCudaErrors( cudaFree ( this->VV ) );
    checkCudaErrors( cudaFree ( this->WW ) );
    checkCudaErrors( cudaFree ( this->UV ) );
    checkCudaErrors( cudaFree ( this->UW ) );
    checkCudaErrors( cudaFree ( this->VW ) );
    checkCudaErrors( cudaFree ( this->T  ) );
    checkCudaErrors( cudaFree ( this->p  ) );
}

TurbulenceAnalyzer::TurbulenceAnalyzer(SPtr<DataBase> dataBase, uint analyzeStartIter)
    : U ( nullptr ),
      V ( nullptr ),
      W ( nullptr ),
      UU( nullptr ),
      VV( nullptr ),
      WW( nullptr ),
      UV( nullptr ),
      UW( nullptr ),
      VW( nullptr ),
      T ( nullptr ),
      p ( nullptr )
{
    this->dataBase = dataBase;

    this->analyzeStartIter = analyzeStartIter;

    this->counter = 0;

    checkCudaErrors( cudaMalloc ( &this->U , sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->V , sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->W , sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->UU, sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->VV, sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->WW, sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->UV, sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->UW, sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->VW, sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->T , sizeof(real) * dataBase->numberOfCells ) );
    checkCudaErrors( cudaMalloc ( &this->p , sizeof(real) * dataBase->numberOfCells ) );

    h_U.resize( dataBase->numberOfCells );
    h_V.resize( dataBase->numberOfCells ); 
    h_W.resize( dataBase->numberOfCells );
    h_UU.resize( dataBase->numberOfCells );
    h_VV.resize( dataBase->numberOfCells );
    h_WW.resize( dataBase->numberOfCells );
    h_UV.resize( dataBase->numberOfCells );
    h_UW.resize( dataBase->numberOfCells );
    h_VW.resize( dataBase->numberOfCells );
    h_T.resize( dataBase->numberOfCells );
    h_p.resize( dataBase->numberOfCells );
}

void TurbulenceAnalyzer::writeToFile(std::string filename)
{
    //writeTurbulenceVtkXML( this->dataBase, this->toHostStruct(), 0, filename );
}

TurbulenceAnalyzerStruct TurbulenceAnalyzer::toStruct()
{
    TurbulenceAnalyzerStruct turbulenceAnalyzer;

    turbulenceAnalyzer.U  = this->U;
    turbulenceAnalyzer.V  = this->V;
    turbulenceAnalyzer.W  = this->W;

    turbulenceAnalyzer.UU = this->UU;
    turbulenceAnalyzer.VV = this->VV;
    turbulenceAnalyzer.WW = this->WW;

    turbulenceAnalyzer.UV = this->UV;
    turbulenceAnalyzer.UW = this->UW;
    turbulenceAnalyzer.VW = this->VW;

    turbulenceAnalyzer.T  = this->T;
    turbulenceAnalyzer.p  = this->p;

    return turbulenceAnalyzer;
}

void TurbulenceAnalyzer::download()
{
    checkCudaErrors( cudaMemcpy( this->h_U.data() , this->U , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_V.data() , this->V , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_W.data() , this->W , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_UU.data(), this->UU, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_VV.data(), this->VV, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_WW.data(), this->WW, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_UV.data(), this->UV, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_UW.data(), this->UW, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_VW.data(), this->VW, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_T.data() , this->T , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy( this->h_p.data() , this->p , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
}


