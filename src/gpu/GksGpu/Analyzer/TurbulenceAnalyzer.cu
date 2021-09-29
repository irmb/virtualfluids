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

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

namespace GksGpu {

__global__                 void turbulenceKernel  ( DataBaseStruct dataBase, TurbulenceAnalyzerStruct turbulenceAnalyzer, Parameters parameters, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void turbulenceFunction( DataBaseStruct dataBase, TurbulenceAnalyzerStruct turbulenceAnalyzer, Parameters parameters, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool TurbulenceAnalyzer::run(uint iter, Parameters parameters)
{
    if( iter < this->analyzeStartIter ) return false;

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

    return true;
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

    readCellData(cellIndex, dataBase, cons);

    PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

    //////////////////////////////////////////////////////////////////////////

    if( turbulenceAnalyzer.U  ) turbulenceAnalyzer.U [ cellIndex ] += prim.U;
    if( turbulenceAnalyzer.V  ) turbulenceAnalyzer.V [ cellIndex ] += prim.V;
    if( turbulenceAnalyzer.W  ) turbulenceAnalyzer.W [ cellIndex ] += prim.W;

    if( turbulenceAnalyzer.UU ) turbulenceAnalyzer.UU[ cellIndex ] += prim.U * prim.U;
    if( turbulenceAnalyzer.VV ) turbulenceAnalyzer.VV[ cellIndex ] += prim.V * prim.V;
    if( turbulenceAnalyzer.WW ) turbulenceAnalyzer.WW[ cellIndex ] += prim.W * prim.W;

    if( turbulenceAnalyzer.UV ) turbulenceAnalyzer.UV[ cellIndex ] += prim.U * prim.V;
    if( turbulenceAnalyzer.UW ) turbulenceAnalyzer.UW[ cellIndex ] += prim.U * prim.W;
    if( turbulenceAnalyzer.VW ) turbulenceAnalyzer.VW[ cellIndex ] += prim.V * prim.W;

#ifdef USE_PASSIVE_SCALAR
    if( turbulenceAnalyzer.T  ) turbulenceAnalyzer.T [ cellIndex ] += getT(prim);
#else
    if( turbulenceAnalyzer.T  ) turbulenceAnalyzer.T [ cellIndex ] +=   c1o1 / prim.lambda;
#endif

    if( turbulenceAnalyzer.TT ) turbulenceAnalyzer.TT[ cellIndex ] += ( c1o1 / prim.lambda ) * ( c1o1 / prim.lambda );
    if( turbulenceAnalyzer.p  ) turbulenceAnalyzer.p [ cellIndex ] += c1o2 * prim.rho / prim.lambda;

    //////////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TurbulenceAnalyzer::~TurbulenceAnalyzer()
{
    this->free();
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
      TT( nullptr ),
      p ( nullptr ),
      collect_U ( true  ),
      collect_V ( true  ),
      collect_W ( true  ),
      collect_UU( false ),
      collect_VV( false ),
      collect_WW( false ),
      collect_UV( false ),
      collect_UW( false ),
      collect_VW( false ),
      collect_T ( true  ),
      collect_TT( false ),
      collect_p ( false )
{
    this->dataBase = dataBase;

    this->analyzeStartIter = analyzeStartIter;

    this->counter = 0;

    this->allocate();
}

void TurbulenceAnalyzer::free()
{
    if( this->U  ) checkCudaErrors( cudaFree ( this->U  ) );
    if( this->V  ) checkCudaErrors( cudaFree ( this->V  ) );
    if( this->W  ) checkCudaErrors( cudaFree ( this->W  ) );
    if( this->UU ) checkCudaErrors( cudaFree ( this->UU ) );
    if( this->VV ) checkCudaErrors( cudaFree ( this->VV ) );
    if( this->WW ) checkCudaErrors( cudaFree ( this->WW ) );
    if( this->UV ) checkCudaErrors( cudaFree ( this->UV ) );
    if( this->UW ) checkCudaErrors( cudaFree ( this->UW ) );
    if( this->VW ) checkCudaErrors( cudaFree ( this->VW ) );
    if( this->T  ) checkCudaErrors( cudaFree ( this->T  ) );
    if( this->TT ) checkCudaErrors( cudaFree ( this->TT ) );
    if( this->p  ) checkCudaErrors( cudaFree ( this->p  ) );

    h_U.clear ( );
    h_V.clear ( );
    h_W.clear ( );
    h_UU.clear( );
    h_VV.clear( );
    h_WW.clear( );
    h_UV.clear( );
    h_UW.clear( );
    h_VW.clear( );
    h_T.clear ( );
    h_TT.clear( );
    h_p.clear ( );
}

void TurbulenceAnalyzer::allocate()
{
    this->free();

    if( collect_U  ) checkCudaErrors( cudaMalloc ( &this->U , sizeof(real) * dataBase->numberOfCells ) );
    if( collect_V  ) checkCudaErrors( cudaMalloc ( &this->V , sizeof(real) * dataBase->numberOfCells ) );
    if( collect_W  ) checkCudaErrors( cudaMalloc ( &this->W , sizeof(real) * dataBase->numberOfCells ) );
    if( collect_UU ) checkCudaErrors( cudaMalloc ( &this->UU, sizeof(real) * dataBase->numberOfCells ) );
    if( collect_VV ) checkCudaErrors( cudaMalloc ( &this->VV, sizeof(real) * dataBase->numberOfCells ) );
    if( collect_WW ) checkCudaErrors( cudaMalloc ( &this->WW, sizeof(real) * dataBase->numberOfCells ) );
    if( collect_UV ) checkCudaErrors( cudaMalloc ( &this->UV, sizeof(real) * dataBase->numberOfCells ) );
    if( collect_UW ) checkCudaErrors( cudaMalloc ( &this->UW, sizeof(real) * dataBase->numberOfCells ) );
    if( collect_VW ) checkCudaErrors( cudaMalloc ( &this->VW, sizeof(real) * dataBase->numberOfCells ) );
    if( collect_T  ) checkCudaErrors( cudaMalloc ( &this->T , sizeof(real) * dataBase->numberOfCells ) );
    if( collect_TT ) checkCudaErrors( cudaMalloc ( &this->TT, sizeof(real) * dataBase->numberOfCells ) );
    if( collect_p  ) checkCudaErrors( cudaMalloc ( &this->p , sizeof(real) * dataBase->numberOfCells ) );

    if( collect_U  ) h_U.resize ( dataBase->numberOfCells );
    if( collect_V  ) h_V.resize ( dataBase->numberOfCells ); 
    if( collect_W  ) h_W.resize ( dataBase->numberOfCells );
    if( collect_UU ) h_UU.resize( dataBase->numberOfCells );
    if( collect_VV ) h_VV.resize( dataBase->numberOfCells );
    if( collect_WW ) h_WW.resize( dataBase->numberOfCells );
    if( collect_UV ) h_UV.resize( dataBase->numberOfCells );
    if( collect_UW ) h_UW.resize( dataBase->numberOfCells );
    if( collect_VW ) h_VW.resize( dataBase->numberOfCells );
    if( collect_T  ) h_T.resize ( dataBase->numberOfCells );
    if( collect_TT ) h_TT.resize( dataBase->numberOfCells );
    if( collect_p  ) h_p.resize ( dataBase->numberOfCells );
}

void TurbulenceAnalyzer::writeRestartFile(std::string filename)
{
    this->download(false);

    //////////////////////////////////////////////////////////////////////////

    filename += ".rst";

    *logging::out << logging::Logger::INFO_HIGH << "Writing restart file " << filename << " ... ";

    std::ofstream file;

	file.open( filename.c_str(), std::ios::binary );

	if (!file.is_open()) {
		throw std::runtime_error("\nFile cannot be opened.\n\nERROR!\n\n\n");
        return;
	}

    //////////////////////////////////////////////////////////////////////////

    file.write( (char*) &this->counter, sizeof( uint ) );

    file.write( (char*) &dataBase->numberOfLevels, sizeof( uint ) );
    file.write( (char*) &dataBase->numberOfCells,  sizeof( uint ) );
    file.write( (char*) &dataBase->numberOfFaces,  sizeof( uint ) );

    file.write( (char*) &this->collect_U , sizeof( bool ) );
    file.write( (char*) &this->collect_V , sizeof( bool ) );
    file.write( (char*) &this->collect_W , sizeof( bool ) );
    file.write( (char*) &this->collect_UU, sizeof( bool ) );
    file.write( (char*) &this->collect_VV, sizeof( bool ) );
    file.write( (char*) &this->collect_WW, sizeof( bool ) );
    file.write( (char*) &this->collect_UV, sizeof( bool ) );
    file.write( (char*) &this->collect_UW, sizeof( bool ) );
    file.write( (char*) &this->collect_VW, sizeof( bool ) );
    file.write( (char*) &this->collect_T , sizeof( bool ) );
    file.write( (char*) &this->collect_TT, sizeof( bool ) );
    file.write( (char*) &this->collect_p , sizeof( bool ) );

    if( collect_U  ) file.write( (char*) this->h_U.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_V  ) file.write( (char*) this->h_V.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_W  ) file.write( (char*) this->h_W.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_UU ) file.write( (char*) this->h_UU.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_VV ) file.write( (char*) this->h_VV.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_WW ) file.write( (char*) this->h_WW.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_UV ) file.write( (char*) this->h_UV.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_UW ) file.write( (char*) this->h_UW.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_VW ) file.write( (char*) this->h_VW.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_T  ) file.write( (char*) this->h_T.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_TT ) file.write( (char*) this->h_TT.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_p  ) file.write( (char*) this->h_p.data() , dataBase->numberOfCells * sizeof( real ) );

    //////////////////////////////////////////////////////////////////////////

    file.close();

    *logging::out << logging::Logger::INFO_HIGH << "done!\n";
}

void TurbulenceAnalyzer::readRestartFile(std::string filename)
{
    filename += ".rst";

    *logging::out << logging::Logger::INFO_HIGH << "Reading restart file " << filename << " ... ";
	
    std::ifstream file;

	file.open( filename.c_str(), std::ios::binary );

	if (!file.is_open()) {
		throw std::runtime_error("\nFile cannot be opened.\n\nERROR!\n\n\n");
        return;
	}

    //////////////////////////////////////////////////////////////////////////

    file.read( (char*) &this->counter, sizeof( uint ) );

    uint numberOfLevelsRead;
    uint numberOfCellsRead;
    uint numberOfFacesRead;
    
    file.read( (char*) &numberOfLevelsRead, sizeof( uint ) );
    file.read( (char*) &numberOfCellsRead,  sizeof( uint ) );
    file.read( (char*) &numberOfFacesRead,  sizeof( uint ) );

    if( numberOfLevelsRead != dataBase->numberOfLevels ||
        numberOfCellsRead  != dataBase->numberOfCells  ||
        numberOfFacesRead  != dataBase->numberOfFaces  ){
    
        *logging::out << logging::Logger::INFO_HIGH << "\n";
        *logging::out << logging::Logger::INFO_HIGH << "Levels: " << numberOfLevelsRead << " vs. " << dataBase->numberOfLevels << "\n";
        *logging::out << logging::Logger::INFO_HIGH << "Cells:  " << numberOfCellsRead  << " vs. " << dataBase->numberOfCells  << "\n";
        *logging::out << logging::Logger::INFO_HIGH << "Faces:  " << numberOfFacesRead  << " vs. " << dataBase->numberOfFaces  << "\n";

        file.close();

        throw std::runtime_error("\nERROR: Restart file does not match current setup");
    }

    //////////////////////////////////////////////////////////////////////////

    file.read( (char*) &this->collect_U , sizeof( bool ) );
    file.read( (char*) &this->collect_V , sizeof( bool ) );
    file.read( (char*) &this->collect_W , sizeof( bool ) );
    file.read( (char*) &this->collect_UU, sizeof( bool ) );
    file.read( (char*) &this->collect_VV, sizeof( bool ) );
    file.read( (char*) &this->collect_WW, sizeof( bool ) );
    file.read( (char*) &this->collect_UV, sizeof( bool ) );
    file.read( (char*) &this->collect_UW, sizeof( bool ) );
    file.read( (char*) &this->collect_VW, sizeof( bool ) );
    file.read( (char*) &this->collect_T , sizeof( bool ) );
    file.read( (char*) &this->collect_TT, sizeof( bool ) );
    file.read( (char*) &this->collect_p , sizeof( bool ) );

    if( collect_U  ) file.read( (char*) this->h_U.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_V  ) file.read( (char*) this->h_V.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_W  ) file.read( (char*) this->h_W.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_UU ) file.read( (char*) this->h_UU.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_VV ) file.read( (char*) this->h_VV.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_WW ) file.read( (char*) this->h_WW.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_UV ) file.read( (char*) this->h_UV.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_UW ) file.read( (char*) this->h_UW.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_VW ) file.read( (char*) this->h_VW.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_T  ) file.read( (char*) this->h_T.data() , dataBase->numberOfCells * sizeof( real ) );
    if( collect_TT ) file.read( (char*) this->h_TT.data(), dataBase->numberOfCells * sizeof( real ) );
    if( collect_p  ) file.read( (char*) this->h_p.data() , dataBase->numberOfCells * sizeof( real ) );

    //////////////////////////////////////////////////////////////////////////

    file.close();

    this->upload();

    *logging::out << logging::Logger::INFO_HIGH << "done!\n";
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
    turbulenceAnalyzer.TT = this->TT;
    turbulenceAnalyzer.p  = this->p;

    return turbulenceAnalyzer;
}

void TurbulenceAnalyzer::download(bool normalize)
{
    if( collect_U  ) checkCudaErrors( cudaMemcpy( this->h_U.data() , this->U , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_V  ) checkCudaErrors( cudaMemcpy( this->h_V.data() , this->V , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_W  ) checkCudaErrors( cudaMemcpy( this->h_W.data() , this->W , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_UU ) checkCudaErrors( cudaMemcpy( this->h_UU.data(), this->UU, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_VV ) checkCudaErrors( cudaMemcpy( this->h_VV.data(), this->VV, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_WW ) checkCudaErrors( cudaMemcpy( this->h_WW.data(), this->WW, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_UV ) checkCudaErrors( cudaMemcpy( this->h_UV.data(), this->UV, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_UW ) checkCudaErrors( cudaMemcpy( this->h_UW.data(), this->UW, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_VW ) checkCudaErrors( cudaMemcpy( this->h_VW.data(), this->VW, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_T  ) checkCudaErrors( cudaMemcpy( this->h_T.data() , this->T , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_TT ) checkCudaErrors( cudaMemcpy( this->h_TT.data(), this->TT, sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );
    if( collect_p  ) checkCudaErrors( cudaMemcpy( this->h_p.data() , this->p , sizeof(real) * dataBase->numberOfCells, cudaMemcpyDeviceToHost ) );

    if(normalize)
    {
        for( uint cellIndex = 0; cellIndex < dataBase->numberOfCells; cellIndex++ )
        {
            if( collect_U  ) this->h_U [ cellIndex ] /= real(this->counter);
            if( collect_V  ) this->h_V [ cellIndex ] /= real(this->counter);
            if( collect_W  ) this->h_W [ cellIndex ] /= real(this->counter);
            if( collect_UU ) this->h_UU[ cellIndex ] /= real(this->counter);
            if( collect_VV ) this->h_VV[ cellIndex ] /= real(this->counter);
            if( collect_WW ) this->h_WW[ cellIndex ] /= real(this->counter);
            if( collect_UV ) this->h_UV[ cellIndex ] /= real(this->counter);
            if( collect_UW ) this->h_UW[ cellIndex ] /= real(this->counter);
            if( collect_VW ) this->h_VW[ cellIndex ] /= real(this->counter);
            if( collect_T  ) this->h_T [ cellIndex ] /= real(this->counter);
            if( collect_TT ) this->h_TT[ cellIndex ] /= real(this->counter);
            if( collect_p  ) this->h_p [ cellIndex ] /= real(this->counter);

            if( collect_UU ) this->h_UU[ cellIndex ] -= this->h_U[ cellIndex ] * this->h_U[ cellIndex ];
            if( collect_VV ) this->h_VV[ cellIndex ] -= this->h_V[ cellIndex ] * this->h_V[ cellIndex ];
            if( collect_WW ) this->h_WW[ cellIndex ] -= this->h_W[ cellIndex ] * this->h_W[ cellIndex ];

            if( collect_UV ) this->h_UV[ cellIndex ] -= this->h_U[ cellIndex ] * this->h_V[ cellIndex ];
            if( collect_UW ) this->h_UW[ cellIndex ] -= this->h_U[ cellIndex ] * this->h_W[ cellIndex ];
            if( collect_VW ) this->h_VW[ cellIndex ] -= this->h_V[ cellIndex ] * this->h_W[ cellIndex ];
        
            if( collect_TT ) this->h_TT[ cellIndex ] -= this->h_T[ cellIndex ] * this->h_T[ cellIndex ];
        }
    }
}

void TurbulenceAnalyzer::upload()
{
    if( collect_U  ) checkCudaErrors( cudaMemcpy( this->U , this->h_U.data() , sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_V  ) checkCudaErrors( cudaMemcpy( this->V , this->h_V.data() , sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_W  ) checkCudaErrors( cudaMemcpy( this->W , this->h_W.data() , sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_UU ) checkCudaErrors( cudaMemcpy( this->UU, this->h_UU.data(), sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_VV ) checkCudaErrors( cudaMemcpy( this->VV, this->h_VV.data(), sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_WW ) checkCudaErrors( cudaMemcpy( this->WW, this->h_WW.data(), sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_UV ) checkCudaErrors( cudaMemcpy( this->UV, this->h_UV.data(), sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_UW ) checkCudaErrors( cudaMemcpy( this->UW, this->h_UW.data(), sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_VW ) checkCudaErrors( cudaMemcpy( this->VW, this->h_VW.data(), sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_T  ) checkCudaErrors( cudaMemcpy( this->T , this->h_T.data() , sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_TT ) checkCudaErrors( cudaMemcpy( this->TT, this->h_TT.data(), sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
    if( collect_p  ) checkCudaErrors( cudaMemcpy( this->p , this->h_p.data() , sizeof(real) * dataBase->numberOfCells, cudaMemcpyHostToDevice ) );
}

} // namespace GksGpu


