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

#include "Parameter/Parameter.h"
// includes, kernels
#include "GPU/GPU_Kernels.cuh"
#include "GPU/constant.h"

__global__                 void enstrophyKernel  ( real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* enstrophy, uint* isFluid, uint size_Mat );

__host__ __device__ inline void enstrophyFunction( real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB,            real* enstrophy, uint* isFluid, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool EnstrophyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

	int lev = 0;
	int size_Mat = this->para->getParD(lev)->size_Mat_SP;
	
	thrust::device_vector<real> enstrophy( size_Mat, zero );
    thrust::device_vector<uint> isFluid  ( size_Mat, 0);

	unsigned int numberOfThreads = 128;
    int Grid = (size_Mat / numberOfThreads)+1;
    int Grid1, Grid2;
    if (Grid>512)
    {
       Grid1 = 512;
       Grid2 = (Grid/Grid1)+1;
    } 
    else
    {
       Grid1 = 1;
       Grid2 = Grid;
    }
    dim3 grid(Grid1, Grid2);
    dim3 threads(numberOfThreads, 1, 1 );

    LBCalcMacCompSP27<<< grid, threads >>> (para->getParD(lev)->vx_SP,
										    para->getParD(lev)->vy_SP,
										    para->getParD(lev)->vz_SP,
										    para->getParD(lev)->rho_SP,
										    para->getParD(lev)->press_SP,
										    para->getParD(lev)->geoSP,
										    para->getParD(lev)->neighborX_SP,
										    para->getParD(lev)->neighborY_SP,
										    para->getParD(lev)->neighborZ_SP,
										    para->getParD(lev)->size_Mat_SP,
										    para->getParD(lev)->d0SP.f[0],
										    para->getParD(lev)->evenOrOdd); 
	//cudaDeviceSynchronize();
	getLastCudaError("LBCalcMacSP27 execution failed"); 

	enstrophyKernel <<< grid, threads >>> ( para->getParD(lev)->vx_SP,
											para->getParD(lev)->vy_SP, 
											para->getParD(lev)->vz_SP, 
											para->getParD(lev)->rho_SP, 
											para->getParD(lev)->neighborX_SP,
											para->getParD(lev)->neighborY_SP,
											para->getParD(lev)->neighborZ_SP,
											para->getParD(lev)->neighborWSB_SP,
											para->getParD(lev)->geoSP,
											enstrophy.data().get(), 
                                            isFluid.data().get(),
											size_Mat);
	cudaDeviceSynchronize(); 
	getLastCudaError("enstrophyKernel execution failed");

	real EnstrophyTmp       = thrust::reduce(enstrophy.begin(), enstrophy.end(), zero, thrust::plus<real>());
    uint numberOfFluidNodes = thrust::reduce(isFluid.begin(),   isFluid.end(),   0,    thrust::plus<uint>());

	this->enstrophyTimeSeries.push_back( EnstrophyTmp / real(numberOfFluidNodes) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void enstrophyKernel(real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* enstrophy, uint* isFluid, uint size_Mat)
{
    //////////////////////////////////////////////////////////////////////////
    const uint x = threadIdx.x;  // Globaler x-Index 
    const uint y = blockIdx.x;   // Globaler y-Index 
    const uint z = blockIdx.y;   // Globaler z-Index 

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint index = nx*(ny*z + y) + x;
	////////////////////////////////////////////////////////////////////////////////
    //printf("%d\n", index);

    if( index >= size_Mat) return;

	unsigned int BC;
	BC = geo[index];
	if (BC != GEO_FLUID) return;

    enstrophyFunction( veloX, veloY, veloZ, rho, neighborX, neighborY, neighborZ, neighborWSB, enstrophy, isFluid, index );
}

__host__ __device__ void enstrophyFunction(real* veloX, real* veloY, real* veloZ, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, real* enstrophy, uint* isFluid, uint index)
{
	//////////////////////////////////////////////////////////////////////////////
	//neighbor index
	uint k = index;
	uint kPx = neighborX[k];
	uint kPy = neighborY[k];
	uint kPz = neighborZ[k];
	uint kMxyz = neighborWSB[k];
	uint kMx = neighborZ[neighborY[kMxyz]];
	uint kMy = neighborZ[neighborX[kMxyz]];
	uint kMz = neighborY[neighborX[kMxyz]];
    //////////////////////////////////////////////////////////////////////////
    uint kPx2 = neighborX[kPx];
    uint kPy2 = neighborY[kPy];
    uint kPz2 = neighborZ[kPz];
    uint kMxWSB = neighborWSB[kMx];
    uint kMyWSB = neighborWSB[kMy];
    uint kMzWSB = neighborWSB[kMz];
    uint kMx2 = neighborZ[neighborY[kMxWSB]];
    uint kMy2 = neighborZ[neighborX[kMyWSB]];
    uint kMz2 = neighborY[neighborX[kMzWSB]];
    //////////////////////////////////////////////////////////////////////////
    uint kPx3 = neighborX[kPx2];
    uint kPy3 = neighborY[kPy2];
    uint kPz3 = neighborZ[kPz2];
    uint kMx2WSB = neighborWSB[kMx2];
    uint kMy2WSB = neighborWSB[kMy2];
    uint kMz2WSB = neighborWSB[kMz2];
    uint kMx3 = neighborZ[neighborY[kMx2WSB]];
    uint kMy3 = neighborZ[neighborX[kMy2WSB]];
    uint kMz3 = neighborY[neighborX[kMz2WSB]];
    //////////////////////////////////////////////////////////////////////////
    uint kPx4 = neighborX[kPx3];
    uint kPy4 = neighborY[kPy3];
    uint kPz4 = neighborZ[kPz3];
    uint kMx3WSB = neighborWSB[kMx3];
    uint kMy3WSB = neighborWSB[kMy3];
    uint kMz3WSB = neighborWSB[kMz3];
    uint kMx4 = neighborZ[neighborY[kMx3WSB]];
    uint kMy4 = neighborZ[neighborX[kMy3WSB]];
    uint kMz4 = neighborY[neighborX[kMz3WSB]];
    //////////////////////////////////////////////////////////////////////////
	//getVeloX//
	real veloXNeighborPx = veloX[kPx];
	real veloXNeighborMx = veloX[kMx];
	real veloXNeighborPy = veloX[kPy];
	real veloXNeighborMy = veloX[kMy];
	real veloXNeighborPz = veloX[kPz];
	real veloXNeighborMz = veloX[kMz];
	//getVeloY//
	real veloYNeighborPx = veloY[kPx];
	real veloYNeighborMx = veloY[kMx];
	real veloYNeighborPy = veloY[kPy];
	real veloYNeighborMy = veloY[kMy];
	real veloYNeighborPz = veloY[kPz];
	real veloYNeighborMz = veloY[kMz];
	//getVeloZ//
	real veloZNeighborPx = veloZ[kPx];
	real veloZNeighborMx = veloZ[kMx];
	real veloZNeighborPy = veloZ[kPy];
	real veloZNeighborMy = veloZ[kMy];
	real veloZNeighborPz = veloZ[kPz];
	real veloZNeighborMz = veloZ[kMz];
	//////////////////////////////////////////////////////////////////////////////
	//getVeloX//
	real veloXNeighborPx2 = veloX[kPx2];
	real veloXNeighborMx2 = veloX[kMx2];
	real veloXNeighborPy2 = veloX[kPy2];
	real veloXNeighborMy2 = veloX[kMy2];
	real veloXNeighborPz2 = veloX[kPz2];
	real veloXNeighborMz2 = veloX[kMz2];
	//getVeloY//
	real veloYNeighborPx2 = veloY[kPx2];
	real veloYNeighborMx2 = veloY[kMx2];
	real veloYNeighborPy2 = veloY[kPy2];
	real veloYNeighborMy2 = veloY[kMy2];
	real veloYNeighborPz2 = veloY[kPz2];
	real veloYNeighborMz2 = veloY[kMz2];
	//getVeloZ//
	real veloZNeighborPx2 = veloZ[kPx2];
	real veloZNeighborMx2 = veloZ[kMx2];
	real veloZNeighborPy2 = veloZ[kPy2];
	real veloZNeighborMy2 = veloZ[kMy2];
	real veloZNeighborPz2 = veloZ[kPz2];
	real veloZNeighborMz2 = veloZ[kMz2];
	//////////////////////////////////////////////////////////////////////////////
	//getVeloX//
	real veloXNeighborPx3 = veloX[kPx3];
	real veloXNeighborMx3 = veloX[kMx3];
	real veloXNeighborPy3 = veloX[kPy3];
	real veloXNeighborMy3 = veloX[kMy3];
	real veloXNeighborPz3 = veloX[kPz3];
	real veloXNeighborMz3 = veloX[kMz3];
	//getVeloY//
	real veloYNeighborPx3 = veloY[kPx3];
	real veloYNeighborMx3 = veloY[kMx3];
	real veloYNeighborPy3 = veloY[kPy3];
	real veloYNeighborMy3 = veloY[kMy3];
	real veloYNeighborPz3 = veloY[kPz3];
	real veloYNeighborMz3 = veloY[kMz3];
	//getVeloZ//
	real veloZNeighborPx3 = veloZ[kPx3];
	real veloZNeighborMx3 = veloZ[kMx3];
	real veloZNeighborPy3 = veloZ[kPy3];
	real veloZNeighborMy3 = veloZ[kMy3];
	real veloZNeighborPz3 = veloZ[kPz3];
	real veloZNeighborMz3 = veloZ[kMz3];
	//////////////////////////////////////////////////////////////////////////////
	//getVeloX//
	real veloXNeighborPx4 = veloX[kPx4];
	real veloXNeighborMx4 = veloX[kMx4];
	real veloXNeighborPy4 = veloX[kPy4];
	real veloXNeighborMy4 = veloX[kMy4];
	real veloXNeighborPz4 = veloX[kPz4];
	real veloXNeighborMz4 = veloX[kMz4];
	//getVeloY//
	real veloYNeighborPx4 = veloY[kPx4];
	real veloYNeighborMx4 = veloY[kMx4];
	real veloYNeighborPy4 = veloY[kPy4];
	real veloYNeighborMy4 = veloY[kMy4];
	real veloYNeighborPz4 = veloY[kPz4];
	real veloYNeighborMz4 = veloY[kMz4];
	//getVeloZ//
	real veloZNeighborPx4 = veloZ[kPx4];
	real veloZNeighborMx4 = veloZ[kMx4];
	real veloZNeighborPy4 = veloZ[kPy4];
	real veloZNeighborMy4 = veloZ[kMy4];
	real veloZNeighborPz4 = veloZ[kPz4];
	real veloZNeighborMz4 = veloZ[kMz4];
	//////////////////////////////////////////////////////////////////////////////
	real dxvx = zero;
	real dyvx = zero;
	real dzvx = zero;
	real dxvy = zero;
	real dyvy = zero;
	real dzvy = zero;
	real dxvz = zero;
	real dyvz = zero;
	real dzvz = zero;
    //////////////////////////////////////////////////////////////////////////

	//dxvy = (veloYNeighborPx - veloYNeighborMx) / two;
	//dxvz = (veloZNeighborPx - veloZNeighborMx) / two;

	//dyvx = (veloXNeighborPy - veloXNeighborMy) / two;
	//dyvz = (veloZNeighborPy - veloZNeighborMy) / two;

	//dzvx = (veloXNeighborPz - veloXNeighborMz) / two;
	//dzvy = (veloYNeighborPz - veloYNeighborMz) / two;

    //////////////////////////////////////////////////////////////////////////

	//dxvy = ( ( eight * veloYNeighborPx - eight * veloYNeighborMx ) - ( veloYNeighborPx2 - veloYNeighborMx2) ) / twelve;
	//dxvz = ( ( eight * veloZNeighborPx - eight * veloZNeighborMx ) - ( veloZNeighborPx2 - veloZNeighborMx2) ) / twelve;

	//dyvx = ( ( eight * veloXNeighborPy - eight * veloXNeighborMy ) - ( veloXNeighborPy2 - veloXNeighborMy2) ) / twelve;
	//dyvz = ( ( eight * veloZNeighborPy - eight * veloZNeighborMy ) - ( veloZNeighborPy2 - veloZNeighborMy2) ) / twelve;

	//dzvx = ( ( eight * veloXNeighborPz - eight * veloXNeighborMz ) - ( veloXNeighborPz2 - veloXNeighborMz2) ) / twelve;
	//dzvy = ( ( eight * veloYNeighborPz - eight * veloYNeighborMz ) - ( veloYNeighborPz2 - veloYNeighborMz2) ) / twelve;

    //////////////////////////////////////////////////////////////////////////

	//dxvy = ( (five * nine) * ( veloYNeighborPx - veloYNeighborMx ) - nine * ( veloYNeighborPx2 - veloYNeighborMx2) + ( veloYNeighborPx3 - veloYNeighborMx3) ) / (six * ten);
	//dxvz = ( (five * nine) * ( veloZNeighborPx - veloZNeighborMx ) - nine * ( veloZNeighborPx2 - veloZNeighborMx2) + ( veloZNeighborPx3 - veloZNeighborMx3) ) / (six * ten);

	//dyvx = ( (five * nine) * ( veloXNeighborPy - veloXNeighborMy ) - nine * ( veloXNeighborPy2 - veloXNeighborMy2) + ( veloXNeighborPy3 - veloXNeighborMy3) ) / (six * ten);
	//dyvz = ( (five * nine) * ( veloZNeighborPy - veloZNeighborMy ) - nine * ( veloZNeighborPy2 - veloZNeighborMy2) + ( veloZNeighborPy3 - veloZNeighborMy3) ) / (six * ten);

	//dzvx = ( (five * nine) * ( veloXNeighborPz - veloXNeighborMz ) - nine * ( veloXNeighborPz2 - veloXNeighborMz2) + ( veloXNeighborPz3 - veloXNeighborMz3) ) / (six * ten);
	//dzvy = ( (five * nine) * ( veloYNeighborPz - veloYNeighborMz ) - nine * ( veloYNeighborPz2 - veloYNeighborMz2) + ( veloYNeighborPz3 - veloYNeighborMz3) ) / (six * ten);

    //////////////////////////////////////////////////////////////////////////

	dxvy = ( (twentyeight * eight) * ( veloYNeighborPx - veloYNeighborMx ) - (seven * eight) * ( veloYNeighborPx2 - veloYNeighborMx2) + (eight * four * c1o3) * ( veloYNeighborPx3 - veloYNeighborMx3) - ( veloYNeighborPx4 - veloYNeighborMx4) ) / (seven * ten * four);
	dxvz = ( (twentyeight * eight) * ( veloZNeighborPx - veloZNeighborMx ) - (seven * eight) * ( veloZNeighborPx2 - veloZNeighborMx2) + (eight * four * c1o3) * ( veloZNeighborPx3 - veloZNeighborMx3) - ( veloZNeighborPx4 - veloZNeighborMx4) ) / (seven * ten * four);

	dyvx = ( (twentyeight * eight) * ( veloXNeighborPy - veloXNeighborMy ) - (seven * eight) * ( veloXNeighborPy2 - veloXNeighborMy2) + (eight * four * c1o3) * ( veloXNeighborPy3 - veloXNeighborMy3) - ( veloXNeighborPy4 - veloXNeighborMy4) ) / (seven * ten * four);
	dyvz = ( (twentyeight * eight) * ( veloZNeighborPy - veloZNeighborMy ) - (seven * eight) * ( veloZNeighborPy2 - veloZNeighborMy2) + (eight * four * c1o3) * ( veloZNeighborPy3 - veloZNeighborMy3) - ( veloZNeighborPy4 - veloZNeighborMy4) ) / (seven * ten * four);

	dzvx = ( (twentyeight * eight) * ( veloXNeighborPz - veloXNeighborMz ) - (seven * eight) * ( veloXNeighborPz2 - veloXNeighborMz2) + (eight * four * c1o3) * ( veloXNeighborPz3 - veloXNeighborMz3) - ( veloXNeighborPz4 - veloXNeighborMz4) ) / (seven * ten * four);
	dzvy = ( (twentyeight * eight) * ( veloYNeighborPz - veloYNeighborMz ) - (seven * eight) * ( veloYNeighborPz2 - veloYNeighborMz2) + (eight * four * c1o3) * ( veloYNeighborPz3 - veloYNeighborMz3) - ( veloYNeighborPz4 - veloYNeighborMz4) ) / (seven * ten * four);

    //////////////////////////////////////////////////////////////////////////

	real tmpX = dyvz - dzvy;
	real tmpY = dzvx - dxvz;
	real tmpZ = dxvy - dyvx;
    //////////////////////////////////////////////////////////////////////////

    isFluid[ index ] = 1;

    enstrophy[ index ] = c1o2 * (rho[index] + one) * ( tmpX*tmpX + tmpY*tmpY + tmpZ*tmpZ );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EnstrophyAnalyzer::EnstrophyAnalyzer(SPtr<Parameter> para, uint analyzeIter)
{
	this->para = para;
	this->analyzeIter = analyzeIter;
}

void EnstrophyAnalyzer::writeToFile( std::string filename )
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "EnstrophyAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + "_EnstrophyData.dat" );

    for( auto& EKin : this->enstrophyTimeSeries )
        file << std::setprecision(15) << EKin << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


