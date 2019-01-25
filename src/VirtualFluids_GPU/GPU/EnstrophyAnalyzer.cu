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
	//getVeloLocal//
	real veloLocalX = veloX[k];
	real veloLocalY = veloY[k];
	real veloLocalZ = veloZ[k];
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

	dxvy = (veloYNeighborPx - veloYNeighborMx) / two;
	dxvz = (veloZNeighborPx - veloZNeighborMx) / two;

	dyvx = (veloXNeighborPy - veloXNeighborMy) / two;
	dyvz = (veloZNeighborPy - veloZNeighborMy) / two;

	dzvx = (veloXNeighborPz - veloXNeighborMz) / two;
	dzvy = (veloYNeighborPz - veloYNeighborMz) / two;

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
        file << EKin << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


