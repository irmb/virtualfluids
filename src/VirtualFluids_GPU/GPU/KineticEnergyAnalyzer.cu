#include "KineticEnergyAnalyzer.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cmath>
#include <sstream>

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>

#include <iomanip>

#include "Core/Logger/Logger.h"

#include "Parameter\Parameter.h"
// includes, kernels
#include "GPU/GPU_Kernels.cuh"
#include "GPU/constant.h"

__global__                 void kineticEnergyKernel  (real* vx, real* vy, real* vz, real* rho, real* kineticEnergy, uint size_Mat);

__host__ __device__ inline void kineticEnergyFunction(real* vx, real* vy, real* vz, real* rho, real* kineticEnergy, uint index);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool KineticEnergyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

	int lev = 0;
	int size_Mat = this->para->getParD(lev)->size_Mat_SP;

    //thrust::device_vector<real> kineticEnergy(size_Mat, zero);
	thrust::device_vector<real> kineticEnergy(size_Mat);

	//printf("%d \n", size_Mat);
	//printf("%d \n", kineticEnergy.begin());
	//printf("%d \n\n", kineticEnergy.end());

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
     getLastCudaError("LBCalcMacSP27 execution failed"); 

	 kineticEnergyKernel <<< grid, threads >>> (para->getParD(lev)->vx_SP, 
												para->getParD(lev)->vy_SP, 
												para->getParD(lev)->vz_SP, 
												para->getParD(lev)->rho_SP, 
												kineticEnergy.data().get(), 
												size_Mat);
	 cudaDeviceSynchronize();

	 getLastCudaError("kineticEnergyKernel execution failed");


	 //thrust::host_vector<real> test = kineticEnergy;

	 //for (auto val : test) printf("%f\n", val);

	 real EKin = thrust::reduce(kineticEnergy.begin(), kineticEnergy.end(), zero, thrust::plus<real>());// / real(size_Mat);

    this->kineticEnergyTimeSeries.push_back( EKin );

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void kineticEnergyKernel(real* vx, real* vy, real* vz, real* rho, real* kineticEnergy, uint size_Mat)
{
	////////////////////////////////////////////////////////////////////////////////
	uint  index;                   // Zugriff auf arrays im device
									   //
	uint tx = threadIdx.x;     // Thread index = lokaler i index
	uint by = blockIdx.x;      // Block index x
	uint bz = blockIdx.y;      // Block index y
	uint  x = tx + STARTOFFX;  // Globaler x-Index 
	uint  y = by + STARTOFFY;  // Globaler y-Index 
	uint  z = bz + STARTOFFZ;  // Globaler z-Index 

	const unsigned sizeX = blockDim.x;
	const unsigned sizeY = gridDim.x;
	const unsigned nx = sizeX + 2 * STARTOFFX;
	const unsigned ny = sizeY + 2 * STARTOFFY;

	index = nx*(ny*z + y) + x;
	////////////////////////////////////////////////////////////////////////////////

    if( index >= size_Mat) return;

    kineticEnergyFunction( vx, vy, vz, rho, kineticEnergy, index );
}

__host__ __device__ void kineticEnergyFunction(real* vx, real* vy, real* vz, real* rho, real* kineticEnergy, uint index)
{
    kineticEnergy[ index ] = c1o2 * ( vx[index] * vx[index] + vy[index] * vy[index] + vz[index] * vz[index] ) * (rho[index] + one);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KineticEnergyAnalyzer::KineticEnergyAnalyzer(SPtr<Parameter> para, uint analyzeIter)
{
    this->para = para;
    this->analyzeIter = analyzeIter;
}

void KineticEnergyAnalyzer::writeToFile(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "KineticEnergyAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + "_KineticEnergyData.dat" );

    for( auto& EKin : this->kineticEnergyTimeSeries )
        file << EKin << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


