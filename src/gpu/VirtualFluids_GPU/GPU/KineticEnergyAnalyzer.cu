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

//#include "Core/Logger/Logger.h"

#include "Parameter/Parameter.h"
// includes, kernels
#include "GPU/GPU_Kernels.cuh"
#include <lbm/constants/NumericConstants.h>

using namespace vf::lbm::constant;
using namespace vf::lbm::dir;

__global__                 void kineticEnergyKernel  (real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, uint size_Mat);

__host__ __device__ inline void kineticEnergyFunction(real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, uint index);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool KineticEnergyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

	int lev = 0;
	int size_Mat = (int)this->para->getParD(lev)->numberOfNodes;

    thrust::device_vector<real> kineticEnergy(size_Mat, c0o1);
    thrust::device_vector<uint> isFluid      (size_Mat, 0);

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

    LBCalcMacCompSP27<<< grid, threads >>> (para->getParD(lev)->velocityX,
											para->getParD(lev)->velocityY,
											para->getParD(lev)->velocityZ,
											para->getParD(lev)->rho,
											para->getParD(lev)->pressure,
											para->getParD(lev)->typeOfGridNode,
											para->getParD(lev)->neighborX,
											para->getParD(lev)->neighborY,
											para->getParD(lev)->neighborZ,
											para->getParD(lev)->numberOfNodes,
											para->getParD(lev)->distributions.f[0],
											para->getParD(lev)->isEvenTimestep); 
    getLastCudaError("LBCalcMacSP27 execution failed"); 

	kineticEnergyKernel <<< grid, threads >>> ( para->getParD(lev)->velocityX, 
											    para->getParD(lev)->velocityY, 
												para->getParD(lev)->velocityZ, 
												para->getParD(lev)->rho, 
											    para->getParD(lev)->neighborX,
											    para->getParD(lev)->neighborY,
											    para->getParD(lev)->neighborZ,
											    para->getParD(lev)->neighborInverse,
											    para->getParD(lev)->typeOfGridNode,
												kineticEnergy.data().get(), 
                                                isFluid.data().get(),
												size_Mat);
	cudaDeviceSynchronize();

	 getLastCudaError("kineticEnergyKernel execution failed");

	 real EKin               = thrust::reduce(kineticEnergy.begin(), kineticEnergy.end(), c0o1, thrust::plus<real>());
     uint numberOfFluidNodes = thrust::reduce(isFluid.begin(),       isFluid.end(),       0,    thrust::plus<uint>());

    //std::cout << "EKin " << EKin << "   " << numberOfFluidNodes << std::endl;

    this->kineticEnergyTimeSeries.push_back( EKin / real(numberOfFluidNodes) );

    //TODO: Should this function probably return nothing?
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void kineticEnergyKernel(real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, uint size_Mat)
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

    //if( index % 34 == 0 || index % 34 == 33 ) return;

    if( index >= size_Mat) return;

	unsigned int BC;
	BC = geo[index];
	if (BC != GEO_FLUID) return;

    kineticEnergyFunction( vx, vy, vz, rho, neighborX, neighborY, neighborZ, neighborWSB, geo, kineticEnergy, isFluid, index );
}

__host__ __device__ void kineticEnergyFunction(real* vx, real* vy, real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ, uint* neighborWSB, uint* geo, real* kineticEnergy, uint* isFluid, uint index)
{
    //////////////////////////////////////////////////////////////////////////////
	//neighbor index                                
	uint k     = index;                             
	uint kPx   = neighborX[k];                      if( geo[ kPx   ] != GEO_FLUID ) return;
	uint kPy   = neighborY[k];                      if( geo[ kPy   ] != GEO_FLUID ) return;
	uint kPz   = neighborZ[k];                      if( geo[ kPz   ] != GEO_FLUID ) return;
	uint kMxyz = neighborWSB[k];                    if( geo[ kMxyz ] != GEO_FLUID ) return;
	uint kMx   = neighborZ[neighborY[kMxyz]];       if( geo[ kMx   ] != GEO_FLUID ) return;
	uint kMy   = neighborZ[neighborX[kMxyz]];       if( geo[ kMy   ] != GEO_FLUID ) return;
	uint kMz   = neighborY[neighborX[kMxyz]];       if( geo[ kMz   ] != GEO_FLUID ) return;
    //////////////////////////////////////////////////////////////////////////

    isFluid[ index ] = 1;

    kineticEnergy[ index ] = c1o2 * ( vx[index] * vx[index] + vy[index] * vy[index] + vz[index] * vz[index] ) * (rho[index] + c1o1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KineticEnergyAnalyzer::KineticEnergyAnalyzer(SPtr<Parameter> para, uint analyzeIter)
{
    this->para = para;
    this->analyzeIter = analyzeIter;
}

void KineticEnergyAnalyzer::writeToFile(std::string filename)
{
    //*logging::out << logging::Logger::INFO_INTERMEDIATE << "KineticEnergyAnalyzer::writeToFile( " << filename << " )" << "\n";
	std::cout << "KineticEnergyAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + "_KineticEnergyData.dat" );

    for( auto& EKin : this->kineticEnergyTimeSeries )
        file << std::setprecision(15) << EKin << std::endl;

    file.close();

    //*logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
	std::cout << "done!\n";
}


