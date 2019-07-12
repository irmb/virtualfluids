#include "KineticEnergyAnalyzer.h"

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

__global__                 void kineticEnergyKernel  ( DataBaseStruct dataBase, real* kineticEnergy, uint startIndex, uint numberOfEntities );

__host__ __device__ inline void kineticEnergyFunction( DataBaseStruct dataBase, real* kineticEnergy, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool KineticEnergyAnalyzer::run(uint iter)
{
    if( iter % this->analyzeIter != 0 ) return false;

    thrust::device_vector<real> kineticEnergy( this->dataBase->perLevelCount[ 0 ].numberOfBulkCells );

    CudaUtility::CudaGrid grid( dataBase->perLevelCount[ 0 ].numberOfBulkCells, 32 );

    runKernel( kineticEnergyKernel,
               kineticEnergyFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               kineticEnergy.data().get(),
               dataBase->perLevelCount[ 0 ].startOfCells );

    getLastCudaError("KineticEnergyAnalyzer::run(uint iter)");

    real EKin = thrust::reduce( kineticEnergy.begin(), kineticEnergy.end(), c0o1, thrust::plus<real>() )
              / real(dataBase->perLevelCount[ 0 ].numberOfBulkCells);

    this->kineticEnergyTimeSeries.push_back( EKin );

    //*logging::out << logging::Logger::INFO_HIGH << "EKin = " << EKin << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void kineticEnergyKernel(DataBaseStruct dataBase, real* kineticEnergy, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    kineticEnergyFunction( dataBase, kineticEnergy, startIndex, index );
}

__host__ __device__ void kineticEnergyFunction(DataBaseStruct dataBase, real* kineticEnergy, uint startIndex, uint index)
{
    uint cellIndex = startIndex + index;

    //////////////////////////////////////////////////////////////////////////

    ConservedVariables cons;

    cons.rho  = dataBase.data[ RHO__(cellIndex, dataBase.numberOfCells) ];
    cons.rhoU = dataBase.data[ RHO_U(cellIndex, dataBase.numberOfCells) ];
    cons.rhoV = dataBase.data[ RHO_V(cellIndex, dataBase.numberOfCells) ];
    cons.rhoW = dataBase.data[ RHO_W(cellIndex, dataBase.numberOfCells) ];

    //////////////////////////////////////////////////////////////////////////

    kineticEnergy[ cellIndex ] = c1o2 * ( cons.rhoU * cons.rhoU + cons.rhoV * cons.rhoV + cons.rhoW * cons.rhoW ) / cons.rho;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KineticEnergyAnalyzer::KineticEnergyAnalyzer(SPtr<DataBase> dataBase, uint analyzeIter, uint outputIter)
{
    this->dataBase = dataBase;

    this->analyzeIter = analyzeIter;
    this->outputIter  = outputIter;
}

void KineticEnergyAnalyzer::writeToFile(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "KineticEnergyAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + ".dat" );

    for( auto& EKin : this->kineticEnergyTimeSeries )
        file << std::setprecision(15) << EKin << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


