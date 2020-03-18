#include "HeatFluxAnalyzer.h"

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

#include "GksGpu/BoundaryConditions/BoundaryCondition.h"

#include "FlowStateData/AccessDeviceData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"

#include "FluxComputation/SutherlandsLaw.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

__global__                 void heatFluxKernel  ( DataBaseStruct  dataBase, GksGpu::BoundaryConditionStruct  boundaryCondition, Parameters  parameters, real* heatFlux, uint startIndex, uint numberOfEntities );
__host__ __device__ inline void heatFluxFunction( DataBaseStruct& dataBase, GksGpu::BoundaryConditionStruct& boundaryCondition, Parameters& parameters, real* heatFlux, uint startIndex, uint index );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool HeatFluxAnalyzer::run(uint iter, Parameters parameters)
{
    if( iter % this->analyzeIter != 0 ) return false;

    uint numberOfCells = this->boundaryCondition->numberOfCellsPerLevel[ dataBase->numberOfLevels - 1 ];

    thrust::device_vector<real> heatFlux( numberOfCells );

    CudaUtility::CudaGrid grid( numberOfCells, 32 );

    for( uint level = 0; level < dataBase->numberOfLevels - 1; level++ ) parameters.dx *= c1o2; 

    runKernel( heatFluxKernel,
               heatFluxFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               boundaryCondition->toStruct(),
               parameters,
               heatFlux.data().get(),
               boundaryCondition->startOfCellsPerLevel[ dataBase->numberOfLevels - 1 ] );

    getLastCudaError("HeatFluxAnalyzer::run(uint iter)");

    real q = thrust::reduce( heatFlux.begin(), heatFlux.end(), c0o1, thrust::plus<real>() ) * parameters.dx * parameters.dx;

    real qIdeal = c1o4 * (parameters.K + c5o1) * ( parameters.mu / parameters.Pr ) * ( c1o1 / lambdaHot - c1o1 / lambdaCold );

    this->heatFluxTimeSeries.push_back( q / qIdeal );

    if( iter % this->outputIter == 0 ) *logging::out << logging::Logger::INFO_HIGH << "q = " << q / qIdeal << "\n";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void heatFluxKernel(DataBaseStruct  dataBase, GksGpu::BoundaryConditionStruct  boundaryCondition, Parameters  parameters, real* heatFlux, uint startIndex, uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    heatFluxFunction( dataBase, boundaryCondition, parameters, heatFlux, startIndex, index );
}

__host__ __device__ void heatFluxFunction(DataBaseStruct& dataBase, GksGpu::BoundaryConditionStruct& boundaryCondition, Parameters& parameters, real* heatFlux, uint startIndex, uint index)
{
    uint ghostCellIndex  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIndex = boundaryCondition.domainCells[ startIndex + index ];

    if( isCellProperties( dataBase.cellProperties[ domainCellIndex ], CELL_PROPERTIES_GHOST ) )
    {
        heatFlux[ startIndex + index ] = c0o1;
        return;
    }

    //////////////////////////////////////////////////////////////////////////

    ConservedVariables ghostCons;

    readCellData(ghostCellIndex, dataBase, ghostCons);

    ConservedVariables domainCons;

    readCellData(domainCellIndex, dataBase, domainCons);

    PrimitiveVariables ghostPrim  = toPrimitiveVariables(ghostCons,  parameters.K);
    PrimitiveVariables domainPrim = toPrimitiveVariables(domainCons, parameters.K);

    //////////////////////////////////////////////////////////////////////////

    real lambda = c1o2 * (ghostPrim.lambda + domainPrim.lambda);

    real r   = parameters.lambdaRef / lambda;

    real mu;
    if ( parameters.viscosityModel == ViscosityModel::constant ){
        mu = parameters.mu;
    }
    else if ( parameters.viscosityModel == ViscosityModel::sutherlandsLaw ){
        mu = sutherlandsLaw( parameters, r );
    }

    heatFlux[ startIndex + index ] = c1o4 * (parameters.K + c5o1) * ( mu / parameters.Pr ) / parameters.dx * ( c1o1 / domainPrim.lambda - c1o1 / ghostPrim.lambda );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HeatFluxAnalyzer::HeatFluxAnalyzer( SPtr<DataBase> dataBase, SPtr<GksGpu::BoundaryCondition> boundaryCondition, uint analyzeIter, uint outputIter, real lambdaHot, real lambdaCold, real L )
    : dataBase(dataBase), 
      boundaryCondition(boundaryCondition), 
      analyzeIter(analyzeIter), 
      outputIter(outputIter), 
      lambdaHot(lambdaHot), 
      lambdaCold(lambdaCold), 
      L(L)
{
}

void HeatFluxAnalyzer::writeToFile(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "HeatFluxAnalyzer::writeToFile( " << filename << " )" << "\n";

    std::ofstream file;

    file.open(filename + ".dat" );

    for( auto& EKin : this->heatFluxTimeSeries )
        file << std::setprecision(15) << EKin << std::endl;

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}


