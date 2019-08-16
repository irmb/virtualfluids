#include "ConcreteHeatFlux.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#include <thrust/host_vector.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Core/PointerDefinitions.h"
#include "Core/RealConstants.h"
#include "Core/Logger/Logger.h"

#include "DataBase/DataBase.h"
#include "DataBase/DataBaseStruct.h"

#include "Definitions/MemoryAccessPattern.h"
#include "Definitions/PassiveScalar.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"
#include "FlowStateData/ThermalDependencies.cuh"

#include "FluxComputation/Moments.cuh"
#include "FluxComputation/ApplyFlux.cuh"
#include "FluxComputation/Transformation.cuh"
#include "FluxComputation/AssembleFlux.cuh"
#include "FluxComputation/ExpansionCoefficients.cuh"

#include "CudaUtility/CudaRunKernel.hpp"

//////////////////////////////////////////////////////////////////////////

__global__                 void boundaryConditionKernel  ( const DataBaseStruct dataBase, 
                                                           const ConcreteHeatFluxStruct boundaryCondition, 
                                                           const Parameters parameters,
                                                           const uint startIndex,
                                                           const uint numberOfEntities );

__host__ __device__ inline void boundaryConditionFunction( const DataBaseStruct& dataBase, 
                                                           const ConcreteHeatFluxStruct& boundaryCondition, 
                                                           const Parameters& parameters,
                                                           const uint startIndex,
                                                           const uint index );

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void ConcreteHeatFlux::runBoundaryConditionKernel(const SPtr<DataBase> dataBase, 
                                          const Parameters parameters, 
                                          const uint level)
{    
    CudaUtility::CudaGrid grid( this->numberOfCellsPerLevel[ level ], 32 );

    runKernel( boundaryConditionKernel,
               boundaryConditionFunction,
               dataBase->getDeviceType(), grid, 
               dataBase->toStruct(),
               this->toStruct(),
               parameters,
               this->startOfCellsPerLevel[ level ] );

    cudaDeviceSynchronize();

    getLastCudaError("HeatFlux::runBoundaryConditionKernel( const SPtr<DataBase> dataBase, const Parameters parameters, const uint level )");
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

__global__ void boundaryConditionKernel(const DataBaseStruct dataBase, 
                                        const ConcreteHeatFluxStruct boundaryCondition, 
                                        const Parameters parameters,
                                        const uint startIndex,
                                        const uint numberOfEntities)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if( index >= numberOfEntities ) return;

    boundaryConditionFunction( dataBase, boundaryCondition, parameters, startIndex, index );
}

__host__ __device__ inline void boundaryConditionFunction(const DataBaseStruct& dataBase, 
                                                          const ConcreteHeatFluxStruct& boundaryCondition, 
                                                          const Parameters& parameters,
                                                          const uint startIndex,
                                                          const uint index)
{
    uint ghostCellIdx  = boundaryCondition.ghostCells [ startIndex + index ];
    uint domainCellIdx = boundaryCondition.domainCells[ startIndex + index ];

    real dx = boundaryCondition.L / real(boundaryCondition.numberOfPoints + 1);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    PrimitiveVariables domainCellPrim;
    {
        ConservedVariables domainCellData;
        readCellData(domainCellIdx, dataBase, domainCellData);
        domainCellPrim = toPrimitiveVariables(domainCellData, parameters.K);
    }

    real TF = getT(domainCellPrim);

    for( uint i = 0; i < boundaryCondition.numberOfPoints; i++ )
    {
        uint finiteDifferenceIndex = index * boundaryCondition.numberOfPoints + i;

        real T0 = boundaryCondition.temperatures[ finiteDifferenceIndex ];

        real Tn;
        if( i == 0 )
            Tn = TF;
        else
            Tn = boundaryCondition.temperatures[ finiteDifferenceIndex - 1 ];

        real Tp;
        if( i == boundaryCondition.numberOfPoints - 1 )
            Tp = boundaryCondition.ambientTemperature;
        else
            Tp = boundaryCondition.temperatures[ finiteDifferenceIndex + 1 ];

        real dTdxx = ( Tp + Tn - c2o1 * T0 ) / ( dx * dx );

        boundaryCondition.temperatures[ finiteDifferenceIndex ] += parameters.dt * boundaryCondition.temperatureConductivity * dTdxx;
    }

    ConservedVariables flux;

    {
        real T0 = boundaryCondition.temperatures[ index * boundaryCondition.numberOfPoints     ];
        real T1 = boundaryCondition.temperatures[ index * boundaryCondition.numberOfPoints + 1 ];
        real T2 = boundaryCondition.temperatures[ index * boundaryCondition.numberOfPoints + 2 ];

        real k = boundaryCondition.temperatureConductivity * boundaryCondition.density * boundaryCondition.specificHeatCapacity;

        flux.rhoE = - k * ( - c3o2 * T0 + c2o1 * T1 - c1o2 * T2 ) / parameters.dx;
    }

    flux = (parameters.dt * parameters.dx * parameters.dx) * flux;

    applyFluxToNegCell(dataBase, domainCellIdx, flux, 'a', parameters);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

ConcreteHeatFlux::~ConcreteHeatFlux()
{
    checkCudaErrors( cudaFree( this->temperatures ) );
}

ConcreteHeatFlux::ConcreteHeatFlux(SPtr<DataBase> dataBase, uint numberOfPoints, real temperatureConductivity, real density, real specificHeatCapacity, real L, real ambientTemperature)
    : BoundaryCondition( dataBase )
{
    this->numberOfPoints = numberOfPoints;

    this->temperatureConductivity = temperatureConductivity;
    this->density                 = density;
    this->specificHeatCapacity    = specificHeatCapacity;

    this->L = L;
    this->ambientTemperature = ambientTemperature;

    this->temperatures = nullptr;
}

void ConcreteHeatFlux::init()
{
    checkCudaErrors( cudaMalloc( &this->temperatures, sizeof(real) * numberOfPoints * this->numberOfCells ) );

    // initialize values
    thrust::device_ptr<real> dev_ptr(this->temperatures);
    thrust::fill(dev_ptr, dev_ptr + numberOfPoints * this->numberOfCells, this->ambientTemperature);
}

bool ConcreteHeatFlux::isWall()
{
    return true;
}

bool ConcreteHeatFlux::isInsulated()
{
    return true;
}

bool ConcreteHeatFlux::isFluxBC()
{
    return true;
}

bool ConcreteHeatFlux::secondCellsNeeded()
{
    return false;
}

void ConcreteHeatFlux::writeVTKFile(SPtr<DataBase> dataBase, Parameters& parameters, std::string filename)
{
    std::vector<uint> ghostCells (this->numberOfCells);
    std::vector<uint> domainCells(this->numberOfCells);

    std::vector<real> temperatures( this->numberOfCells * this->numberOfPoints );

    checkCudaErrors( cudaMemcpy(ghostCells.data() , this->ghostCells , sizeof(uint) * this->numberOfCells, cudaMemcpyDeviceToHost ) );
    checkCudaErrors( cudaMemcpy(domainCells.data(), this->domainCells, sizeof(uint) * this->numberOfCells, cudaMemcpyDeviceToHost ) );

    checkCudaErrors( cudaMemcpy(temperatures.data(), this->temperatures, sizeof(real) * this->numberOfCells * this->numberOfPoints, cudaMemcpyDeviceToHost ) );


    std::vector<Vec3> nodes;
    std::vector< std::array<uint, 8> > cells;

    for( uint index = 0; index < this->numberOfCells; index ++ )
    {
        real dx = this->L / real(this->numberOfPoints + 1);

        Vec3 displacement = dataBase->nodeCoordinates[ dataBase->cellToNode[ ghostCells [index] ][0] ]
                          - dataBase->nodeCoordinates[ dataBase->cellToNode[ domainCells[index] ][0] ];
        real dn = displacement.length();
        displacement = ( c1o1 / displacement.length() ) * displacement;

        char direction = 'z';
        if ( std::abs(displacement.x) > std::abs(displacement.y) && std::abs(displacement.x) > std::abs(displacement.z) ) direction = 'x';
        if ( std::abs(displacement.y) > std::abs(displacement.x) && std::abs(displacement.y) > std::abs(displacement.z) ) direction = 'y';

        Vec3 dn1, dn2, dn3, dn4;
        if( direction == 'x' )
        {
            dn1.y =  c1o2*dn; dn1.z =  c1o2*dn; 
            dn2.y = -c1o2*dn; dn2.z =  c1o2*dn; 
            dn3.y = -c1o2*dn; dn3.z = -c1o2*dn; 
            dn4.y =  c1o2*dn; dn4.z = -c1o2*dn;
        }
        if( direction == 'y' )
        {
            dn1.x =  c1o2*dn; dn1.z =  c1o2*dn; 
            dn2.x = -c1o2*dn; dn2.z =  c1o2*dn; 
            dn3.x = -c1o2*dn; dn3.z = -c1o2*dn; 
            dn4.x =  c1o2*dn; dn4.z = -c1o2*dn;
        }
        if( direction == 'z' )
        {
            dn1.x =  c1o2*dn; dn1.y =  c1o2*dn; 
            dn2.x = -c1o2*dn; dn2.y =  c1o2*dn; 
            dn3.x = -c1o2*dn; dn3.y = -c1o2*dn; 
            dn4.x =  c1o2*dn; dn4.y = -c1o2*dn;
        }

        //std::cout << "(" << dn1.x << "," << dn1.y << "," << dn1.z << ")" << std::endl;

        Vec3 faceCenter;
        for( uint i = 0; i < 8; i++ )
        {
            faceCenter = faceCenter + c1o16 * dataBase->nodeCoordinates[ dataBase->cellToNode[ ghostCells [index] ][i] ];
            faceCenter = faceCenter + c1o16 * dataBase->nodeCoordinates[ dataBase->cellToNode[ domainCells[index] ][i] ];
        }

        uint nodeStartNumber = nodes.size();

        nodes.push_back( faceCenter + dn1 );
        nodes.push_back( faceCenter + dn2 );
        nodes.push_back( faceCenter + dn3 );
        nodes.push_back( faceCenter + dn4 );

        for( uint i = 1; i <= this->numberOfPoints; i++ )
        {
            nodes.push_back( faceCenter + real(i) * dx * displacement + dn1 );
            nodes.push_back( faceCenter + real(i) * dx * displacement + dn2 );
            nodes.push_back( faceCenter + real(i) * dx * displacement + dn3 );
            nodes.push_back( faceCenter + real(i) * dx * displacement + dn4 );
        }

        nodes.push_back( faceCenter + this->L * displacement + dn1 );
        nodes.push_back( faceCenter + this->L * displacement + dn2 );
        nodes.push_back( faceCenter + this->L * displacement + dn3 );
        nodes.push_back( faceCenter + this->L * displacement + dn4 );

        for( uint i = 0; i <= this->numberOfPoints; i++ )
        {
            cells.push_back({});
            cells.back()[0] = nodeStartNumber + (i    ) * 4    ;
            cells.back()[1] = nodeStartNumber + (i    ) * 4 + 1;
            cells.back()[2] = nodeStartNumber + (i    ) * 4 + 2;
            cells.back()[3] = nodeStartNumber + (i    ) * 4 + 3;
            cells.back()[4] = nodeStartNumber + (i + 1) * 4    ;
            cells.back()[5] = nodeStartNumber + (i + 1) * 4 + 1;
            cells.back()[6] = nodeStartNumber + (i + 1) * 4 + 2;
            cells.back()[7] = nodeStartNumber + (i + 1) * 4 + 3;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "ConcreteHeatFlux::writeVTK( " << filename << " )" << "\n";
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "    nCells = " << this->numberOfCells << " )" << "\n";


    std::ofstream file;

    file.open(filename + ".vtk");

    file << "# vtk DataFile Version 3.0\n";
    file << "by MeshGenerator\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << nodes.size() << " float" << std::endl;

    for (auto node : nodes){
        file << node.x << " " << node.y << " " << node.z << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELLS " << cells.size() << " " << cells.size() * 9 << std::endl;

    for ( auto cell : cells ){

        file << 8 << " ";

        for( auto node : cell ) file << node << " ";

        file << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "CELL_TYPES " << cells.size() << std::endl;

    for ( uint cellIdx = 0; cellIdx < cells.size(); cellIdx++ ){
        file << 12 << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file << "\nPOINT_DATA " << nodes.size() << std::endl;

    file << "FIELD Label " << 1 << std::endl;

    //////////////////////////////////////////////////////////////////////////

    file << "T 1 " << nodes.size() << " double" << std::endl;

    for ( uint index = 0; index < this->numberOfCells; index++ )
    {
        ConservedVariables cons;

        cons.rho  = dataBase->dataHost[ RHO__(domainCells[ index ], dataBase->numberOfCells) ];
        cons.rhoU = dataBase->dataHost[ RHO_U(domainCells[ index ], dataBase->numberOfCells) ];
        cons.rhoV = dataBase->dataHost[ RHO_V(domainCells[ index ], dataBase->numberOfCells) ];
        cons.rhoW = dataBase->dataHost[ RHO_W(domainCells[ index ], dataBase->numberOfCells) ];
        cons.rhoE = dataBase->dataHost[ RHO_E(domainCells[ index ], dataBase->numberOfCells) ];
#ifdef USE_PASSIVE_SCALAR
        cons.rhoS_1 = dataBase->dataHost[ RHO_S_1(domainCells[ index ], dataBase->numberOfCells) ];
        cons.rhoS_2 = dataBase->dataHost[ RHO_S_2(domainCells[ index ], dataBase->numberOfCells) ];
#endif // USE_PASSIVE_SCALAR

        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);
        
#ifdef USE_PASSIVE_SCALAR
        real T = getT(prim);
#else // USE_PASSIVE_SCALAR
        real T = 1.0 / prim.lambda;
#endif // USE_PASSIVE_SCALAR

        file << T << std::endl;
        file << T << std::endl;
        file << T << std::endl;
        file << T << std::endl;

        for( uint i = 0; i < this->numberOfPoints; i++ )
        {
            file << temperatures[ numberOfPoints * index + i ] << std::endl;
            file << temperatures[ numberOfPoints * index + i ] << std::endl;
            file << temperatures[ numberOfPoints * index + i ] << std::endl;
            file << temperatures[ numberOfPoints * index + i ] << std::endl;
        }

        file << this->ambientTemperature << std::endl;
        file << this->ambientTemperature << std::endl;
        file << this->ambientTemperature << std::endl;
        file << this->ambientTemperature << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!" << "\n";
}

