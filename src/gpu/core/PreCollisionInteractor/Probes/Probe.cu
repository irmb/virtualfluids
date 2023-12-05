//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file Probe.h
//! \author Henry Korb, Henrik Asmuth
//=======================================================================================

#include "Probe.h"

#include <cmath>
#include <cuda_runtime.h>
#include <cuda.h>
#include <filesystem>
#include <helper_cuda.h>
#include <memory>

#include <basics/constants/NumericConstants.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <basics/StringUtilities/StringUtil.h>

#include "gpu/core/DataStructureInitializer/GridProvider.h"
#include "gpu/core/Cuda/CudaMemoryManager.h"
#include "gpu/core/Utilities/GeometryUtils.h"
#include "gpu/core/LBM/GPUHelperFunctions/KernelUtilities.h"
#include "gpu/core/Output/FilePartCalculator.h"
#include "gpu/core/Parameter/Parameter.h"

using namespace vf::basics::constant;

__host__ __device__ int calcArrayIndex(int node, int nNodes, int timestep, int nTimesteps, int array)
{
    return node+nNodes*(timestep+nTimesteps*array);
}

uint calcOldTimestep(uint currentTimestep, uint lastTimestepInOldSeries)
{
    return currentTimestep > 0 ? currentTimestep - 1 : lastTimestepInOldSeries; 
}

__device__ void calculatePointwiseQuantities(uint oldTimestepInTimeseries, uint timestepInTimeseries, uint timestepInAverage,
                                             uint nTimesteps, real* quantityArray, bool* quantities,
                                             uint* quantityArrayOffsets, uint nPoints, uint node, real vx, real vy, real vz,
                                             real rho)
{
    //"https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm"
    // also has extensions for higher order and covariances
    int n = timestepInAverage+1;
    real inv_n = 1/real(n);

    if (quantities[int(Statistic::Instantaneous)])
    {
        uint arrOff = quantityArrayOffsets[int(Statistic::Instantaneous)];
        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 0)] = vx;
        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 1)] = vy;
        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 2)] = vz;
        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 3)] = rho;
    }

    if (quantities[int(Statistic::Means)])
    {
        uint arrOff = quantityArrayOffsets[int(Statistic::Means)];
        real vx_m_old  = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 0)];
        real vy_m_old  = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 1)];
        real vz_m_old  = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 2)];
        real rho_m_old = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 3)];

        real vx_m_new  = ((n - 1) * vx_m_old  + vx)  * inv_n;
        real vy_m_new  = ((n - 1) * vy_m_old  + vy)  * inv_n;
        real vz_m_new  = ((n - 1) * vz_m_old  + vz)  * inv_n;
        real rho_m_new = ((n - 1) * rho_m_old + rho) * inv_n;

        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 0)] = vx_m_new;
        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 1)] = vy_m_new;
        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 2)] = vz_m_new;
        quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 3)] = rho_m_new;

        if (quantities[int(Statistic::Variances)])
        {
            arrOff = quantityArrayOffsets[int(Statistic::Variances)];

            real vx_var_old  = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 0)];
            real vy_var_old  = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 1)];
            real vz_var_old  = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 2)];
            real rho_var_old = quantityArray[calcArrayIndex(node, nPoints, oldTimestepInTimeseries, nTimesteps, arrOff + 3)];

            real vx_var_new  = ((n - 1) * (vx_var_old)  + (vx  - vx_m_old)  * (vx  - vx_m_new))  * inv_n;
            real vy_var_new  = ((n - 1) * (vy_var_old)  + (vy  - vy_m_old)  * (vy  - vy_m_new))  * inv_n;
            real vz_var_new  = ((n - 1) * (vz_var_old)  + (vz  - vz_m_old)  * (vz  - vz_m_new))  * inv_n;
            real rho_var_new = ((n - 1) * (rho_var_old) + (rho - rho_m_old) * (rho - rho_m_new)) * inv_n;

            quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 0)] = vx_var_new;
            quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 1)] = vy_var_new;
            quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 2)] = vz_var_new;
            quantityArray[calcArrayIndex(node, nPoints, timestepInTimeseries, nTimesteps, arrOff + 3)] = rho_var_new;
        }
    }
}

__global__ void calcQuantitiesKernel(uint* pointIndices, uint nPoints, uint oldTimestepInTimeseries,
                                     uint timestepInTimeseries, uint timestepInAverage, uint nTimesteps, real* vx, real* vy,
                                     real* vz, real* rho, uint* neighborX, uint* neighborY, uint* neighborZ,
                                     bool* quantities, uint* quantityArrayOffsets, real* quantityArray)
{
    const uint node = vf::gpu::getNodeIndex();

    if(node>=nPoints) return;

    // Get indices of neighbor nodes. 
    // node referring to BSW cell as seen from probe point
    uint k = pointIndices[node];
    real u_interpX, u_interpY, u_interpZ, rho_interp;

    u_interpX = vx[k];
    u_interpY = vy[k];
    u_interpZ = vz[k];
    rho_interp = rho[k];

    calculatePointwiseQuantities(oldTimestepInTimeseries, timestepInTimeseries, timestepInAverage, nTimesteps, quantityArray, quantities, quantityArrayOffsets, nPoints, node, u_interpX, u_interpY, u_interpZ, rho_interp);
}

__global__ void interpAndCalcQuantitiesKernel(uint* pointIndices, uint nPoints, uint oldTimestepInTimeseries,
                                              uint timestepInTimeseries, uint timestepInAverage, uint nTimesteps,
                                              real* distX, real* distY, real* distZ, real* vx, real* vy, real* vz, real* rho,
                                              uint* neighborX, uint* neighborY, uint* neighborZ, bool* quantities,
                                              uint* quantityArrayOffsets, real* quantityArray)
{
    const uint node = vf::gpu::getNodeIndex();

    if(node>=nPoints) return;

    // Get indices of neighbor nodes. 
    // node referring to BSW cell as seen from probe point
    uint k = pointIndices[node];
    real u_interpX, u_interpY, u_interpZ, rho_interp;

    uint ke, kn, kt, kne, kte, ktn, ktne;
    getNeighborIndicesOfBSW(  k, ke, kn, kt, kne, kte, ktn, ktne, neighborX, neighborY, neighborZ);

    // Trilinear interpolation of macroscopic quantities to probe point
    real dW, dE, dN, dS, dT, dB;
    getInterpolationWeights(dW, dE, dN, dS, dT, dB, distX[node], distY[node], distZ[node]);

    u_interpX  = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vx );
    u_interpY  = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vy );
    u_interpZ  = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vz );
    rho_interp = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, rho );

    calculatePointwiseQuantities(oldTimestepInTimeseries, timestepInTimeseries, timestepInAverage, nTimesteps, quantityArray, quantities, quantityArrayOffsets, nPoints, node, u_interpX, u_interpY, u_interpZ, rho_interp);

}

bool Probe::getHasDeviceQuantityArray(){ return this->hasDeviceQuantityArray; }

real Probe::getNondimensionalConversionFactor(int level){ return c1o1; }

void Probe::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaMemoryManager)
{
    using std::placeholders::_1;
    this->velocityRatio      = std::bind(&Parameter::getScaledVelocityRatio,        para, _1); 
    this->densityRatio       = std::bind(&Parameter::getScaledDensityRatio,         para, _1);
    this->forceRatio         = std::bind(&Parameter::getScaledForceRatio,           para, _1);
    this->stressRatio        = std::bind(&Parameter::getScaledStressRatio,          para, _1);
    this->viscosityRatio     = std::bind(&Parameter::getScaledViscosityRatio,       para, _1);
    this->nondimensional     = std::bind(&Probe::getNondimensionalConversionFactor, this, _1);

    probeParams.resize(para->getMaxLevel()+1);

    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        std::vector<int> probeIndices_level;
        std::vector<real> distX_level;
        std::vector<real> distY_level;
        std::vector<real> distZ_level;        
        std::vector<real> pointCoordsX_level;
        std::vector<real> pointCoordsY_level;
        std::vector<real> pointCoordsZ_level;
        
        this->findPoints(para, gridProvider, probeIndices_level, distX_level, distY_level, distZ_level,      
                       pointCoordsX_level, pointCoordsY_level, pointCoordsZ_level,
                       level);
        
        this->addProbeStruct(para, cudaMemoryManager, probeIndices_level, 
                            distX_level, distY_level, distZ_level, 
                            pointCoordsX_level, pointCoordsY_level, pointCoordsZ_level, 
                            level);

        if(this->outputTimeSeries) timeseriesFileNames.push_back(this->writeTimeseriesHeader(para, level));
    }
}

void Probe::addProbeStruct( Parameter* para, CudaMemoryManager* cudaMemoryManager, std::vector<int>& probeIndices,
                            std::vector<real>& distX, std::vector<real>& distY, std::vector<real>& distZ,   
                            std::vector<real>& pointCoordsX, std::vector<real>& pointCoordsY, std::vector<real>& pointCoordsZ,
                            int level)
{
    probeParams[level] = std::make_shared<ProbeStruct>();
    probeParams[level]->nTimesteps = this->getNumberOfTimestepsInTimeseries(para, level);
    probeParams[level]->nPoints  = uint(pointCoordsX.size()); // Note, need to have both nPoints and nIndices because they differ in PlanarAverage
    probeParams[level]->nIndices = uint(probeIndices.size());

    probeParams[level]->pointCoordsX = (real*)malloc(probeParams[level]->nPoints * sizeof(real));
    probeParams[level]->pointCoordsY = (real*)malloc(probeParams[level]->nPoints * sizeof(real));
    probeParams[level]->pointCoordsZ = (real*)malloc(probeParams[level]->nPoints * sizeof(real));

    std::copy(pointCoordsX.begin(), pointCoordsX.end(), probeParams[level]->pointCoordsX);
    std::copy(pointCoordsY.begin(), pointCoordsY.end(), probeParams[level]->pointCoordsY);
    std::copy(pointCoordsZ.begin(), pointCoordsZ.end(), probeParams[level]->pointCoordsZ);

    // Note, dist only needed for kernels that do interpolate
    if( !distX.empty() && !distY.empty() && !distZ.empty() )
    {
        probeParams[level]->hasDistances=true;
        cudaMemoryManager->cudaAllocProbeDistances(this, level);
        std::copy(distX.begin(), distX.end(), probeParams[level]->distXH);
        std::copy(distY.begin(), distY.end(), probeParams[level]->distYH);
        std::copy(distZ.begin(), distZ.end(), probeParams[level]->distZH);
        cudaMemoryManager->cudaCopyProbeDistancesHtoD(this, level);
    }  
    
    cudaMemoryManager->cudaAllocProbeIndices(this, level);
    std::copy(probeIndices.begin(), probeIndices.end(), probeParams[level]->pointIndicesH);
    cudaMemoryManager->cudaCopyProbeIndicesHtoD(this, level);

    uint arrOffset = 0;

    cudaMemoryManager->cudaAllocProbeQuantitiesAndOffsets(this, level);

    for( int var=0; var<int(Statistic::LAST); var++)
    {
        if(this->quantities[var])
        {
            probeParams[level]->quantitiesH[var] = true;
            probeParams[level]->arrayOffsetsH[var] = arrOffset;
            arrOffset += uint( this->getPostProcessingVariables(static_cast<Statistic>(var)).size() ); 
        }
    }
    
    cudaMemoryManager->cudaCopyProbeQuantitiesAndOffsetsHtoD(this, level);

    probeParams[level]->nArrays = arrOffset;

    cudaMemoryManager->cudaAllocProbeQuantityArray(this, level);

    std::fill_n(probeParams[level]->quantitiesArrayH, probeParams[level]->nArrays*probeParams[level]->nPoints*probeParams[level]->nTimesteps, c0o1);

    if(this->hasDeviceQuantityArray)
        cudaMemoryManager->cudaCopyProbeQuantityArrayHtoD(this, level);

}

void Probe::interact(Parameter* para, CudaMemoryManager* cudaMemoryManager, int level, uint t)
{
    uint t_level = para->getTimeStep(level, t, false);

    SPtr<ProbeStruct> probeStruct = this->getProbeStruct(level);

    //!Skip empty probes
    if(probeStruct->nPoints==0) return;

    //! if tAvg==1 the probe will be evaluated in every sub-timestep of each respective level
    //! else, the probe will only be evaluated in each synchronous time step tAvg

    uint level_coefficient = exp2(level);

    uint tAvg_level = this->tAvg == 1 ? this->tAvg : this->tAvg * level_coefficient;
    uint tOut_level = this->tOut * level_coefficient;
    uint tStartOut_level = this->tStartOut * level_coefficient;
    uint tStartAvg_level = this->tStartAvg * level_coefficient;

    uint tAfterStartAvg = t_level - tStartAvg_level;
    uint tAfterStartOut = t_level - tStartOut_level;

    if( (t > this->tStartAvg) && (tAfterStartAvg % tAvg_level == 0))
    {
        this->calculateQuantities(probeStruct, para, t_level, level);

        if(t > this->tStartTmpAveraging) probeStruct->timestepInTimeAverage++;
        if(this->outputTimeSeries && (t_level >= tStartOut_level)) probeStruct->timestepInTimeseries++;
    }

    //! output only in synchronous timesteps
    if( (t > this->tStartOut) && (tAfterStartOut % tOut_level == 0) )
    {   
        if(this->hasDeviceQuantityArray)
            cudaMemoryManager->cudaCopyProbeQuantityArrayDtoH(this, level);
        this->write(para, level, t);
        
        if(level == 0&& !this->outputTimeSeries) this->writeParallelFile(para, t);

        if(this->outputTimeSeries)
        {
            probeStruct->lastTimestepInOldTimeseries = probeStruct->timestepInTimeseries > 0 ? probeStruct->timestepInTimeseries - 1: 0;
            probeStruct->timestepInTimeseries = 0;
        }
    }
}

void Probe::free(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    for(int level=0; level<=para->getMaxLevel(); level++)
    {   
        if(this->probeParams[level]->hasDistances)
            cudaMemoryManager->cudaFreeProbeDistances(this, level);
        cudaMemoryManager->cudaFreeProbeIndices(this, level);
        cudaMemoryManager->cudaFreeProbeQuantityArray(this, level);
        cudaMemoryManager->cudaFreeProbeQuantitiesAndOffsets(this, level);
    }
}


void Probe::addStatistic(Statistic variable)
{
    if (!this->isAvailableStatistic(variable)) throw std::runtime_error("Probe::addStatistic(): Statistic not available for this probe type!");

    this->quantities[int(variable)] = true;
    switch(variable)
    {
        case Statistic::Variances: 
            this->addStatistic(Statistic::Means); break;

        default: break;
    }
}

template<typename T>
std::string nameComponent(std::string name, T value)
{
    return "_" + name + "_" + StringUtil::toString<T>(value); 
}

std::string Probe::makeParallelFileName(int id, int t)
{
    return this->probeName + "_bin" + nameComponent<int>("ID", id) + nameComponent<int>("t", t) + ".vtk"; 

}

std::string Probe::makeGridFileName(int level, int id, int t, uint part)
{
    return this->probeName + "_bin" + nameComponent<int>("lev", level)
                                    + nameComponent<int>("ID", id)
                                    + nameComponent<int>("Part", part)
                                    + nameComponent<int>("t", t) + ".vtk";
}

std::string Probe::makeTimeseriesFileName(int level, int id)
{
    return this->probeName + "_timeseries" + nameComponent<int>("lev", level)
                                    + nameComponent<int>("ID", id)
                                    + ".txt";
}
void Probe::addAllAvailableStatistics()
{
    for( int var=0; var < int(Statistic::LAST); var++)
    {
        if(this->isAvailableStatistic(static_cast<Statistic>(var))) 
            this->addStatistic(static_cast<Statistic>(var));
    }
}

void Probe::write(Parameter* para, int level, int t)
{
    if(this->outputTimeSeries)
    {
        this->appendTimeseriesFile(para, level, t);
    }
    else
    {
        int t_write = this->fileNameLU ? t: t/this->tOut; 

        const uint numberOfParts = this->getProbeStruct(level)->nPoints / FilePartCalculator::limitOfNodesForVTK + 1;

        std::vector<std::string> fnames;
        for (uint i = 1; i <= numberOfParts; i++)
        {
            this->writeGridFile(para, level, t_write, i);
        }
    }

}

void Probe::writeParallelFile(Parameter* para, int t)
{
    int t_write = this->fileNameLU ? t: t/this->tOut; 
    std::string filename = this->outputPath + this->makeParallelFileName(para->getMyProcessID(), t_write);

    std::vector<std::string> nodedatanames = this->getVarNames();
    std::vector<std::string> cellNames;

    getWriter()->writeParallelFile(filename, fileNamesForCollectionFile, nodedatanames, cellNames);

    this->fileNamesForCollectionFile.clear();
}

void Probe::writeGridFile(Parameter* para, int level, int t, uint part)
{
    std::string fname = this->outputPath + this->makeGridFileName(level, para->getMyProcessID(), t, part);

    std::vector< UbTupleFloat3 > nodes;
    std::vector< std::string > nodedatanames = this->getVarNames();

    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    SPtr<ProbeStruct> probeStruct = this->getProbeStruct(level);

    uint startpos = (part-1) * FilePartCalculator::limitOfNodesForVTK;
    uint sizeOfNodes = std::min(FilePartCalculator::limitOfNodesForVTK, probeStruct->nPoints - startpos);
    uint endpos = startpos + sizeOfNodes;

    //////////////////////////////////////////////////////////////////////////
    nodes.resize(sizeOfNodes);

    for (uint pos = startpos; pos < endpos; pos++)
    {
        nodes[pos-startpos] = makeUbTuple(  float(probeStruct->pointCoordsX[pos]),
                                            float(probeStruct->pointCoordsY[pos]),
                                            float(probeStruct->pointCoordsZ[pos]));
    }

    for( auto it=nodedata.begin(); it!=nodedata.end(); it++) it->resize(sizeOfNodes);

    uint arrLen = probeStruct->nPoints;
    int nTimesteps = probeStruct->nTimesteps;
    int timestep = probeStruct->timestepInTimeseries;

    for( int var=0; var < int(Statistic::LAST); var++)
    {           
        if(this->quantities[var])
        {

            Statistic statistic = static_cast<Statistic>(var);
            std::vector<PostProcessingVariable> postProcessingVariables = this->getPostProcessingVariables(statistic);
            uint n_arrs = uint(postProcessingVariables.size());

            uint arrOff = probeStruct->arrayOffsetsH[var];

            for(uint arr=0; arr<n_arrs; arr++)
            {
                real coeff = postProcessingVariables[arr].conversionFactor(level);
                
                for (uint pos = startpos; pos < endpos; pos++)
                {
                    nodedata[arrOff+arr][pos-startpos] = double(probeStruct->quantitiesArrayH[calcArrayIndex(pos, arrLen, timestep, nTimesteps, arrOff+arr)]*coeff);
                }
            }
        }
    }
    std::string fullName = getWriter()->writeNodesWithNodeData(fname, nodes, nodedatanames, nodedata);
    this->fileNamesForCollectionFile.push_back(fullName.substr(fullName.find_last_of('/') + 1));
}

std::string Probe::writeTimeseriesHeader(Parameter* para, int level)
{
/*
File Layout:
TimeseriesOutput
Quantities: Quant1 Quant2 Quant3
Positions:
point1.x, point1.y, point1.z
point2.x, point2.y, point2.z
...
t0 point1.quant1 point2.quant1 ... point1.quant2 point2.quant2 ...
t1 point1.quant1 point2.quant1 ... point1.quant2 point2.quant2 ...
*/
    auto probeStruct = this->getProbeStruct(level);
    std::filesystem::create_directories(this->outputPath);
    std::string fname = this->outputPath + this->makeTimeseriesFileName(level, para->getMyProcessID());
    std::ofstream out(fname.c_str(), std::ios::out | std::ios::binary);

    if(!out.is_open()) throw std::runtime_error("Could not open timeseries file " + fname + "!");

    out << "TimeseriesOutput \n";
    out << "Quantities: ";
    for(const std::string& name : getVarNames())
        out << name << ", ";
    out << "\n";
    out << "Number of points in this file: \n";
    out << probeStruct->nPoints << "\n";
    out << "Positions: x, y, z\n";
    for( uint i=0; i<probeStruct->nPoints; i++)
        out << probeStruct->pointCoordsX[i] << ", " << probeStruct->pointCoordsY[i] << ", " << probeStruct->pointCoordsZ[i] << "\n";

    out.close();

    return fname;
}

void Probe::appendTimeseriesFile(Parameter* para, int level, int t)
{
    std::ofstream out(this->timeseriesFileNames[level], std::ios::app | std::ios::binary);

    // uint t_level = para->getTimeStep(level, t, false);
    uint tAvg_level = this->tAvg==1 ? this->tAvg: this->tAvg*exp2(-level);

    real dt = para->getTimeRatio()*tAvg_level;
    auto probeStruct = this->getProbeStruct(level);

    real t_start = ( t-this->tOut )*para->getTimeRatio();

    int vals_per_timestep = probeStruct->nPoints*probeStruct->nArrays+1;

    real* timestep_array = (real*) malloc(vals_per_timestep*sizeof(real));

    for(uint timestep=0; timestep<probeStruct->timestepInTimeseries; timestep++)
    {
        int val = 0;
        timestep_array[val] = t_start+timestep*dt;
        val++;

        for( int var=0; var < int(Statistic::LAST); var++)
        {           
            if(!this->quantities[var]) continue;
            
            Statistic statistic = static_cast<Statistic>(var);
            std::vector<PostProcessingVariable> postProcessingVariables = this->getPostProcessingVariables(statistic);
            uint n_arrs = uint(postProcessingVariables.size());

            uint arrOff = probeStruct->arrayOffsetsH[var];

            for(uint arr=0; arr<n_arrs; arr++)
            {
                real coeff = postProcessingVariables[arr].conversionFactor(level);
                for(uint point=0; point<probeStruct->nPoints; point++)
                {
                    timestep_array[val] = probeStruct->quantitiesArrayH[calcArrayIndex(point, probeStruct->nPoints, timestep, probeStruct->nTimesteps, arrOff+arr)]*coeff;
                    val++;
                }
            }
            
        }
        out.write((char*) timestep_array, sizeof(real)*vals_per_timestep);
    }
    out.close();
}



std::vector<std::string> Probe::getVarNames()
{
    std::vector<std::string> varNames;
    for( int statistic=0; statistic < int(Statistic::LAST); statistic++)
    {
        if(this->quantities[statistic])
        {
            std::vector<PostProcessingVariable> postProcessingVariables = this->getPostProcessingVariables(static_cast<Statistic>(statistic));            
            for(size_t i = 0; i<postProcessingVariables.size(); i++) 
                varNames.push_back(postProcessingVariables[i].name);
        }
    }
    return varNames;
}
