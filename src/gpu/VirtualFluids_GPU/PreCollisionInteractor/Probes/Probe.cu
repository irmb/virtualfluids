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

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "VirtualFluids_GPU/GPU/GeometryUtils.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include <Core/StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"


__device__ void calculatePointwiseQuantities(uint n, real* quantityArray, bool* quantities, uint* quantityArrayOffsets, uint nPoints, uint node, real vx, real vy, real vz, real rho)
{
    //"https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm"
    // also has extensions for higher order and covariances
    real inv_n = 1/real(n);

    if(quantities[int(Statistic::Instantaneous)])
    {
        uint arrOff = quantityArrayOffsets[int(Statistic::Instantaneous)];
        quantityArray[(arrOff+0)*nPoints+node] = vx;
        quantityArray[(arrOff+1)*nPoints+node] = vy;
        quantityArray[(arrOff+2)*nPoints+node] = vz;
        quantityArray[(arrOff+3)*nPoints+node] = rho;
    }

    if(quantities[int(Statistic::Means)])
    {
        
        uint arrOff = quantityArrayOffsets[int(Statistic::Means)];
        real vx_m_old  = quantityArray[(arrOff+0)*nPoints+node];
        real vy_m_old  = quantityArray[(arrOff+1)*nPoints+node];
        real vz_m_old  = quantityArray[(arrOff+2)*nPoints+node];
        real rho_m_old = quantityArray[(arrOff+3)*nPoints+node];

        real vx_m_new  = ( (n-1)*vx_m_old + vx  )*inv_n;
        real vy_m_new  = ( (n-1)*vy_m_old + vy  )*inv_n;
        real vz_m_new  = ( (n-1)*vz_m_old + vz  )*inv_n;
        real rho_m_new = ( (n-1)*rho_m_old+ rho )*inv_n;

        quantityArray[(arrOff+0)*nPoints+node] = vx_m_new;
        quantityArray[(arrOff+1)*nPoints+node] = vy_m_new;
        quantityArray[(arrOff+2)*nPoints+node] = vz_m_new;
        quantityArray[(arrOff+3)*nPoints+node] = rho_m_new;
    
        if(quantities[int(Statistic::Variances)])
        {
            arrOff = quantityArrayOffsets[int(Statistic::Variances)];

            real vx_var_old  = quantityArray[(arrOff+0)*nPoints+node];
            real vy_var_old  = quantityArray[(arrOff+1)*nPoints+node];
            real vz_var_old  = quantityArray[(arrOff+2)*nPoints+node];
            real rho_var_old = quantityArray[(arrOff+3)*nPoints+node];

            real vx_var_new  = ( (n-1)*(vx_var_old )+(vx  - vx_m_old )*(vx  - vx_m_new ) )*inv_n;
            real vy_var_new  = ( (n-1)*(vy_var_old )+(vy  - vy_m_old )*(vy  - vy_m_new ) )*inv_n;
            real vz_var_new  = ( (n-1)*(vz_var_old )+(vz  - vz_m_old )*(vz  - vz_m_new ) )*inv_n;
            real rho_var_new = ( (n-1)*(rho_var_old)+(rho - rho_m_old)*(rho - rho_m_new) )*inv_n;

            quantityArray[(arrOff+0)*nPoints+node] = vx_var_new;
            quantityArray[(arrOff+1)*nPoints+node] = vy_var_new;
            quantityArray[(arrOff+2)*nPoints+node] = vz_var_new;
            quantityArray[(arrOff+3)*nPoints+node] = rho_var_new; 
        }
    }
}

__global__ void calcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                )
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

    if(node>=nPoints) return;

    // Get indices of neighbor nodes. 
    // node referring to BSW cell as seen from probe point
    uint k = pointIndices[node];
    real u_interpX, u_interpY, u_interpZ, rho_interp;

    u_interpX = vx[k];
    u_interpY = vy[k];
    u_interpZ = vz[k];
    rho_interp = rho[k];

    calculatePointwiseQuantities(n, quantityArray, quantities, quantityArrayOffsets, nPoints, node, u_interpX, u_interpY, u_interpZ, rho_interp);

}

__global__ void interpAndCalcQuantitiesKernel(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* distX, real* distY, real* distZ,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray
                                )
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

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

    calculatePointwiseQuantities(n, quantityArray, quantities, quantityArrayOffsets, nPoints, node, u_interpX, u_interpY, u_interpZ, rho_interp);

}

bool Probe::getHasDeviceQuantityArray(){ return this->hasDeviceQuantityArray; }

void Probe::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    this->velocityRatio      = para->getVelocityRatio();
    this->densityRatio       = para->getDensityRatio();
    this->forceRatio         = para->getForceRatio();
    this->stressRatio        = para->getDensityRatio()*pow(para->getVelocityRatio(), 2.0);
    this->accelerationRatio = para->getVelocityRatio()/para->getTimeRatio();

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
        
        this->addProbeStruct(cudaManager, probeIndices_level, 
                            distX_level, distY_level, distZ_level, 
                            pointCoordsX_level, pointCoordsY_level, pointCoordsZ_level, 
                            level);
    }
}

void Probe::addProbeStruct(CudaMemoryManager* cudaManager, std::vector<int>& probeIndices,
                                      std::vector<real>& distX, std::vector<real>& distY, std::vector<real>& distZ,   
                                      std::vector<real>& pointCoordsX, std::vector<real>& pointCoordsY, std::vector<real>& pointCoordsZ,
                                      int level)
{
    probeParams[level] = SPtr<ProbeStruct>(new ProbeStruct);
    probeParams[level]->vals = 1;
    probeParams[level]->nPoints  = uint(pointCoordsX.size()); // Note, need to have both nPoints and nIndices because they differ in PlanarAverage
    probeParams[level]->nIndices = uint(probeIndices.size());

    probeParams[level]->pointCoordsX = (real*)malloc(probeParams[level]->nPoints*sizeof(real));
    probeParams[level]->pointCoordsY = (real*)malloc(probeParams[level]->nPoints*sizeof(real));
    probeParams[level]->pointCoordsZ = (real*)malloc(probeParams[level]->nPoints*sizeof(real));

    std::copy(pointCoordsX.begin(), pointCoordsX.end(), probeParams[level]->pointCoordsX);
    std::copy(pointCoordsY.begin(), pointCoordsY.end(), probeParams[level]->pointCoordsY);
    std::copy(pointCoordsZ.begin(), pointCoordsZ.end(), probeParams[level]->pointCoordsZ);

    // Note, dist only needed for kernels that do interpolate
    if( distX.size()>0 && distY.size()>0 && distZ.size()>0 )
    {
        probeParams[level]->hasDistances=true;
        cudaManager->cudaAllocProbeDistances(this, level);
        std::copy(distX.begin(), distX.end(), probeParams[level]->distXH);
        std::copy(distY.begin(), distY.end(), probeParams[level]->distYH);
        std::copy(distZ.begin(), distZ.end(), probeParams[level]->distZH);
        cudaManager->cudaCopyProbeDistancesHtoD(this, level);
    }  
    
    cudaManager->cudaAllocProbeIndices(this, level);
    std::copy(probeIndices.begin(), probeIndices.end(), probeParams[level]->pointIndicesH);
    cudaManager->cudaCopyProbeIndicesHtoD(this, level);

    uint arrOffset = 0;

    cudaManager->cudaAllocProbeQuantitiesAndOffsets(this, level);

    for( int var=0; var<int(Statistic::LAST); var++)
    {
        if(this->quantities[var])
        {
            probeParams[level]->quantitiesH[var] = true;
            probeParams[level]->arrayOffsetsH[var] = arrOffset;
            arrOffset += uint( this->getPostProcessingVariables(static_cast<Statistic>(var)).size() ); 
        }
    }
    
    cudaManager->cudaCopyProbeQuantitiesAndOffsetsHtoD(this, level);

    probeParams[level]->nArrays = arrOffset;

    cudaManager->cudaAllocProbeQuantityArray(this, level);

    for(uint arr=0; arr<probeParams[level]->nArrays; arr++)
    {
        for( uint point=0; point<probeParams[level]->nPoints; point++)
        {
            probeParams[level]->quantitiesArrayH[arr*probeParams[level]->nPoints+point] = 0.0f;
        }
    }
    if(this->hasDeviceQuantityArray)
        cudaManager->cudaCopyProbeQuantityArrayHtoD(this, level);
}

void Probe::interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t)
{
    if(max(int(t) - int(this->tStartAvg), -1) % this->tAvg==0)
    {
        SPtr<ProbeStruct> probeStruct = this->getProbeStruct(level);

        this->calculateQuantities(probeStruct, para, t, level);
        if(t>=this->tStartTmpAveraging) probeStruct->vals++;
    }

    if(max(int(t) - int(this->tStartOut), -1) % this->tOut == 0)
    {
        if(this->hasDeviceQuantityArray)
            cudaManager->cudaCopyProbeQuantityArrayDtoH(this, level);
        this->write(para, level, t);
    }
}

void Probe::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    for(int level=0; level<=para->getMaxLevel(); level++)
    {   
        if(this->probeParams[level]->hasDistances)
            cudaManager->cudaFreeProbeDistances(this, level);
        cudaManager->cudaFreeProbeIndices(this, level);
        cudaManager->cudaFreeProbeQuantityArray(this, level);
        cudaManager->cudaFreeProbeQuantitiesAndOffsets(this, level);
    }
}

void Probe::addStatistic(Statistic variable)
{
    assert(this->isAvailableStatistic(variable));

    this->quantities[int(variable)] = true;
    switch(variable)
    {
        case Statistic::Variances: 
            this->addStatistic(Statistic::Means); break;

        default: break;
    }
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
    int t_write = this->fileNameLU ? t: t/this->tOut; 

    const uint numberOfParts = this->getProbeStruct(level)->nPoints / para->getlimitOfNodesForVTK() + 1;

    std::vector<std::string> fnames;
    for (uint i = 1; i <= numberOfParts; i++)
	{
        std::string fname = this->probeName + "_bin_lev_" + StringUtil::toString<int>(level)
                                         + "_ID_" + StringUtil::toString<int>(para->getMyID())
                                         + "_Part_" + StringUtil::toString<int>(i);
        if(!this->outputTimeSeries) fname += "_t_" + StringUtil::toString<int>(t_write);
        fname += ".vtk";
		fnames.push_back(fname);
        this->fileNamesForCollectionFile.push_back(fname);
    }
    this->writeGridFiles(para, level, fnames, t);

    if(level == 0 && !this->outputTimeSeries) this->writeCollectionFile(para, t);
}

void Probe::writeCollectionFile(Parameter* para, int t)
{
    int t_write = this->fileNameLU ? t: t/this->tOut; 
    std::string filename = this->probeName + "_bin_ID_" + StringUtil::toString<int>(para->getMyID()) 
                                           + "_t_" + StringUtil::toString<int>(t_write) 
                                           + ".vtk";

    std::ofstream file;

    file.open(this->outputPath + "/" + filename + ".pvtu" );

    //////////////////////////////////////////////////////////////////////////
    
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
    file << "  <PUnstructuredGrid GhostLevel=\"1\">" << std::endl;

    file << "    <PPointData>" << std::endl;

    for(std::string varName: this->getVarNames()) //TODO
    {
        file << "       <DataArray type=\"Float64\" Name=\""<< varName << "\" /> " << std::endl;
    }
    file << "    </PPointData>" << std::endl;

    file << "    <PPoints>" << std::endl;
    file << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << std::endl;
    file << "    </PPoints>" << std::endl;

    for( auto& fname : this->fileNamesForCollectionFile )
    {
        const auto filenameWithoutPath=fname.substr( fname.find_last_of('/') + 1 );
        file << "    <Piece Source=\"" << filenameWithoutPath << ".bin.vtu\"/>" << std::endl;
    }

    file << "  </PUnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;

    //////////////////////////////////////////////////////////////////////////

    file.close();

    this->fileNamesForCollectionFile.clear();
}

void Probe::writeGridFiles(Parameter* para, int level, std::vector<std::string>& fnames, int t)
{
    std::vector< UbTupleFloat3 > nodes;
    std::vector< std::string > nodedatanames = this->getVarNames();

    uint startpos = 0;
    uint endpos = 0;
    uint sizeOfNodes = 0;
    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    SPtr<ProbeStruct> probeStruct = this->getProbeStruct(level);

    for (uint part = 0; part < fnames.size(); part++)
    {        
        startpos = part * para->getlimitOfNodesForVTK();
        uint nDataPoints = this->outputTimeSeries? probeStruct->vals: probeStruct->nPoints;
        sizeOfNodes = min(para->getlimitOfNodesForVTK(), nDataPoints - startpos);
        endpos = startpos + sizeOfNodes;

        //////////////////////////////////////////////////////////////////////////
        nodes.resize(sizeOfNodes);

        for (uint pos = startpos; pos < endpos; pos++)
        {
            nodes[pos-startpos] = makeUbTuple(  float(probeStruct->pointCoordsX[pos]),
                                                float(probeStruct->pointCoordsY[pos]),
                                                float(probeStruct->pointCoordsZ[pos]));
        }

        for( auto it=nodedata.begin(); it!=nodedata.end(); it++) it->resize(sizeOfNodes);

        for( int var=0; var < int(Statistic::LAST); var++){           
            if(this->quantities[var])
            {
                Statistic statistic = static_cast<Statistic>(var);
                real coeff;

                std::vector<PostProcessingVariable> postProcessingVariables = this->getPostProcessingVariables(statistic);
                uint n_arrs = uint(postProcessingVariables.size());

                uint arrOff = probeStruct->arrayOffsetsH[var];
                uint arrLen = probeStruct->nPoints;

                for(uint arr=0; arr<n_arrs; arr++)
                {
                    coeff = postProcessingVariables[arr].conversionFactor;
                    
                    for (uint pos = startpos; pos < endpos; pos++)
                    {
                        nodedata[arrOff+arr][pos-startpos] = double(probeStruct->quantitiesArrayH[(arrOff+arr)*arrLen+pos]*coeff);
                    }
                }
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(this->outputPath + "/" + fnames[part], nodes, nodedatanames, nodedata);
    }
}

std::vector<std::string> Probe::getVarNames()
{
    std::vector<std::string> varNames;
    for( int statistic=0; statistic < int(Statistic::LAST); statistic++)
    {
        if(this->quantities[statistic])
        {
            std::vector<PostProcessingVariable> postProcessingVariables = this->getPostProcessingVariables(static_cast<Statistic>(statistic));            
            for(int i = 0; i<postProcessingVariables.size(); i++) 
                varNames.push_back(postProcessingVariables[i].name);
        }
    }
    return varNames;
}