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


std::vector<std::string> getPostProcessingVariableNames(PostProcessingVariable variable)
{
    std::vector<std::string> varNames;
    switch (variable)
    {
    case PostProcessingVariable::Instantaneous:
        varNames.push_back("vx");
        varNames.push_back("vy");
        varNames.push_back("vz");
        varNames.push_back("rho");
        break;
    case PostProcessingVariable::Means:
        varNames.push_back("vx_mean");
        varNames.push_back("vy_mean");
        varNames.push_back("vz_mean");
        varNames.push_back("rho_mean");
        break;
    case PostProcessingVariable::Variances:
        varNames.push_back("vx_var");
        varNames.push_back("vy_var");
        varNames.push_back("vz_var");
        varNames.push_back("rho_var");
        break;
    default:
        break;
    }
    return varNames;
}

__device__ void calculateQuantities(uint n, real* quantityArray, bool* quantities, uint* quantityArrayOffsets, uint nPoints, uint node, real vx, real vy, real vz, real rho)
{
    //"https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm"
    // also has extensions for higher order and covariances
    real inv_n = 1/real(n);

    if(quantities[int(PostProcessingVariable::Instantaneous)])
    {
        uint arrOff = quantityArrayOffsets[int(PostProcessingVariable::Instantaneous)];
        quantityArray[(arrOff+0)*nPoints+node] = vx;
        quantityArray[(arrOff+1)*nPoints+node] = vy;
        quantityArray[(arrOff+2)*nPoints+node] = vz;
        quantityArray[(arrOff+3)*nPoints+node] = rho;
    }

    if(quantities[int(PostProcessingVariable::Means)])
    {
        
        uint arrOff = quantityArrayOffsets[int(PostProcessingVariable::Means)];
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
    
        if(quantities[int(PostProcessingVariable::Variances)])
        {
            arrOff = quantityArrayOffsets[int(PostProcessingVariable::Variances)];

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

__global__ void interpQuantities(   uint* pointIndices,
                                    uint nPoints, uint n,
                                    real* distX, real* distY, real* distZ,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    bool* quantities,
                                    uint* quantityArrayOffsets, real* quantityArray,
                                    bool interpolate
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

    if(interpolate)
    {
        uint ke, kn, kt, kne, kte, ktn, ktne;
        getNeighborIndicesOfBSW(  k, ke, kn, kt, kne, kte, ktn, ktne, neighborX, neighborY, neighborZ);

        // Trilinear interpolation of macroscopic quantities to probe point
        real dW, dE, dN, dS, dT, dB;
        getInterpolationWeights(dW, dE, dN, dS, dT, dB, distX[node], distY[node], distZ[node]);


        u_interpX  = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vx );
        u_interpY  = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vy );
        u_interpZ  = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vz );
        rho_interp = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, rho );
    }
    else
    {
        u_interpX = vx[k];
        u_interpY = vy[k];
        u_interpZ = vz[k];
        rho_interp = rho[k];
    }

    calculateQuantities(n, quantityArray, quantities, quantityArrayOffsets, nPoints, node, u_interpX, u_interpY, u_interpZ, rho_interp);

}

void Probe::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{

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
    probeParams[level]->nPoints = uint(probeIndices.size());

    probeParams[level]->pointCoordsX = (real*)malloc(probeParams[level]->nPoints*sizeof(real));
    probeParams[level]->pointCoordsY = (real*)malloc(probeParams[level]->nPoints*sizeof(real));
    probeParams[level]->pointCoordsZ = (real*)malloc(probeParams[level]->nPoints*sizeof(real));

    std::copy(pointCoordsX.begin(), pointCoordsX.end(), probeParams[level]->pointCoordsX);
    std::copy(pointCoordsY.begin(), pointCoordsY.end(), probeParams[level]->pointCoordsY);
    std::copy(pointCoordsZ.begin(), pointCoordsZ.end(), probeParams[level]->pointCoordsZ);

    // Might have to catch nPoints=0 ?!?!
    cudaManager->cudaAllocProbeDistances(this, level);
    cudaManager->cudaAllocProbeIndices(this, level);

    std::copy(distX.begin(), distX.end(), probeParams[level]->distXH);
    std::copy(distY.begin(), distY.end(), probeParams[level]->distYH);
    std::copy(distZ.begin(), distZ.end(), probeParams[level]->distZH);
    std::copy(probeIndices.begin(), probeIndices.end(), probeParams[level]->pointIndicesH);

    cudaManager->cudaCopyProbeDistancesHtoD(this, level);
    cudaManager->cudaCopyProbeIndicesHtoD(this, level);

    uint arrOffset = 0;

    cudaManager->cudaAllocProbeQuantitiesAndOffsets(this, level);

    for( int var=0; var<int(PostProcessingVariable::LAST); var++){
    if(this->quantities[var])
    {

        probeParams[level]->quantitiesH[var] = true;
        probeParams[level]->arrayOffsetsH[var] = arrOffset;
        arrOffset += uint(getPostProcessingVariableNames(static_cast<PostProcessingVariable>(var)).size());
    }}
    
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
    cudaManager->cudaCopyProbeQuantityArrayHtoD(this, level);
}

void Probe::interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t)
{

    if(t>this->tStartAvg)
    {
        SPtr<ProbeStruct> probeStruct = this->getProbeStruct(level);

        this->calculateQuantities(probeStruct, para, level);
        probeStruct->vals++;

        if(max(int(t) - int(this->tStartOut), -1) % this->tOut == 0)
        {
            cudaManager->cudaCopyProbeQuantityArrayDtoH(this, level);

            this->write(para, level, t);
        }

    }
}

void Probe::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        cudaManager->cudaFreeProbeDistances(this, level);
        cudaManager->cudaFreeProbeIndices(this, level);
        cudaManager->cudaFreeProbeQuantityArray(this, level);
        cudaManager->cudaFreeProbeQuantitiesAndOffsets(this, level);
    }
}

void Probe::addPostProcessingVariable(PostProcessingVariable variable)
{
    this->quantities[int(variable)] = true;
    switch(variable)
    {
        case PostProcessingVariable::Variances: 
            this->addPostProcessingVariable(PostProcessingVariable::Means); break;
        default: break;
    }
}

std::string Probe::makeParallelFileName(int id, int t)
{
    return this->probeName + "_bin_ID_" + StringUtil::toString<int>(id) 
                                           + "_t_" + StringUtil::toString<int>(t) 
                                           + ".vtk";
}

std::string Probe::makeGridFileName(int level, int id, int t, uint part)
{
    return this->probeName + "_bin_lev_" + StringUtil::toString<int>(level)
                                         + "_ID_" + StringUtil::toString<int>(id)
                                         + "_Part_" + StringUtil::toString<int>(part) 
                                         + "_t_" + StringUtil::toString<int>(t) 
                                         + ".vtk";
}

void Probe::write(Parameter* para, int level, int t)
{
    const uint numberOfParts = this->getProbeStruct(level)->nPoints / para->getlimitOfNodesForVTK() + 1;

    std::vector<std::string> fnames;
    for (uint i = 1; i <= numberOfParts; i++)
	{
        this->writeGridFile(para, level, t, i);
    }
    if(level == 0) this->writeParallelFile(para, t);
}

void Probe::writeParallelFile(Parameter* para, int t)
{
    std::string filename = this->outputPath + "/" + this->makeParallelFileName(para->getMyID(), t);

    std::vector<std::string> cellNames;

    getWriter()->writeParallelFile(filename, fileNamesForCollectionFile, varNames, cellNames);

    this->fileNamesForCollectionFile.clear();
}

void Probe::writeGridFile(Parameter* para, int level, int t, uint part)
{
    std::string fname = this->outputPath + "/" + this->makeGridFileName(level, para->getMyID(), t, part);

    std::vector< UbTupleFloat3 > nodes;
    std::vector< std::string > nodedatanames = this->getVarNames();

    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    SPtr<ProbeStruct> probeStruct = this->getProbeStruct(level);

    uint startpos = part * para->getlimitOfNodesForVTK();
    uint sizeOfNodes = min(para->getlimitOfNodesForVTK(), probeStruct->nPoints - startpos);
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

    for( int var=0; var < int(PostProcessingVariable::LAST); var++){
    if(this->quantities[var])
    {
        PostProcessingVariable quantity = static_cast<PostProcessingVariable>(var);
        real coeff;
        uint n_arrs = uint(getPostProcessingVariableNames(quantity).size());

        switch(quantity)
        {
        case PostProcessingVariable::Instantaneous:
            coeff = para->getVelocityRatio();
        break;
        case PostProcessingVariable::Means:
            coeff = para->getVelocityRatio();
        break;
        case PostProcessingVariable::Variances:
            coeff = pow(para->getVelocityRatio(),2);
        break;
        default: break;
        }

        uint arrOff = probeStruct->arrayOffsetsH[var];
        uint arrLen = probeStruct->nPoints;

        for(uint arr=0; arr<n_arrs; arr++)
        {
            for (uint pos = startpos; pos < endpos; pos++)
            {
                nodedata[arrOff+arr][pos-startpos] = double(probeStruct->quantitiesArrayH[(arrOff+arr)*arrLen+pos]*coeff);
            }
        }
    }}
    this->fileNamesForCollectionFile.push_back(getWriter()->writeNodesWithNodeData(fname, nodes, nodedatanames, nodedata));
}

std::vector<std::string> Probe::getVarNames()
{
    std::vector<std::string> varNames;
    for( int var=0; var < int(PostProcessingVariable::LAST); var++){
    if(this->quantities[var])
    {
        std::vector<std::string> names = getPostProcessingVariableNames(static_cast<PostProcessingVariable>(var));
        varNames.insert(varNames.end(), names.begin(), names.end());
    }}
    return varNames;
}
