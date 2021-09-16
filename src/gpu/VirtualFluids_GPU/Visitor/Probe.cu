#include "Probe.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <unordered_set>

#include "VirtualFluids_GPU/GPU/GeometryUtils.h"
#include "Kernel/Utilities/CudaGrid.h"
#include "basics/writer/WbWriterVtkXmlBinary.h"
#include <Core/StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"


__global__ void interpQuantities(   uint* pointIndices,
                                    uint nPoints,
                                    real* distX, real* distY, real* distZ,
                                    real* vx, real* vy, real* vz, real* rho,            
                                    uint* neighborX, uint* neighborY, uint* neighborZ,
                                    // real* vx_point, real* vy_point, real* vz_point, real* rho_point,
                                    PostProcessingVariable* PostProcessingVariables,
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
    uint ke, kn, kt, kne, kte, ktn, ktne;
    getNeighborIndicesBSW(  k, ke, kn, kt, kne, kte, ktn, ktne, neighborX, neighborY, neighborZ);

    // Trilinear interpolation of macroscopic quantities to probe point
    real dW, dE, dN, dS, dT, dB;
    getInterpolationWeights(dW, dE, dN, dS, dT, dB, distX[node], distY[node], distZ[node]);

    // vx_point [node] = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vx );
    // vy_point [node] = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vy );
    // vz_point [node] = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vz );
    // rho_point[node] = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, rho );

    real u_interpX, u_interpY, u_interpZ, rho_interp;

    // printf("k %i, u %f \n",k, vx[k]);
    u_interpX = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vx );
    u_interpY = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vy );
    u_interpZ = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, vz );
    rho_interp = trilinearInterpolation( dW, dE, dN, dS, dT, dB, k, ke, kn, kt, kne, kte, ktn, ktne, rho );

    //TODO change computation of  means and variances to something more stable, see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm

    for( int variableIndex = 0; PostProcessingVariables[variableIndex] != PostProcessingVariable::LAST; variableIndex++)
    {
        PostProcessingVariable variable = PostProcessingVariables[variableIndex];
        uint arrayOffset = quantityArrayOffsets[int(variable)]*nPoints;

        switch(variable)
        {
            case PostProcessingVariable::Means:
            {
                // printf("u_interp: %f \n", u_interpX);
                quantityArray[arrayOffset+node] += u_interpX; arrayOffset += nPoints;
                quantityArray[arrayOffset+node] += u_interpY; arrayOffset += nPoints;
                quantityArray[arrayOffset+node] += u_interpZ; arrayOffset += nPoints;
                quantityArray[arrayOffset+node] += rho_interp;
            } break;
            case PostProcessingVariable::Variances:
            { 
                quantityArray[arrayOffset+node] += pow(u_interpX, 2.f); arrayOffset += nPoints;
                quantityArray[arrayOffset+node] += pow(u_interpY, 2.f); arrayOffset += nPoints;
                quantityArray[arrayOffset+node] += pow(u_interpZ, 2.f); arrayOffset += nPoints;
                quantityArray[arrayOffset+node] += pow(rho_interp, 2.f); 
            } break;
            default: break;
        }
    }
}


void Probe::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{

    probeParams.resize(para->getMaxLevel()+1);

    //Remove double entries
    std::unordered_set<PostProcessingVariable> s(this->postProcessingVariables.begin(), this->postProcessingVariables.end());
    this->postProcessingVariables.assign(s.begin(), s.end());
    this->addPostProcessingVariable(PostProcessingVariable::LAST);

    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        std::vector<int> probeIndices_level;
        std::vector<real> distX_level;
        std::vector<real> distY_level;
        std::vector<real> distZ_level;        
        std::vector<real> pointCoordsX_level;
        std::vector<real> pointCoordsY_level;
        std::vector<real> pointCoordsZ_level;
        real dx = abs(para->getParH(level)->coordX_SP[1]-para->getParH(level)->coordX_SP[para->getParH(level)->neighborX_SP[1]]);
        for(uint j=0; j<para->getParH(level)->size_Mat_SP; j++ )
        {    
            for(uint point=0; point<this->nProbePoints; point++)
            {
                real pointCoordX = this->pointCoordsX[point];
                real pointCoordY = this->pointCoordsY[point];
                real pointCoordZ = this->pointCoordsZ[point];
                real distX = pointCoordX-para->getParH(level)->coordX_SP[j];
                real distY = pointCoordY-para->getParH(level)->coordY_SP[j];
                real distZ = pointCoordZ-para->getParH(level)->coordZ_SP[j];
                if( distX <=dx && distY <=dx && distZ <=dx &&
                    distX >0.f && distY >0.f && distZ >0.f)
                {
                    probeIndices_level.push_back(j);
                    distX_level.push_back( distX/dx );
                    distY_level.push_back( distY/dx );
                    distZ_level.push_back( distZ/dx );
                    pointCoordsX_level.push_back( pointCoordX );
                    pointCoordsY_level.push_back( pointCoordY );
                    pointCoordsZ_level.push_back( pointCoordZ );
                    // printf("Found Point %i, x: %f, y: %f,z: %f, \n For %f %f %f, \n distx: %f, disty: %f, distz: %f \n", j, para->getParH(level)->coordX_SP[j],para->getParH(level)->coordY_SP[j],para->getParH(level)->coordZ_SP[j],
                    // this->pointCoordsX[point], this->pointCoordsY[point], this->pointCoordsZ[point], 
                    // distX, distY, distZ);
                }
            }
        }
        
        probeParams[level] = new ProbeStruct;
        probeParams[level]->nPoints = probeIndices_level.size();
        probeParams[level]->pointCoordsX = pointCoordsX_level.data();
        probeParams[level]->pointCoordsY = pointCoordsY_level.data();
        probeParams[level]->pointCoordsZ = pointCoordsZ_level.data();
        // Might have to catch nPoints=0 ?!?!
        cudaManager->cudaAllocProbeDistances(this, level);
        cudaManager->cudaAllocProbeIndices(this, level);

        std::copy(distX_level.begin(), distX_level.end(), probeParams[level]->distXH);
        std::copy(distY_level.begin(), distY_level.end(), probeParams[level]->distYH);
        std::copy(distZ_level.begin(), distZ_level.end(), probeParams[level]->distZH);
        std::copy(probeIndices_level.begin(), probeIndices_level.end(), probeParams[level]->pointIndicesH);

        cudaManager->cudaCopyProbeDistancesHtoD(this, level);
        cudaManager->cudaCopyProbeIndicesHtoD(this, level);

        uint arrOffset = 0;

        cudaManager->cudaAllocProbeQuantities(this, level);

        for(PostProcessingVariable variable: this->postProcessingVariables)
        {
            probeParams[level]->arrayOffsetsH[int(variable)] = arrOffset;
            switch(variable)
            {
                case PostProcessingVariable::Means: arrOffset += 4; break;                
                case PostProcessingVariable::Variances: arrOffset += 4; break;
                default: break;
            }
        }

        probeParams[level]->nArrays = arrOffset;
        cudaManager->cudaAllocProbeQuantityArray(this, level);
        std::copy(this->postProcessingVariables.begin(), this->postProcessingVariables.end(), probeParams[level]->quantitiesH);
        cudaManager->cudaCopyProbeQuantitiesHtoD(this, level);

        for(int arr=0; arr<probeParams[level]->nArrays; arr++)
        {
            for( int point=0; point<probeParams[level]->nPoints; point++)
            {
                probeParams[level]->quantitiesArrayH[arr*probeParams[level]->nPoints+point] = 0.0f;
            }
        }
        
        cudaManager->cudaCopyProbeQuantityArrayHtoD(this, level);
    }
}


void Probe::visit(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t)
{    
    ProbeStruct* probeStruct = this->getProbeStruct(level);

    vf::gpu::CudaGrid grid = vf::gpu::CudaGrid(probeStruct->nPoints, 128);

    interpQuantities<<<grid.grid, grid.threads>>>(  probeStruct->pointIndicesD, probeStruct->nPoints,
                                                    probeStruct->distXD, probeStruct->distYD, probeStruct->distZD,
                                                    para->getParD(level)->vx_SP, para->getParD(level)->vy_SP, para->getParD(level)->vz_SP, para->getParD(level)->rho_SP, 
                                                    para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP, 
                                                    probeStruct->quantitiesD, probeStruct->arrayOffsetsD, probeStruct->quantitiesArrayD  );
    if(max(int(t) - int(this->tStart), -1) % this->tOut == 0)
    {
        cudaManager->cudaCopyProbeQuantityArrayDtoH(this, level);
        this->write(para, level, t);
    }


}

void Probe::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        cudaManager->cudaFreeProbeDistances(this, level);
        cudaManager->cudaFreeProbeIndices(this, level);
        cudaManager->cudaFreeProbeQuantityArray(this, level);
        cudaManager->cudaFreeProbeQuantities(this, level);

    }
}

void Probe::setProbePointsFromList(std::vector<real> &_pointCoordsX, std::vector<real> &_pointCoordsY, std::vector<real> &_pointCoordsZ)
{
    bool isSameLength = ( (_pointCoordsX.size()==_pointCoordsY.size()) && (_pointCoordsY.size()==_pointCoordsZ.size()));
    assert("Probe: point lists have different lengths" && isSameLength);
    this->pointCoordsX = _pointCoordsX;
    this->pointCoordsY = _pointCoordsY;
    this->pointCoordsZ = _pointCoordsZ;
    this->nProbePoints = _pointCoordsX.size();
    printf("Added list of %u  points \n", this->nProbePoints );
}

void Probe::setProbePointsFromXNormalPlane(real pos_x, real pos0_y, real pos0_z, real pos1_y, real pos1_z, real delta_y, real delta_z)
{
    int n_points_y = int((pos1_y-pos0_y)/delta_y);
    int n_points_z = int((pos1_z-pos0_z)/delta_z);
    int n_points = n_points_y*n_points_z;

    std::vector<real> pointCoordsXtmp, pointCoordsYtmp, pointCoordsZtmp;
    pointCoordsXtmp.reserve(n_points);
    pointCoordsYtmp.reserve(n_points);
    pointCoordsZtmp.reserve(n_points);

    for(int n_y=0; n_y<n_points_y; n_y++)
    {
        for(int n_z=0; n_z<n_points_z; n_z++)
        {
            pointCoordsXtmp.push_back(pos_x);
            pointCoordsYtmp.push_back(pos0_y+delta_y*n_y);
            pointCoordsZtmp.push_back(pos0_z+delta_z*n_z);
        }
    }
    this->setProbePointsFromList(pointCoordsXtmp, pointCoordsYtmp, pointCoordsZtmp);
}

void Probe::addPostProcessingVariable(PostProcessingVariable _variable)
{
    this->postProcessingVariables.push_back(_variable);
    switch(_variable)
    {
        case PostProcessingVariable::Means: break;
        case PostProcessingVariable::Variances: this->postProcessingVariables.push_back(PostProcessingVariable::Means); break;
        default: break;
    }
}

void Probe::write(Parameter* para, int level, int t)
{
    const uint numberOfParts = this->getProbeStruct(level)->nPoints / para->getlimitOfNodesForVTK() + 1;


    std::vector<std::string> fnames;
    for (uint i = 1; i <= numberOfParts; i++)
	{
		fnames.push_back(this->probeName + "_bin_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(para->getMyID()) + "_Part_" + StringUtil::toString<int>(i) + "_t_" + StringUtil::toString<int>(t) + ".vtk");
        this->fileNamesForCollectionFile.push_back(fnames.back());
        this->writeGridFile(para, level, fnames, t);
    }

    this->writeCollectionFile(para, t);


}

void Probe::writeCollectionFile(Parameter* para, int t)
{
    std::string filename = this->probeName + "_bin_ID_" + StringUtil::toString<int>(para->getMyID()) + "_t_" + StringUtil::toString<int>(t) + ".vtk";

    std::ofstream file;

    file.open( filename + ".pvtu" );

    //////////////////////////////////////////////////////////////////////////
    
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
    file << "  <PUnstructuredGrid GhostLevel=\"1\">" << std::endl;

    file << "    <PPointData>" << std::endl;

    for(std::string varName: this->getVarNames())
    {
        file << "       <DataArray type=\"Float32\" Name=\""<< varName << "\" /> " << std::endl;
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

void Probe::writeGridFile(Parameter* para, int level, std::vector<std::string>& fnames, int t)
{
    std::vector< UbTupleFloat3 > nodes;
    std::vector< std::string > nodedatanames = this->getVarNames();

    uint startpos = 0;
    uint endpos = 0;
    uint sizeOfNodes = 0;
    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    real inv_t = 1/(real(max(t,1))*pow(2,level));
    for (uint part = 0; part < fnames.size(); part++)
    {        
        startpos = part * para->getlimitOfNodesForVTK();
        sizeOfNodes = min(para->getlimitOfNodesForVTK(), this->getProbeStruct(level)->nPoints - startpos);
        endpos = startpos + sizeOfNodes;

        //////////////////////////////////////////////////////////////////////////
        nodes.resize(sizeOfNodes);

        for (uint pos = startpos; pos < endpos; pos++)
        {
            nodes[pos-startpos] = makeUbTuple(  float(this->getProbeStruct(level)->pointCoordsX[pos]),
                                                float(this->getProbeStruct(level)->pointCoordsY[pos]),
                                                float(this->getProbeStruct(level)->pointCoordsZ[pos]));
        }

        for( auto it=nodedata.begin(); it!=nodedata.end(); it++) it->resize(sizeOfNodes);

        //TODO maybe change order of loops, could be faster, maybe not important, also still very ugly
        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            int dn = pos-startpos;
            //////////////////////////////////////////////////////////////////////////
            int offset = 0;
            for(PostProcessingVariable variable: this->postProcessingVariables)
            {
                int arrOffset = this->getProbeStruct(level)->arrayOffsetsH[int(variable)]*this->getProbeStruct(level)->nPoints;
                switch(variable)
                {
                    case PostProcessingVariable::Means:
                    {
                        nodedata[offset][dn] = (double)this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*para->getVelocityRatio()*inv_t; arrOffset+=probeParams[level]->nPoints; offset++;
                        nodedata[offset][dn] = (double)this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*para->getVelocityRatio()*inv_t; arrOffset+=probeParams[level]->nPoints; offset++;
                        nodedata[offset][dn] = (double)this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*para->getVelocityRatio()*inv_t; arrOffset+=probeParams[level]->nPoints; offset++;
                        nodedata[offset][dn] = (double)this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*para->getVelocityRatio()*inv_t; arrOffset+=probeParams[level]->nPoints; offset++;
                    } break;
                    case PostProcessingVariable::Variances:
                    {
                        int meansShift = this->getProbeStruct(level)->arrayOffsetsH[int(PostProcessingVariable::Means)]*this->getProbeStruct(level)->nPoints-arrOffset;
                        nodedata[offset][dn] = double(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*inv_t - pow(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+meansShift+pos]*inv_t,2))*pow(para->getVelocityRatio(),2); arrOffset+=probeParams[level]->nPoints; offset++;
                        nodedata[offset][dn] = double(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*inv_t - pow(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+meansShift+pos]*inv_t,2))*pow(para->getVelocityRatio(),2); arrOffset+=probeParams[level]->nPoints; offset++;
                        nodedata[offset][dn] = double(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*inv_t - pow(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+meansShift+pos]*inv_t,2))*pow(para->getVelocityRatio(),2); arrOffset+=probeParams[level]->nPoints; offset++;
                        nodedata[offset][dn] = double(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+pos]*inv_t - pow(this->getProbeStruct(level)->quantitiesArrayH[arrOffset+meansShift+pos]*inv_t,2))*pow(para->getVelocityRatio(),2); arrOffset+=probeParams[level]->nPoints; offset++;
                    } break;
                    default: break;
                }
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(fnames[part], nodes, nodedatanames, nodedata);
    }
}

std::vector<std::string> Probe::getVarNames()
{
    std::vector<std::string> varNames;
    for(PostProcessingVariable variable: this->postProcessingVariables)
    {
        switch(variable)
        {
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
            default: break;
        }
    }
    return varNames;
}
