#include "PrecursorWriter.h"
#include "basics/writer/WbWriterVtkXmlImageBinary.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuda/CudaGrid.h>
#include "Kernel/Utilities/DistributionHelper.cuh"

#include <Core/StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"

using namespace vf::lbm::dir;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO check everything for multiple level
void index1d(int& idx, int y, int z, int ny, int nz)
{
    idx = y+ny*z;
}

void index2d(int idx, int& y, int& z, int ny, int nz)
{
    z = idx/ny;
    y = idx-ny*z;
}

__inline__ __host__ __device__ uint lIndex(const uint component, const uint node, const uint timestep, const uint nComponents, const uint nNodes)
{
    return node+nNodes*(component+timestep*nComponents);
}

__inline__ __host__ __device__ uint lIndex(const uint component, const uint node, const uint nNodes)
{
    return node+component*nNodes;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void fillArrayVelocities(const uint nNodes, 
                                    uint* indices, 
                                    real *precursorData,
                                    real *vx,
                                    real *vy,
                                    real *vz,
                                    real velocityRatio)


{
    const uint node = vf::gpu::getNodeIndex();

    if(node>=nNodes) return;

    precursorData[lIndex(0u, node, nNodes)] = vx[indices[node]]*velocityRatio;
    precursorData[lIndex(1u, node, nNodes)] = vy[indices[node]]*velocityRatio;
    precursorData[lIndex(2u, node, nNodes)] = vz[indices[node]]*velocityRatio;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void fillArrayDistributions( uint nNodes, uint* indices, 
                                        real* precursorData,
                                        real* distributions,
                                        uint* neighborX, uint* neighborY, uint* neighborZ,
                                        bool isEvenTimestep,
                                        unsigned long numberOfLBnodes)
{
    const uint node = vf::gpu::getNodeIndex();

    if(node>=nNodes) return;

    Distributions27 dist = vf::gpu::getDistributionReferences27(distributions, numberOfLBnodes, isEvenTimestep);
    
    ////////////////////////////////////////////////////////////////////////////////
    // ! - Set neighbor indices (necessary for indirect addressing)
    uint k_000 = indices[node];
    // uint k_M00 = neighborX[k_000];
    uint k_0M0 = neighborY[k_000];
    uint k_00M = neighborZ[k_000];
    // uint k_MM0 = neighborY[k_M00];
    // uint k_M0M = neighborZ[k_M00];
    uint k_0MM = neighborZ[k_0M0];
    // uint k_MMM = neighborZ[k_MM0];

    ////////////////////////////////////////////////////////////////////////////////////
    //! - Get local distributions in PX directions
    //!
    precursorData[lIndex(PrecP00, node, nNodes)] = (dist.f[DIR_P00])[k_000];
    precursorData[lIndex(PrecPP0, node, nNodes)] = (dist.f[DIR_PP0])[k_000];
    precursorData[lIndex(PrecPM0, node, nNodes)] = (dist.f[DIR_PM0])[k_0M0];
    precursorData[lIndex(PrecP0P, node, nNodes)] = (dist.f[DIR_P0P])[k_000];
    precursorData[lIndex(PrecP0M, node, nNodes)] = (dist.f[DIR_P0M])[k_00M];
    precursorData[lIndex(PrecPPP, node, nNodes)] = (dist.f[DIR_PPP])[k_000];
    precursorData[lIndex(PrecPMP, node, nNodes)] = (dist.f[DIR_PMP])[k_0M0];
    precursorData[lIndex(PrecPPM, node, nNodes)] = (dist.f[DIR_PPM])[k_00M];
    precursorData[lIndex(PrecPMM, node, nNodes)] = (dist.f[DIR_PMM])[k_0MM];
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrecursorWriter::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    VF_LOG_INFO("PrecursorWriter: Start initializing...");
    VF_LOG_INFO("Writing yz-planes at x={}m every {}. timestep, starting at t={}", this->xPos, this->tSave, this->tStartOut);

    precursorStructs.resize(para->getMaxLevel()+1);
    for(int level=0; level<=para->getMaxLevel(); level++)
    {

        real dx = abs(para->getParH(level)->coordinateX[1]-para->getParH(level)->coordinateX[para->getParH(level)->neighborX[1]]);
        int maxPoints = (int((yMax-yMin)/dx)+1)* (int((zMax-zMin)/dx)+1);

        real lowestY, lowestZ, highestY, highestZ;

        lowestY = para->getParH(level)->coordinateY[para->getParH(level)->numberOfNodes-1];
        highestY = para->getParH(level)->coordinateY[1];        
        
        lowestZ = para->getParH(level)->coordinateZ[para->getParH(level)->numberOfNodes-1];
        highestZ = para->getParH(level)->coordinateZ[1];

        std::vector<uint> indicesOnGrid;
        std::vector<int> indicesOnPlane;
        std::vector<real> coordY, coordZ;

        for(uint j=1; j<para->getParH(level)->numberOfNodes; j++ )
        {
            real pointCoordX = para->getParH(level)->coordinateX[j];
            real pointCoordY = para->getParH(level)->coordinateY[j];
            real pointCoordZ = para->getParH(level)->coordinateZ[j];
            if( para->getParH(level)->typeOfGridNode[j] == GEO_FLUID &&
                pointCoordX < (dx+xPos) && pointCoordX >= xPos       &&
                pointCoordY<=yMax && pointCoordY>=yMin               && 
                pointCoordZ<=zMax && pointCoordZ>=zMin)
            {
                highestY = max(highestY, pointCoordY);
                highestZ = max(highestZ, pointCoordZ);

                lowestY = min(lowestY, pointCoordY);
                lowestZ = min(lowestZ, pointCoordZ);
                indicesOnGrid.push_back(j);    
                coordY.push_back(pointCoordY);            
                coordZ.push_back(pointCoordZ);    
            }
        }
        if(indicesOnGrid.size()==0)
            throw std::runtime_error("PrecursorWriter did not find any points on the grid");

        int ny = int((highestY-lowestY)/dx)+1;
        int nz = int((highestZ-lowestZ)/dx)+1;

        for(uint i=0;i<indicesOnGrid.size(); i++)
        {
                int idxY = int((coordY[i]-lowestY)/dx);
                int idxZ = int((coordZ[i]-lowestZ)/dx);
                int idx;
                index1d(idx, idxY, idxZ, ny, nz);
                indicesOnPlane.push_back(idx);
        }

        precursorStructs[level] = SPtr<PrecursorStruct>(new PrecursorStruct);
        precursorStructs[level]->nPoints = (uint)indicesOnGrid.size();
        precursorStructs[level]->indicesOnPlane = (int*) malloc(precursorStructs[level]->nPoints*sizeof(int));
        precursorStructs[level]->spacing = makeUbTuple(dx, dx, tSave*para->getTimeRatio()*pow(2,-level));
        precursorStructs[level]->origin = makeUbTuple(lowestY, lowestZ);
        precursorStructs[level]->extent = makeUbTuple(0, ny-1, 0, nz-1);
        precursorStructs[level]->nPointsInPlane = ny*nz;
        precursorStructs[level]->timestepsPerFile = min(para->getlimitOfNodesForVTK()/(ny*nz), maxtimestepsPerFile);
        precursorStructs[level]->filesWritten = 0;
        precursorStructs[level]->timestepsBuffered = 0;
        
        switch (outputVariable)
        {
        case OutputVariable::Velocities:
            precursorStructs[level]->nQuantities = 3;
            break;
        case OutputVariable::Distributions:
            precursorStructs[level]->nQuantities = 9;
            break;
        
        default:
            break;
        }

        cudaManager->cudaAllocPrecursorWriter(this, level);
    
        std::copy(indicesOnGrid.begin(), indicesOnGrid.end(), precursorStructs[level]->indicesH);
        std::copy(indicesOnPlane.begin(), indicesOnPlane.end(), precursorStructs[level]->indicesOnPlane);

        cudaManager->cudaCopyPrecursorWriterIndicesHtoD(this, level);

        VF_LOG_INFO("Found {} points in precursor plane on level {}", precursorStructs[level]->nPoints, level);
    }
    VF_LOG_INFO("PrecursorWriter: Done initializing.");
}


void PrecursorWriter::interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t)
{
    uint t_level         = para->getTimeStep(level, t, true);
    uint tStartOut_level = tStartOut*pow(2, level);
    uint tEnd_level      = para->getTimestepEnd()*pow(2, level);

    if(t_level>tStartOut_level && ((t_level-tStartOut_level) % tSave)==0)
    {
        vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, precursorStructs[level]->nPoints);

        if(this->outputVariable==OutputVariable::Velocities)
        {
            fillArrayVelocities<<<grid.grid, grid.threads>>>(   precursorStructs[level]->nPoints, precursorStructs[level]->indicesD, 
                                                                precursorStructs[level]->bufferD, 
                                                                para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                para->getVelocityRatio());
            getLastCudaError("In PrecursorWriter::interact fillArrayVelocities execution failed");
        }
        else if(this->outputVariable==OutputVariable::Distributions)
        {
            fillArrayDistributions<<<grid.grid, grid.threads>>>(precursorStructs[level]->nPoints, precursorStructs[level]->indicesD, 
                                                                precursorStructs[level]->bufferD,
                                                                para->getParD(level)->distributions.f[0],
                                                                para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                para->getEvenOrOdd(level), para->getParD(level)->numberOfNodes);
            getLastCudaError("In PrecursorWriter::interact fillArrayDistributions execution failed");
        }
        cudaManager->cudaCopyPrecursorWriterOutputVariablesDtoH(this, level);

        // switch device buffer and data pointer so precursor data is gathered in buffer and copied from bufferD to bufferH
        real *tmp = precursorStructs[level]->bufferD;
        precursorStructs[level]->bufferD = precursorStructs[level]->dataD;
        precursorStructs[level]->dataD = tmp;

        precursorStructs[level]->timestepsBuffered++;

        if(precursorStructs[level]->timestepsBuffered >= precursorStructs[level]->timestepsPerFile || t == para->getTimestepEnd())
        {
        // switch host buffer and data pointer so precursor data is copied in buffer and written from data

            tmp = precursorStructs[level]->bufferH;
            precursorStructs[level]->bufferH = precursorStructs[level]->dataH;
            precursorStructs[level]->dataH = tmp;

            writeFuture.wait();
            writeFuture = std::async(std::launch::async, [this](Parameter* para, uint level, uint timesteps){ this->write(para, level, timesteps); }, para, level, precursorStructs[level]->timestepsBuffered);
            precursorStructs[level]->timestepsBuffered = 0;
        }
    }
}


void PrecursorWriter::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    writeFuture.wait();
    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        if(getPrecursorStruct(level)->timestepsBuffered>0)
            write(para, level, getPrecursorStruct(level)->timestepsBuffered);

        cudaManager->cudaFreePrecursorWriter(this, level);
    }
}


void PrecursorWriter::write(Parameter* para, int level, uint timestepsBuffered)
{
    std::string fname = this->makeFileName(fileName, level, para->getMyProcessID(), precursorStructs[level]->filesWritten) + getWriter()->getFileExtension();
    std::string wholeName = outputPath + "/" + fname;

    uint nPointsInPlane = precursorStructs[level]->nPointsInPlane;

    int startTime = precursorStructs[level]->filesWritten*precursorStructs[level]->timestepsPerFile;

    UbTupleInt6 extent = makeUbTuple(   val<1>(precursorStructs[level]->extent),    val<2>(precursorStructs[level]->extent), 
                                        val<3>(precursorStructs[level]->extent),    val<4>(precursorStructs[level]->extent), 
                                        startTime,                          startTime+(int)timestepsBuffered-1);

    UbTupleFloat3 origin = makeUbTuple( val<1>(precursorStructs[level]->origin), val<2>(precursorStructs[level]->origin), 0.f);

    std::vector<std::vector<double>> nodedata;
    
    for(uint quant=0; quant<precursorStructs[level]->nQuantities; quant++)
    {
        std::vector<double> doubleArr(nPointsInPlane*timestepsBuffered, NAN);
        for( uint timestep=0; timestep<timestepsBuffered; timestep++)
        {
            for (uint pos=0; pos < precursorStructs[level]->nPoints; pos++)
            {
                int indexOnPlane = precursorStructs[level]->indicesOnPlane[pos]+timestep*nPointsInPlane;
                doubleArr[indexOnPlane] = double(precursorStructs[level]->dataH[lIndex(quant, pos, timestep, precursorStructs[level]->nQuantities, precursorStructs[level]->nPoints)]);
            }
        }
        nodedata.push_back(doubleArr);
    }

    std::vector<std::vector<double>> celldata;
    getWriter()->writeData(wholeName, nodedatanames, celldatanames, nodedata, celldata, extent, origin, precursorStructs[level]->spacing, extent);
    precursorStructs[level]->filesWritten++;
}

std::string PrecursorWriter::makeFileName(std::string fileName, int level, int id, uint filesWritten)
{
    return fileName + "_lev_" + StringUtil::toString<int>(level)
                    + "_ID_" + StringUtil::toString<int>(id)
                    + "_File_" + StringUtil::toString<int>(filesWritten);
}

void PrecursorWriter::getTaggedFluidNodes(Parameter *para, GridProvider* gridProvider)
{
    for(uint level=0; level<(uint)para->getMaxLevel(); level++)
    {
        if(outputVariable==OutputVariable::Velocities)
        {
            std::vector<uint> indices(precursorStructs[level]->indicesH, precursorStructs[level]->indicesH+precursorStructs[level]->nPoints);
            gridProvider->tagFluidNodeIndices(indices, CollisionTemplate::WriteMacroVars, level);
        }
    }
}