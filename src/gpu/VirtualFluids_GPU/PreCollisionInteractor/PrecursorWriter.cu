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
            if( pointCoordX < (dx+xPos) && pointCoordX >= xPos &&
                pointCoordY<=yMax && pointCoordY>=yMin && 
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
        assert("PrecursorWriter did not find any points on the grid"&& indicesOnGrid.size()==0);
        int ny = int((highestY-lowestY)/dx)+1;
        int nz = int((highestZ-lowestZ)/dx)+1;

        for(uint i=0;i<indicesOnGrid.size(); i++)
        {
                int idxY = int((coordY[i]-lowestY)/dx);
                int idxZ = int((coordZ[i]-lowestZ)/dx);
                int idx;
                index1d(idx, idxY, idxZ, ny, nz);
                indicesOnPlane.push_back(idx);
                // printf("idx %d, idy %d, idz %d, ny %d, nz %d\n", idx, idxY, idxZ, ny, nz);
        }

        precursorStructs[level] = SPtr<PrecursorStruct>(new PrecursorStruct);
        precursorStructs[level]->nPoints = (uint)indicesOnGrid.size();
        precursorStructs[level]->indicesOnPlane = (int*) malloc(precursorStructs[level]->nPoints*sizeof(int));
        precursorStructs[level]->spacing = makeUbTuple(dx, dx, tSave*para->getTimeRatio());
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

        // printf("points %zu points on plane %zu \n",  indicesOnGrid.size(),  indicesOnPlane.size());

        cudaManager->cudaAllocPrecursorWriter(this, level);
    
        std::copy(indicesOnGrid.begin(), indicesOnGrid.end(), precursorStructs[level]->indicesH);
        std::copy(indicesOnPlane.begin(), indicesOnPlane.end(), precursorStructs[level]->indicesOnPlane);

        cudaManager->cudaCopyPrecursorWriterIndicesHtoD(this, level);
    }
}


void PrecursorWriter::interact(Parameter* para, CudaMemoryManager* cudaManager, int level, uint t)
{
    if(t>tStartOut ? ((t-tStartOut) % tSave)==0 : false)
    {
        SPtr<PrecursorStruct> precursorStruct = precursorStructs[level];
        vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, precursorStruct->nPoints);

        if(this->outputVariable==OutputVariable::Velocities)
        {
            fillArrayVelocities<<<grid.grid, grid.threads>>>(   precursorStruct->nPoints, precursorStruct->indicesD, 
                                                                precursorStruct->dataD, 
                                                                para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                para->getVelocityRatio());
            getLastCudaError("In PrecursorWriter::interact fillArrayVelocities execution failed");
        }
        else if(this->outputVariable==OutputVariable::Distributions)
        {
            fillArrayDistributions<<<grid.grid, grid.threads>>>(precursorStruct->nPoints, precursorStruct->indicesD, 
                                                                precursorStruct->dataD,
                                                                para->getParD(level)->distributions.f[0],
                                                                para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                para->getEvenOrOdd(level), para->getParD(level)->numberOfNodes);
            getLastCudaError("In PrecursorWriter::interact fillArrayDistributions execution failed");
        }

        cudaManager->cudaCopyPrecursorWriterOutputVariablesDtoH(this, level);

        precursorStruct->timestepsBuffered++;

        if(precursorStruct->timestepsBuffered >= precursorStruct->timestepsPerFile)
            this->write(para, level);
    }
}


void PrecursorWriter::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        if(getPrecursorStruct(level)->timestepsBuffered>0)
            write(para, level);

        cudaManager->cudaFreePrecursorWriter(this, level);
    }
}


void PrecursorWriter::write(Parameter* para, int level)
{
    SPtr<PrecursorStruct> precursorStruct = this->getPrecursorStruct(level);
    std::string fname = this->makeFileName(fileName, level, para->getMyProcessID(), precursorStruct->filesWritten) + getWriter()->getFileExtension();
    std::string wholeName = outputPath + "/" + fname;

    uint nPointsInPlane = precursorStruct->nPointsInPlane;

    int startTime = precursorStruct->filesWritten*precursorStruct->timestepsPerFile;

    // printf("points in plane %d, total timesteps %d, ntimesteps %d \n", nPointsInPlane, nTotalTimesteps, nTimesteps);

    UbTupleInt6 extent = makeUbTuple(   val<1>(precursorStruct->extent),    val<2>(precursorStruct->extent), 
                                        val<3>(precursorStruct->extent),    val<4>(precursorStruct->extent), 
                                        startTime,                          startTime+(int)precursorStruct->timestepsBuffered-1);

    UbTupleFloat3 origin = makeUbTuple( val<1>(precursorStruct->origin), val<2>(precursorStruct->origin), 0.f);

    std::vector<std::vector<double>> nodedata;
    
    for(uint quant = 0; quant < precursorStruct->nQuantities; quant++)
    {
        std::vector<double> doubleArr(nPointsInPlane*precursorStruct->timestepsBuffered, NAN);
        for( uint timestep=0; timestep<precursorStruct->timestepsBuffered; timestep++)
        {
            for (uint pos=0; pos < precursorStruct->nPoints; pos++)
            {
                int indexOnPlane = precursorStruct->indicesOnPlane[pos]+timestep*nPointsInPlane;
                doubleArr[indexOnPlane] = double(precursorStruct->dataH[lIndex(quant, pos, timestep, precursorStruct->nQuantities, precursorStruct->nPoints)]);
            }
        }
        nodedata.push_back(doubleArr);
    }

    // if(this->outputVariable==OutputVariable::Velocities)
    // {
    //     std::vector<double> vxDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         vyDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         vzDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN);


    //     for( uint timestep=0; timestep<precursorStruct->timestepsBuffered; timestep++)
    //     {
    //         // printf("offset %d npoints %d buf %d, max%d\n",timestep, precursorStruct->nPoints, precursorStruct->timestepsBuffered, precursorStruct->timestepsPerFile);
    //         for (uint pos = 0; pos < precursorStruct->nPoints; pos++)
    //         {
    //             int indexOnPlane = precursorStruct->indicesOnPlane[pos]+timestep*nPointsInPlane;
    //             // int idx = pos+timestep*precursorStruct->nPoints*3;
    //             // printf("timestep %i, pos %i, iOP %i \n", timestep, pos, indexOnPlane);
    //             // printf("vx %f, vy %f, vz%f nodedata x %f\n", vx[level][timestep][pos], vy[level][timestep][pos], vz[level][timestep][pos], vxDouble[indexOnPlane]);
    //             vxDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(0u, pos, timestep, 3, precursorStruct->nPoints)]);
    //             vyDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(1u, pos, timestep, 3, precursorStruct->nPoints)]);
    //             vzDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(2u, pos, timestep, 3, precursorStruct->nPoints)]);
    //         }
    //     }
    //     nodedata = {vxDouble, vyDouble, vzDouble};
    // }
    // else if(this->outputVariable==OutputVariable::Distributions)
    // {
    //     std::vector<double> fP00Double(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fPP0Double(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fPM0Double(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fP0PDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fP0MDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fPPPDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fPMPDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fPPMDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
    //                         fPMMDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN);

    //     for( uint timestep=0; timestep<precursorStruct->timestepsBuffered; timestep++)
    //     {
    //         for (uint pos = 0; pos < precursorStruct->nPoints; pos++)
    //         {
    //             int indexOnPlane = precursorStruct->indicesOnPlane[pos]+timestep*nPointsInPlane;
    //             // int idx = pos+timestep*precursorStruct->nPoints*9;
    //             // printf("timestep %i, pos %i, iOP %i \n", timestep, pos, indexOnPlane);
    //             // printf("f0 %f\n", double(precursorStruct->dataH[PrecP00*precursorStruct->nPoints+idx]));
    //             fP00Double[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecP00, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fPP0Double[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecPP0, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fPM0Double[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecPM0, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fP0PDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecP0P, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fP0MDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecP0M, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fPPPDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecPPP, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fPMPDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecPMP, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fPPMDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecPPM, pos, timestep, 9, precursorStruct->nPoints)]);
    //             fPMMDouble[indexOnPlane] = double(precursorStruct->dataH[lIndex(PrecPMM, pos, timestep, 9, precursorStruct->nPoints)]);
    //         }
    //     }
    //     nodedata = {fP00Double, fPP0Double, fPM0Double, fP0PDouble, fP0MDouble, fPPPDouble, fPMPDouble, fPPMDouble, fPMMDouble};
    // }

    precursorStruct->timestepsBuffered = 0;

    std::vector<std::vector<double>> celldata;
    getWriter()->writeData(wholeName, nodedatanames, celldatanames, nodedata, celldata, extent, origin, precursorStruct->spacing, extent);
    precursorStruct->filesWritten++;
}

std::string PrecursorWriter::makeFileName(std::string fileName, int level, int id, uint filesWritten)
{
    return fileName + "_lev_" + StringUtil::toString<int>(level)
                    + "_ID_" + StringUtil::toString<int>(id)
                    + "_File_" + StringUtil::toString<int>(filesWritten);
}