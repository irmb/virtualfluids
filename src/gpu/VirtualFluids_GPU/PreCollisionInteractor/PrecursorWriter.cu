#include "PrecursorWriter.h"
#include "basics/writer/WbWriterVtkXmlImageBinary.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuda/CudaGrid.h>

#include <Core/StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"
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

__global__ void fillArray(uint nNodes, uint* indices, 
                            real *precursorVx,
                            real *precursorVy, 
                            real *precursorVz,  
                            real *vx,
                            real *vy,
                            real *vz,
                            real velocityRatio)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

    if(node>=nNodes) return;

    precursorVx[node] = vx[indices[node]]*velocityRatio;
    precursorVy[node] = vy[indices[node]]*velocityRatio;
    precursorVz[node] = vz[indices[node]]*velocityRatio;
}

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
        printf("ny %d nz %d \n", ny, nz);
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

        printf("points %zu points on plane %zu \n",  indicesOnGrid.size(),  indicesOnPlane.size());

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

        fillArray<<<grid.grid, grid.threads>>>(precursorStruct->nPoints, precursorStruct->indicesD, 
                                                precursorStruct->vxD, precursorStruct->vyD, precursorStruct->vzD, 
                                                para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                para->getVelocityRatio());

        cudaManager->cudaCopyPrecursorWriterVelocitiesDtoH(this, level);
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
    std::vector<double> vxDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
                        vyDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN), 
                        vzDouble(nPointsInPlane*precursorStruct->timestepsBuffered, NAN);

    UbTupleInt6 extent = makeUbTuple(   val<1>(precursorStruct->extent),    val<2>(precursorStruct->extent), 
                                        val<3>(precursorStruct->extent),    val<4>(precursorStruct->extent), 
                                        startTime,                          startTime+(int)precursorStruct->timestepsBuffered-1);

    UbTupleFloat3 origin = makeUbTuple( val<1>(precursorStruct->origin), val<2>(precursorStruct->origin), 0.f);
    for( uint timestep=0; timestep<precursorStruct->timestepsBuffered; timestep++)
    {
        // printf("offset %d npoints %d buf %d, max%d\n",timestep, precursorStruct->nPoints, precursorStruct->timestepsBuffered, precursorStruct->timestepsPerFile);
        for (uint pos = 0; pos < precursorStruct->nPoints; pos++)
        {
            int indexOnPlane = precursorStruct->indicesOnPlane[pos]+timestep*nPointsInPlane;
            int idx = pos+timestep*precursorStruct->nPoints;
            // printf("timestep %i, pos %i, iOP %i \n", timestep, pos, indexOnPlane);
            // printf("vx %f, vy %f, vz%f nodedata x %f\n", vx[level][timestep][pos], vy[level][timestep][pos], vz[level][timestep][pos], vxDouble[indexOnPlane]);
            vxDouble[indexOnPlane] = double(precursorStruct->vxH[idx]);
            vyDouble[indexOnPlane] = double(precursorStruct->vyH[idx]);
            vzDouble[indexOnPlane] = double(precursorStruct->vzH[idx]);
        }
    }

    precursorStruct->timestepsBuffered = 0;

    std::vector<std::vector<double>> nodedata = {vxDouble, vyDouble, vzDouble};

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