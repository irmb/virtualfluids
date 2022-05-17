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
                            real *vz)
{
    const uint x = threadIdx.x; 
    const uint y = blockIdx.x;
    const uint z = blockIdx.y;

    const uint nx = blockDim.x;
    const uint ny = gridDim.x;

    const uint node = nx*(ny*z + y) + x;

    if(node>=nNodes) return;

    precursorVx[node] = vx[indices[node]];
    precursorVy[node] = vy[indices[node]];
    precursorVz[node] = vz[indices[node]];
}

void PrecursorWriter::init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)
{
    precursorStructs.resize(para->getMaxLevel()+1);
    vx.resize(para->getMaxLevel()+1);
    vy.resize(para->getMaxLevel()+1);
    vz.resize(para->getMaxLevel()+1);
    for(int level=0; level<=para->getMaxLevel(); level++)
    {

        real dx = abs(para->getParH(level)->coordX_SP[1]-para->getParH(level)->coordX_SP[para->getParH(level)->neighborX_SP[1]]);
        int maxPoints = (int((yMax-yMin)/dx)+1)* (int((zMax-zMin)/dx)+1);

        real lowestY, lowestZ, highestY, highestZ;

        lowestY = para->getParH(level)->coordY_SP[para->getParH(level)->size_Mat_SP-1];
        highestY = para->getParH(level)->coordY_SP[1];        
        
        lowestZ = para->getParH(level)->coordZ_SP[para->getParH(level)->size_Mat_SP-1];
        highestZ = para->getParH(level)->coordZ_SP[1];

        std::vector<uint> indicesOnGrid;
        std::vector<int> indicesOnPlane;
        std::vector<real> coordY, coordZ;

        for(uint j=1; j<para->getParH(level)->size_Mat_SP; j++ )
        {
            real pointCoordX = para->getParH(level)->coordX_SP[j];
            real pointCoordY = para->getParH(level)->coordY_SP[j];
            real pointCoordZ = para->getParH(level)->coordZ_SP[j];
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
        int npoints = indicesOnGrid.size();

        precursorStructs[level] = SPtr<PrecursorStruct>(new PrecursorStruct);
        precursorStructs[level]->nPoints = indicesOnGrid.size();
        precursorStructs[level]->indicesOnPlane = (int*) malloc(precursorStructs[level]->nPoints*sizeof(int));
        cudaManager->cudaAllocPrecursorWriter(this, level);
    
        std::copy(indicesOnGrid.begin(), indicesOnGrid.end(), precursorStructs[level]->indicesH);
        std::copy(indicesOnPlane.begin(), indicesOnPlane.end(), precursorStructs[level]->indicesOnPlane);
        precursorStructs[level]->spacing = makeUbTuple(dx, dx, tSave*para->getTimeRatio());
        precursorStructs[level]->origin = makeUbTuple(lowestY, lowestZ);
        precursorStructs[level]->extent = makeUbTuple(0, ny-1, 0, nz-1);
        precursorStructs[level]->nPointsInPlane = ny*nz;
        precursorStructs[level]->timestepsPerFile = min(para->getlimitOfNodesForVTK()/(ny*nz), maxtimestepsPerFile);
        precursorStructs[level]->filesWritten = 0;
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
                                                para->getParD(level)->vx_SP, para->getParD(level)->vy_SP, para->getParD(level)->vz_SP);

        cudaManager->cudaCopyPrecursorWriterVelocitiesDtoH(this, level);
        
        std::vector<real> new_vx(precursorStruct->vxH,precursorStruct->vxH+precursorStruct->nPoints), 
                          new_vy(precursorStruct->vyH,precursorStruct->vyH+precursorStruct->nPoints), 
                          new_vz(precursorStruct->vzH,precursorStruct->vzH+precursorStruct->nPoints);

        vx[level].push_back(new_vx);
        vy[level].push_back(new_vy);
        vz[level].push_back(new_vz);


        if(vx[level].size() > precursorStruct->timestepsPerFile)
            this->write(para, level);

    }
}


void PrecursorWriter::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        if(vx[level].size()>0)
            write(para, level);

        cudaManager->cudaFreePrecursorWriter(this, level);
    }
}


void PrecursorWriter::write(Parameter* para, int level)
{
    SPtr<PrecursorStruct> precursorStruct = this->getPrecursorStruct(level);
    std::string fname = this->makeFileName(level, para->getMyID(), precursorStruct->filesWritten) + getWriter()->getFileExtension();
    std::string wholeName = outputPath + "/" + fname;

    uint nPointsInPlane = precursorStruct->nPointsInPlane;

    int nTimesteps = vx[level].size();
    int startTime = precursorStruct->filesWritten*precursorStruct->timestepsPerFile;

    // printf("points in plane %d, total timesteps %d, ntimesteps %d \n", nPointsInPlane, nTotalTimesteps, nTimesteps);
    std::vector<double> vxDouble(nPointsInPlane*nTimesteps, NAN), vyDouble(nPointsInPlane*nTimesteps, NAN), vzDouble(nPointsInPlane*nTimesteps, NAN);

    UbTupleInt6 extent = makeUbTuple(   val<1>(precursorStruct->extent),    val<2>(precursorStruct->extent), 
                                        val<3>(precursorStruct->extent),    val<4>(precursorStruct->extent), 
                                        startTime,                          startTime+nTimesteps-1);

    UbTupleFloat3 origin = makeUbTuple( val<1>(precursorStruct->origin), val<1>(precursorStruct->origin), 0.f);

    real coeff = para->getVelocityRatio();
        
    for( uint timestep=0; timestep<nTimesteps; timestep++)
    {
        for (uint pos = 0; pos < this->getPrecursorStruct(level)->nPoints; pos++)
        {

            int indexOnPlane = precursorStruct->indicesOnPlane[pos]+timestep*nPointsInPlane;
            // printf("timestep %i, pos %i, iOP %i \n", timestep, pos, indexOnPlane);
            // printf("vx %f, vy %f, vz%f nodedata x %f\n", vx[level][timestep][pos]*coeff, vy[level][timestep][pos]*coeff, vz[level][timestep][pos]*coeff, vxDouble[indexOnPlane]);
            vxDouble[indexOnPlane] = double(vx[level][timestep][pos]*coeff);
            vyDouble[indexOnPlane] = double(vy[level][timestep][pos]*coeff);
            vzDouble[indexOnPlane] = double(vz[level][timestep][pos]*coeff);
        }
    }

    vx[level].clear();
    vy[level].clear();
    vz[level].clear();

    std::vector<std::vector<double>> nodedata = {vxDouble, vyDouble, vzDouble};

    std::vector<std::vector<double>> celldata;
    getWriter()->writeData(wholeName, nodedatanames, celldatanames, nodedata, celldata, extent, origin, precursorStruct->spacing, extent);
    precursorStruct->filesWritten++;
}

std::string PrecursorWriter::makeFileName(int level, int id, uint filesWritten)
{
    return fileName + "_lev_" + StringUtil::toString<int>(level)
                    + "_ID_" + StringUtil::toString<int>(id)
                    + "_File_" + StringUtil::toString<int>(filesWritten);
}