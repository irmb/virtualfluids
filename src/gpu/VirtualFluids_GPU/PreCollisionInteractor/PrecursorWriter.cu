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
                printf("idx %d, idy %d, idz %d, ny %d, nz %d\n", idx, idxY, idxZ, ny, nz);
        }
        int npoints = indicesOnGrid.size();

        precursorStructs[level] = SPtr<PrecursorStruct>(new PrecursorStruct);
        precursorStructs[level]->nPoints = indicesOnGrid.size();
        precursorStructs[level]->indicesOnPlane = (int*) malloc(precursorStructs[level]->nPoints*sizeof(int));
        cudaManager->cudaAllocPrecursorWriter(this, level);
    
        std::copy(indicesOnGrid.begin(), indicesOnGrid.end(), precursorStructs[level]->indicesH);
        std::copy(indicesOnPlane.begin(), indicesOnPlane.end(), precursorStructs[level]->indicesOnPlane);
        precursorStructs[level]->spacing = makeUbTuple(dx, dx, tSave*para->getTimeRatio());
        precursorStructs[level]->origin = makeUbTuple(lowestY, lowestZ, 0);
        precursorStructs[level]->extent = makeUbTuple(0, ny-1, 0, nz-1);
        precursorStructs[level]->ny = ny;
        precursorStructs[level]->nz = nz;
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

        if(t>this->tStartOut? ((t-tStartOut) % this->tWrite) == 0 : false)
        {
            this->write(para, level);
        }
    }
}


void PrecursorWriter::free(Parameter* para, CudaMemoryManager* cudaManager)
{
    if(vx[0].size()>0)
    {
        for(int level=0; level<=para->getMaxLevel(); level++)
            write(para, level);
    }
    for(int level=0; level<=para->getMaxLevel(); level++)
        cudaManager->cudaFreePrecursorWriter(this, level);
}


void PrecursorWriter::write(Parameter* para, int level)
{
    const uint timestepsPerFile = para->getlimitOfNodesForVTK() / this->getPrecursorStruct(level)->nPoints;

    const uint numberOfParts = this->vx[level].size() / timestepsPerFile + 1;

    for (uint part = 0; part < numberOfParts; part++)
	{
        this->writeGridFile(para, level, part, timestepsPerFile);
    }

    if(level == 0) this->writeParallelFile(para, level, numberOfParts);

    nFilesWritten++;
    pieceSources.clear();
    pieceExtents.clear();
    vx[level].clear();
    vy[level].clear();
    vz[level].clear();
}


void PrecursorWriter::writeGridFile(Parameter* para, int level, uint part, int numberOfTimestepsPerPart)
{
    std::string fname = this->makeGridFileName(level, para->getMyID(), part) + getWriter()->getFileExtension();
    std::string wholeName = outputPath + "/" + fname;

    SPtr<PrecursorStruct> precursorStruct = this->getPrecursorStruct(level);
    uint nPointsInPlane = precursorStruct->ny*precursorStruct->nz;

    int nTotalTimesteps = vx[level].size();
    int startTime = numberOfTimestepsPerPart*part;
    int nTimesteps = min(numberOfTimestepsPerPart, nTotalTimesteps-startTime);

    printf("points in plane %d, total timesteps %d, ntimesteps %d \n", nPointsInPlane, nTotalTimesteps, nTimesteps);
    std::vector<double> vxDouble(nPointsInPlane*nTimesteps, NAN), vyDouble(nPointsInPlane*nTimesteps, NAN), vzDouble(nPointsInPlane*nTimesteps, NAN);


    UbTupleInt6 wholeExtent = makeUbTuple(  val<1>(precursorStruct->extent),    val<2>(precursorStruct->extent), 
                                            val<3>(precursorStruct->extent),    val<4>(precursorStruct->extent), 
                                            0,                                  nTotalTimesteps-1);

    UbTupleInt6 extent = makeUbTuple(   val<1>(precursorStruct->extent),    val<2>(precursorStruct->extent), 
                                        val<3>(precursorStruct->extent),    val<4>(precursorStruct->extent), 
                                        0,                                  nTimesteps-1);

    real coeff = para->getVelocityRatio();
        
    for( uint timestep=startTime; timestep<nTimesteps; timestep++)
    {
        for (uint pos = 0; pos < this->getPrecursorStruct(level)->nPoints; pos++)
        {

            int indexOnPlane = precursorStruct->indicesOnPlane[pos]+(timestep-startTime)*nPointsInPlane;
            // printf("timestep %i, pos %i, iOP %i \n", timestep, pos, indexOnPlane);
            // printf("vx %f, vy %f, vz%f nodedata x %f\n", vx[level][timestep][pos]*coeff, vy[level][timestep][pos]*coeff, vz[level][timestep][pos]*coeff, vxDouble[indexOnPlane]);
            vxDouble[indexOnPlane] = double(vx[level][timestep][pos]*coeff);
            vyDouble[indexOnPlane] = double(vy[level][timestep][pos]*coeff);
            vzDouble[indexOnPlane] = double(vz[level][timestep][pos]*coeff);
        }
    }

    std::vector<std::vector<double>> nodedata;

    nodedata.push_back(vxDouble);
    nodedata.push_back(vyDouble);
    nodedata.push_back(vzDouble);

    std::vector<std::vector<double>> celldata;
    pieceExtents.push_back(extent);
    pieceSources.push_back(fname);
    getWriter()->writeData(wholeName, nodedatanames, celldatanames, nodedata, celldata, wholeExtent, precursorStruct->origin, precursorStruct->spacing, extent);
}


void PrecursorWriter::writeParallelFile(Parameter* para, int level, uint parts)
{
    std::string fname = makeParallelFileName(level, para->getMyID());
    std::string wholeName = outputPath + "/" + fname;

    SPtr<PrecursorStruct> precursorStruct = this->getPrecursorStruct(level);


    UbTupleInt6 wholeExtent = makeUbTuple(val<1>(precursorStruct->extent), val<2>(precursorStruct->extent), 
                                            val<3>(precursorStruct->extent), val<4>(precursorStruct->extent), 
                                            val<5>(pieceExtents.front()),    val<6>(pieceExtents.back()));

    getWriter()->writeParallelFile(wholeName, wholeExtent, this->getPrecursorStruct(level)->origin, this->getPrecursorStruct(level)->spacing, pieceSources, pieceExtents, nodedatanames, celldatanames);
}

std::string PrecursorWriter::makeGridFileName(int level, int id, uint part)
{
    return fileName + "_lev_" + StringUtil::toString<int>(level)
                    + "_ID_" + StringUtil::toString<int>(id)
                    + "_Part_" + StringUtil::toString<int>(part) 
                    + "_File_" + StringUtil::toString<int>(nFilesWritten);
}
std::string PrecursorWriter::makeParallelFileName(int level, int id)
{
    return fileName + "_lev_" + StringUtil::toString<int>(level)
                    + "_ID_" + StringUtil::toString<int>(id)
                    + "_File_" + StringUtil::toString<int>(nFilesWritten);
}