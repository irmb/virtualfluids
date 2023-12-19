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
//! \file PrecursorWriter.cu
//! \ingroup PreCollisionInteractor
//! \author Henrik Asmuth, Henry Korb
//======================================================================================
#include "PrecursorWriter.h"
#include "basics/writer/WbWriterVtkXmlImageBinary.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuda_helper/CudaGrid.h>
#include "Utilities/KernelUtilities.h"

#include "StringUtilities/StringUtil.h"

#include "Parameter/Parameter.h"
#include "DataStructureInitializer/GridProvider.h"
#include "Cuda/CudaMemoryManager.h"
#include "Output/FilePartCalculator.h"

using namespace vf::lbm::dir;
using namespace vf::gpu;

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

__inline__ __host__ __device__ uint linearIdx(const uint component, const uint node, const uint timestep, const uint numberOfComponents, const uint numberOfNodes)
{
    return node+numberOfNodes*(component+numberOfComponents*timestep);
}

__inline__ __host__ __device__ uint linearIdx(const uint component, const uint node, const uint numberOfNodes)
{
    return node+component*numberOfNodes;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void fillArrayVelocities(const uint numberOfPrecursorNodes, 
                                    uint* indices, 
                                    real *precursorData,
                                    real *vx,
                                    real *vy,
                                    real *vz,
                                    real velocityRatio)


{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if(nodeIndex>=numberOfPrecursorNodes) return;

    precursorData[linearIdx(0u, nodeIndex, numberOfPrecursorNodes)] = vx[indices[nodeIndex]]*velocityRatio;
    precursorData[linearIdx(1u, nodeIndex, numberOfPrecursorNodes)] = vy[indices[nodeIndex]]*velocityRatio;
    precursorData[linearIdx(2u, nodeIndex, numberOfPrecursorNodes)] = vz[indices[nodeIndex]]*velocityRatio;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void fillArrayDistributions( uint numberOfPrecursorNodes, 
                                        uint* indices, 
                                        real* precursorData,
                                        real* distributions,
                                        uint* neighborX, uint* neighborY, uint* neighborZ,
                                        bool isEvenTimestep,
                                        unsigned long numberOfLBnodes)
{
    ////////////////////////////////////////////////////////////////////////////////
    //! - Get node index coordinates from threadIdx, blockIdx, blockDim and gridDim.
    //!
    const unsigned nodeIndex = vf::gpu::getNodeIndex();

    if(nodeIndex>=numberOfPrecursorNodes) return;

    Distributions27 dist;
    getPointersToDistributions(dist, distributions, numberOfLBnodes, isEvenTimestep);
    
    ////////////////////////////////////////////////////////////////////////////////
    // ! - Set neighbor indices (necessary for indirect addressing)
    uint k_000 = indices[nodeIndex];
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
    precursorData[linearIdx(PrecP00, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dP00])[k_000];
    precursorData[linearIdx(PrecPP0, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPP0])[k_000];
    precursorData[linearIdx(PrecPM0, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPM0])[k_0M0];
    precursorData[linearIdx(PrecP0P, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dP0P])[k_000];
    precursorData[linearIdx(PrecP0M, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dP0M])[k_00M];
    precursorData[linearIdx(PrecPPP, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPPP])[k_000];
    precursorData[linearIdx(PrecPMP, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPMP])[k_0M0];
    precursorData[linearIdx(PrecPPM, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPPM])[k_00M];
    precursorData[linearIdx(PrecPMM, nodeIndex, numberOfPrecursorNodes)] = (dist.f[dPMM])[k_0MM];
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrecursorWriter::init()
{
    VF_LOG_INFO("PrecursorWriter: Start initializing...");
    VF_LOG_INFO("Writing yz-planes at x={}m every {}. timestep, starting at t={}", this->xPos, this->tSave, this->tStartOut);

    precursorStructs.resize(para->getMaxLevel()+1);
    for(int level=0; level<=para->getMaxLevel(); level++)
    {

        real dx = abs(para->getParH(level)->coordinateX[1]-para->getParH(level)->coordinateX[para->getParH(level)->neighborX[1]]);

        real lowestY, lowestZ, highestY, highestZ;

        lowestY = para->getParH(level)->coordinateY[para->getParH(level)->numberOfNodes-1];
        highestY = para->getParH(level)->coordinateY[1];        
        
        lowestZ = para->getParH(level)->coordinateZ[para->getParH(level)->numberOfNodes-1];
        highestZ = para->getParH(level)->coordinateZ[1];

        std::vector<uint> indicesOnGrid;
        std::vector<int> indicesOnPlane;
        std::vector<real> coordY, coordZ;

        for(size_t pos = 1; pos < para->getParH(level)->numberOfNodes; pos++ )
        {
            real pointCoordX = para->getParH(level)->coordinateX[pos];
            real pointCoordY = para->getParH(level)->coordinateY[pos];
            real pointCoordZ = para->getParH(level)->coordinateZ[pos];
            if( para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID &&
                pointCoordX < (dx+xPos) && pointCoordX >= xPos       &&
                pointCoordY<=yMax && pointCoordY>=yMin               && 
                pointCoordZ<=zMax && pointCoordZ>=zMin)
            {
                highestY = max(highestY, pointCoordY);
                highestZ = max(highestZ, pointCoordZ);

                lowestY = min(lowestY, pointCoordY);
                lowestZ = min(lowestZ, pointCoordZ);
                indicesOnGrid.push_back((uint)pos);    
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
        precursorStructs[level]->numberOfPointsInBC = (uint)indicesOnGrid.size();
        precursorStructs[level]->indicesOnPlane = (int*) malloc(precursorStructs[level]->numberOfPointsInBC*sizeof(int));
        precursorStructs[level]->spacing = makeUbTuple(dx, dx, tSave*para->getTimeRatio()*pow(2,-level));
        precursorStructs[level]->origin = makeUbTuple(lowestY, lowestZ);
        precursorStructs[level]->extent = makeUbTuple(0, ny-1, 0, nz-1);
        precursorStructs[level]->numberOfPointsInData = ny*nz;
        precursorStructs[level]->numberOfTimestepsPerFile = min(FilePartCalculator::limitOfNodesForVTK/(ny*nz), maxtimestepsPerFile);
        precursorStructs[level]->numberOfFilesWritten = 0;
        precursorStructs[level]->numberOfTimestepsBuffered = 0;
        
        switch (outputVariable)
        {
        case OutputVariable::Velocities:
            precursorStructs[level]->numberOfQuantities = 3;
            break;
        case OutputVariable::Distributions:
            precursorStructs[level]->numberOfQuantities = 9;
            break;
        
        default:
            break;
        }

        cudaMemoryManager->cudaAllocPrecursorWriter(this, level);
    
        std::copy(indicesOnGrid.begin(), indicesOnGrid.end(), precursorStructs[level]->indicesH);
        std::copy(indicesOnPlane.begin(), indicesOnPlane.end(), precursorStructs[level]->indicesOnPlane);

        cudaMemoryManager->cudaCopyPrecursorWriterIndicesHtoD(this, level);

        VF_LOG_INFO("Found {} points in precursor plane on level {}", precursorStructs[level]->numberOfPointsInBC, level);
    }
    VF_LOG_INFO("PrecursorWriter: Done initializing.");
}


void PrecursorWriter::interact(int level, uint t)
{
    uint t_level         = para->getTimeStep(level, t, true);
    uint tStartOut_level = tStartOut*pow(2, level);

    if(t_level>tStartOut_level && ((t_level-tStartOut_level) % tSave)==0)
    {
        vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, precursorStructs[level]->numberOfPointsInBC);

        if(this->outputVariable==OutputVariable::Velocities)
        {
            fillArrayVelocities<<<grid.grid, grid.threads>>>(   precursorStructs[level]->numberOfPointsInBC, precursorStructs[level]->indicesD, 
                                                                precursorStructs[level]->bufferD, 
                                                                para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                para->getVelocityRatio());
            getLastCudaError("In PrecursorWriter::interact fillArrayVelocities execution failed");
        }
        else if(this->outputVariable==OutputVariable::Distributions)
        {
            fillArrayDistributions<<<grid.grid, grid.threads>>>(precursorStructs[level]->numberOfPointsInBC, precursorStructs[level]->indicesD, 
                                                                precursorStructs[level]->bufferD,
                                                                para->getParD(level)->distributions.f[0],
                                                                para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                para->getEvenOrOdd(level), para->getParD(level)->numberOfNodes);
            getLastCudaError("In PrecursorWriter::interact fillArrayDistributions execution failed");
        }
        cudaMemoryManager->cudaCopyPrecursorWriterOutputVariablesDtoH(this, level);

        // switch device buffer and data pointer so precursor data is gathered in buffer and copied from bufferD to bufferH
        real *tmp = precursorStructs[level]->bufferD;
        precursorStructs[level]->bufferD = precursorStructs[level]->dataD;
        precursorStructs[level]->dataD = tmp;

        precursorStructs[level]->numberOfTimestepsBuffered++;

        if(precursorStructs[level]->numberOfTimestepsBuffered >= precursorStructs[level]->numberOfTimestepsPerFile || t == para->getTimestepEnd())
        {
        // switch host buffer and data pointer so precursor data is copied in buffer and written from data

            tmp = precursorStructs[level]->bufferH;
            precursorStructs[level]->bufferH = precursorStructs[level]->dataH;
            precursorStructs[level]->dataH = tmp;

            writeFuture.wait();
            writeFuture = std::async(std::launch::async, [this](uint level, uint timesteps){ this->write(level, timesteps); }, level, precursorStructs[level]->numberOfTimestepsBuffered);
            precursorStructs[level]->numberOfTimestepsBuffered = 0;
        }
    }
}


PrecursorWriter::~PrecursorWriter()
{
    writeFuture.wait();
    for(int level=0; level<=para->getMaxLevel(); level++)
    {
        if(getPrecursorStruct(level)->numberOfTimestepsBuffered>0)
            write(level, getPrecursorStruct(level)->numberOfTimestepsBuffered);

        cudaMemoryManager->cudaFreePrecursorWriter(this, level);
    }
}


void PrecursorWriter::write(int level, uint numberOfTimestepsBuffered)
{
    std::string fname = this->makeFileName(fileName, level, para->getMyProcessID(), precursorStructs[level]->numberOfFilesWritten) + getWriter()->getFileExtension();
    std::string wholeName = outputPath + "/" + fname;

    uint numberOfPointsInData = precursorStructs[level]->numberOfPointsInData;

    int startTime = precursorStructs[level]->numberOfFilesWritten*precursorStructs[level]->numberOfTimestepsPerFile;

    UbTupleInt6 extent = makeUbTuple(   val<1>(precursorStructs[level]->extent),    val<2>(precursorStructs[level]->extent), 
                                        val<3>(precursorStructs[level]->extent),    val<4>(precursorStructs[level]->extent), 
                                        startTime,                          startTime+(int)numberOfTimestepsBuffered-1);

    UbTupleFloat3 origin = makeUbTuple( val<1>(precursorStructs[level]->origin), val<2>(precursorStructs[level]->origin), 0.f);

    std::vector<std::vector<double>> nodedata;
    
    for(uint quant=0; quant<precursorStructs[level]->numberOfQuantities; quant++)
    {
        std::vector<double> doubleArr(numberOfPointsInData*numberOfTimestepsBuffered, NAN);
        for( uint timestep=0; timestep<numberOfTimestepsBuffered; timestep++)
        {
            for (uint pos=0; pos < precursorStructs[level]->numberOfPointsInBC; pos++)
            {
                int indexOnPlane = precursorStructs[level]->indicesOnPlane[pos]+timestep*numberOfPointsInData;
                doubleArr[indexOnPlane] = double(precursorStructs[level]->dataH[linearIdx(quant, pos, timestep, precursorStructs[level]->numberOfQuantities, precursorStructs[level]->numberOfPointsInBC)]);
            }
        }
        nodedata.push_back(doubleArr);
    }

    std::vector<std::vector<double>> celldata;
    getWriter()->writeData(wholeName, nodedatanames, celldatanames, nodedata, celldata, extent, origin, precursorStructs[level]->spacing, extent, this->writePrecision);
    precursorStructs[level]->numberOfFilesWritten++;
}

std::string PrecursorWriter::makeFileName(std::string fileName, int level, int id, uint numberOfFilesWritten)
{
    return fileName + "_lev_" + StringUtil::toString<int>(level)
                    + "_ID_" + StringUtil::toString<int>(id)
                    + "_File_" + StringUtil::toString<int>(numberOfFilesWritten);
}

void PrecursorWriter::getTaggedFluidNodes(GridProvider* gridProvider)
{
    for(uint level=0; level<(uint)para->getMaxLevel(); level++)
    {
        if(outputVariable==OutputVariable::Velocities)
        {
            std::vector<uint> indices(precursorStructs[level]->indicesH, precursorStructs[level]->indicesH+precursorStructs[level]->numberOfPointsInBC);
            gridProvider->tagFluidNodeIndices(indices, CollisionTemplate::WriteMacroVars, level);
        }
    }
}