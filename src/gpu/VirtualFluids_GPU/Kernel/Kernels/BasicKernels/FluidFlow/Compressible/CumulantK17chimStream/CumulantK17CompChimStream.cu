#include "CumulantK17CompChimStream.h"

#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "CumulantK17CompChimStream_Device.cuh"

#include <cuda.h>

std::shared_ptr<CumulantK17CompChimStream> CumulantK17CompChimStream::getNewInstance(std::shared_ptr<Parameter> para,
                                                                               int level)
{
    return std::shared_ptr<CumulantK17CompChimStream>(new CumulantK17CompChimStream(para, level));
}

void CumulantK17CompChimStream::run()
{
    LB_Kernel_CumulantK17CompChimStream <<< cudaGrid.grid, cudaGrid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        level,
        para->getForcesDev(),
        para->getQuadricLimitersDev(),
        para->getParD(level)->isEvenTimestep,
        para->getParD(level)->fluidNodeIndices,
        para->getParD(level)->numberOfFluidNodes);
    getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

void CumulantK17CompChimStream::runOnIndices(const unsigned int *indices, unsigned int size_indices, int streamIndex)
{
    cudaStream_t stream = para->getStreamManager()->getStream(CudaStreamManager::StreamIndex::bulkStreamIndex);

    LB_Kernel_CumulantK17CompChimStream<<< cudaGrid.grid, cudaGrid.threads, 0, stream>>>(
        para->getParD(level)->omega, 
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ, 
        para->getParD(level)->distributions.f[0], 
        para->getParD(level)->numberOfNodes, 
        level,
        para->getForcesDev(), 
        para->getQuadricLimitersDev(), 
        para->getParD(level)->isEvenTimestep,
        indices,
        size_indices);
    getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
    
}

CumulantK17CompChimStream::CumulantK17CompChimStream(std::shared_ptr<Parameter> para, int level): KernelImp(para, level)
{
    myPreProcessorTypes.push_back(InitCompSP27);
    myKernelGroup = BasicKernel;
    this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
    this->kernelUsesFluidNodeIndices = true;
}

