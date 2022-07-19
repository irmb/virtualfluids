#include "CumulantK17CompChimRedesigned.h"

#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "CumulantK17CompChimRedesigned_Device.cuh"

#include <cuda.h>

std::shared_ptr<CumulantK17CompChimRedesigned> CumulantK17CompChimRedesigned::getNewInstance(std::shared_ptr<Parameter> para,
                                                                               int level)
{
    return std::shared_ptr<CumulantK17CompChimRedesigned>(new CumulantK17CompChimRedesigned(para, level));
}

void CumulantK17CompChimRedesigned::run()
{
    LB_Kernel_CumulantK17CompChimRedesigned <<< cudaGrid.grid, cudaGrid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        level,
        para->getForcesDev(),
        para->getQuadricLimitersDev(),
        para->getParD(level)->rho,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->isEvenTimestep,
        para->getParD(level)->fluidNodeIndices,
        para->getParD(level)->numberOfFluidNodes);
    getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

void CumulantK17CompChimRedesigned::runOnIndices(const unsigned int *indices, unsigned int size_indices, int streamIndex)
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager()->getStream(streamIndex);

    LB_Kernel_CumulantK17CompChimRedesigned<<< cudaGrid.grid, cudaGrid.threads, 0, stream>>>(
        para->getParD(level)->omega, 
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ, 
        para->getParD(level)->distributions.f[0], 
        para->getParD(level)->numberOfNodes, 
        level,
        para->getForcesDev(), 
        para->getQuadricLimitersDev(),
        para->getParD(level)->rho,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->isEvenTimestep,
        indices,
        size_indices);
    getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
    
}

CumulantK17CompChimRedesigned::CumulantK17CompChimRedesigned(std::shared_ptr<Parameter> para, int level): KernelImp(para, level)
{
    myPreProcessorTypes.push_back(InitCompSP27);
    myKernelGroup = BasicKernel;
    this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
}

