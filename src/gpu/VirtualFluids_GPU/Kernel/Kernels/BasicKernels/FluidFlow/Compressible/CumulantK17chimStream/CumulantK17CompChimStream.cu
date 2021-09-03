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
    dim3 grid, threads;
    std::tie(grid, threads) =
        *calcGridDimensions(para->getParD(level)->numberOfFluidNodes, para->getParD(level)->numberofthreads);

	LB_Kernel_CumulantK17CompChimStream <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->neighborX_SP,
		para->getParD(level)->neighborY_SP,
		para->getParD(level)->neighborZ_SP,
		para->getParD(level)->d0SP.f[0],
		para->getParD(level)->size_Mat_SP,
		level,
		para->getForcesDev(),
        para->getQuadricLimitersDev(),
		para->getParD(level)->evenOrOdd,
        para->getParD(level)->fluidNodeIndices,
		para->getParD(level)->numberOfFluidNodes);
	getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

void CumulantK17CompChimStream::runOnIndices(const unsigned int *indices, unsigned int size_indices, int streamIndex)
{
    dim3 grid, threads;
    std::tie(grid, threads) =
        *calcGridDimensions(para->getParD(level)->numberOfFluidNodes, para->getParD(level)->numberofthreads);

    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager()->getStream(streamIndex);

    LB_Kernel_CumulantK17CompChimStream<<<grid, threads, 0, stream>>>(
        para->getParD(level)->omega, 
	    para->getParD(level)->neighborX_SP, 
	    para->getParD(level)->neighborY_SP,
        para->getParD(level)->neighborZ_SP, 
	    para->getParD(level)->d0SP.f[0], 
	    para->getParD(level)->size_Mat_SP, 
	    level,
        para->getForcesDev(), 
	    para->getQuadricLimitersDev(), 
	    para->getParD(level)->evenOrOdd,
        indices,
	    size_indices);
    getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
    
}

CumulantK17CompChimStream::CumulantK17CompChimStream(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

