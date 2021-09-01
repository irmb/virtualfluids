#include "CumulantK17CompChimSparse.h"

#include "Parameter/Parameter.h"
#include "CumulantK17CompChimSparse_Device.cuh"

#include <cuda.h>

std::shared_ptr<CumulantK17CompChimSparse> CumulantK17CompChimSparse::getNewInstance(std::shared_ptr<Parameter> para,
                                                                               int level)
{
    return std::shared_ptr<CumulantK17CompChimSparse>(new CumulantK17CompChimSparse(para, level));
}

void CumulantK17CompChimSparse::run()
{
    dim3 grid, threads;
    std::tie(grid, threads) =
        *calcGridDimensions(para->getParD(level)->numberOfFluidNodes, para->getParD(level)->numberofthreads);

	LB_Kernel_CumulantK17CompChimSparse <<< grid, threads >>>(
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

void CumulantK17CompChimSparse::runOnIndices(const unsigned int *indices, unsigned int size_indices, int streamIndex)
{
    dim3 grid, threads;
    std::tie(grid, threads) =
        *calcGridDimensions(para->getParD(level)->numberOfFluidNodes, para->getParD(level)->numberofthreads);

    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager().getStream(streamIndex);

    LB_Kernel_CumulantK17CompChimSparse<<<grid, threads, 0, stream>>>(
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

CumulantK17CompChimSparse::CumulantK17CompChimSparse(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

