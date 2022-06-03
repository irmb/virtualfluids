#include "CumulantK17Comp.h"

#include "Parameter/Parameter.h"
#include "CumulantK17Comp_Device.cuh"
#include "cuda/CudaGrid.h"

std::shared_ptr<CumulantK17Comp> CumulantK17Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK17Comp>(new CumulantK17Comp(para,level));
}

void CumulantK17Comp::run()
{
	LB_Kernel_CumulantK17Comp <<< cudaGrid.grid, cudaGrid.threads >>>(para->getParD(level)->omega,
													para->getParD(level)->geoSP,
													para->getParD(level)->neighborX_SP,
													para->getParD(level)->neighborY_SP,
													para->getParD(level)->neighborZ_SP,
													para->getParD(level)->d0SP.f[0],
													para->getParD(level)->size_Mat_SP,
													level,
													para->getForcesDev(),
                                                    para->getQuadricLimitersDev(),
													para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_CumulantK17Comp execution failed");
}

CumulantK17Comp::CumulantK17Comp(std::shared_ptr<Parameter> para, int level): KernelImp(para, level)
{
	myPreProcessorTypes.push_back(InitCompSP27);
	myKernelGroup = BasicKernel;
	this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->size_Mat_SP);
}