#include "CumulantK17Comp.h"

#include "Parameter/Parameter.h"
#include "CumulantK17Comp_Device.cuh"

std::shared_ptr<CumulantK17Comp> CumulantK17Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK17Comp>(new CumulantK17Comp(para,level));
}

void CumulantK17Comp::run()
{
    dim3 grid, threads;
    std::tie(grid, threads) =
        *calcGridDimensions(para->getParD(level)->size_Mat_SP, para->getParD(level)->numberofthreads);

	LB_Kernel_CumulantK17Comp <<< grid, threads >>>(para->getParD(level)->omega,
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

CumulantK17Comp::CumulantK17Comp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}