#include "CascadeIncompSP27.h"

#include "CascadeIncompSP27_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<CascadeIncompSP27> CascadeIncompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CascadeIncompSP27>(new CascadeIncompSP27(para, level));
}

void CascadeIncompSP27::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_Cascade_Incomp_SP_27 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_Cascade_Incomp_SP_27 execution failed");
}

CascadeIncompSP27::CascadeIncompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	myKernelGroup = BasicKernel;
}

CascadeIncompSP27::CascadeIncompSP27()
{
}
