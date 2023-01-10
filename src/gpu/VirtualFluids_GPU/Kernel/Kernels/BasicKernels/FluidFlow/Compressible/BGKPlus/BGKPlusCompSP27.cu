#include "BGKPlusCompSP27.h"

#include "BGKPlusCompSP27_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<BGKPlusCompSP27> BGKPlusCompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<BGKPlusCompSP27>(new BGKPlusCompSP27(para, level));
}

void BGKPlusCompSP27::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_BGK_Plus_Comp_SP_27 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_BGK_Plus_Comp_SP_27 execution failed");
}

BGKPlusCompSP27::BGKPlusCompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

BGKPlusCompSP27::BGKPlusCompSP27()
{
}
