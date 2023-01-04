#include "CumulantK20Comp.h"

#include "CumulantK20Comp_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<CumulantK20Comp> CumulantK20Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK20Comp>(new CumulantK20Comp(para, level));
}

void CumulantK20Comp::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_CumulantK20Comp <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->g6.g[0],
        para->getParD(level)->numberOfNodes,
        level,
        para->getForcesDev(),
        para->getQuadricLimitersDev(),
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_CumulantK20Comp execution failed");
}

CumulantK20Comp::CumulantK20Comp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);
	myPreProcessorTypes.push_back(InitF3);

	myKernelGroup = F3Kernel;
}