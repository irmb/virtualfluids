#include "CumulantK18Comp.h"

#include "CumulantK18Comp_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<CumulantK18Comp> CumulantK18Comp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK18Comp>(new CumulantK18Comp(para, level));
}

void CumulantK18Comp::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_CumulantK18Comp <<< grid.grid, grid.threads >>>(
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
    getLastCudaError("LB_Kernel_CumulantK18Comp execution failed");
}

CumulantK18Comp::CumulantK18Comp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);
	myPreProcessorTypes.push_back(InitF3);

	
}