#include "CumulantAll4CompSP27.h"

#include "CumulantAll4CompSP27_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<CumulantAll4CompSP27> CumulantAll4CompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantAll4CompSP27>(new CumulantAll4CompSP27(para, level));
}

void CumulantAll4CompSP27::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_Cumulant_D3Q27All4 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        level,
        para->getForcesDev(),
        para->getQuadricLimitersDev(),
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_Cumulant_D3Q27All4 execution failed");
}

CumulantAll4CompSP27::CumulantAll4CompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}