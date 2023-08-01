#include "ADIncomp7.h"

#include "ADIncomp7_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<ADIncomp7> ADIncomp7::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<ADIncomp7>(new ADIncomp7(para, level));
}

void ADIncomp7::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_AD_Incomp_7<<< grid.grid, grid.threads >>>(
        para->getParD(level)->diffusivity, 
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY, 
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->distributionsAD7.f[0], 
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_AD_Incomp_7 execution failed");
}

ADIncomp7::ADIncomp7(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitIncompAD7);

}

ADIncomp7::ADIncomp7()
{
}
