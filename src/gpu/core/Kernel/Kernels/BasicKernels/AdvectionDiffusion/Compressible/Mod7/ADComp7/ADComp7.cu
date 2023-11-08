#include "ADComp7.h"

#include "ADComp7_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<ADComp7> ADComp7::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<ADComp7>(new ADComp7(para, level));
}

void ADComp7::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_AD_Comp_7<<< grid.grid, grid.threads >>>(
        para->getParD(level)->diffusivity,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0], 
        para->getParD(level)->distributionsAD7.f[0], 
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_AD_Comp_7 execution failed");
}

ADComp7::ADComp7(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompAD7);

}

ADComp7::ADComp7()
{
}
