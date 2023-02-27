#include "ADIncomp27.h"

#include "ADIncomp27_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<ADIncomp27> ADIncomp27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<ADIncomp27>(new ADIncomp27(para, level));
}

void ADIncomp27::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_AD_Incomp_27<<< grid.grid, grid.threads >>>(
        para->getParD(level)->diffusivity, 
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX, 
        para->getParD(level)->neighborY, 
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0], 
        para->getParD(level)->distributionsAD27.f[0], 
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_AD_Incomp_27 execution failed");
}

ADIncomp27::ADIncomp27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitIncompAD27);

}

ADIncomp27::ADIncomp27()
{
}
