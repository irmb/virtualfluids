#include "CumulantK15SpongeComp.h"

#include "CumulantK15SpongeComp_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<CumulantK15SpongeComp> CumulantK15SpongeComp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK15SpongeComp>(new CumulantK15SpongeComp(para, level));
}

void CumulantK15SpongeComp::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_CumulantK15SpongeComp <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->coordinateX,
        para->getParD(level)->coordinateY,
        para->getParD(level)->coordinateZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_CumulantK15SpongeComp execution failed");
}

CumulantK15SpongeComp::CumulantK15SpongeComp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

CumulantK15SpongeComp::CumulantK15SpongeComp()
{
}
