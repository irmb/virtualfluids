#include "B15CompressibleNavierStokesBGKplus.h"

#include "B15CompressibleNavierStokesBGKplus_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<B15CompressibleNavierStokesBGKplus> B15CompressibleNavierStokesBGKplus::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<B15CompressibleNavierStokesBGKplus>(new B15CompressibleNavierStokesBGKplus(para, level));
}

void B15CompressibleNavierStokesBGKplus::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    B15CompressibleNavierStokesBGKplus_Device<<<grid.grid, grid.threads>>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("B15CompressibleNavierStokesBGKplus_Device execution failed");
}

B15CompressibleNavierStokesBGKplus::B15CompressibleNavierStokesBGKplus(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

B15CompressibleNavierStokesBGKplus::B15CompressibleNavierStokesBGKplus()
{
}
