#include "M02IncompressibleNavierStokes.h"

#include "M02IncompressibleNavierStokes_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<M02IncompressibleNavierStokes> M02IncompressibleNavierStokes::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<M02IncompressibleNavierStokes>(new M02IncompressibleNavierStokes(para, level));
}

void M02IncompressibleNavierStokes::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    M02IncompressibleNavierStokes_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_MRT_Incomp_SP_27 execution failed");
}

M02IncompressibleNavierStokes::M02IncompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	
}

M02IncompressibleNavierStokes::M02IncompressibleNavierStokes()
{
}
