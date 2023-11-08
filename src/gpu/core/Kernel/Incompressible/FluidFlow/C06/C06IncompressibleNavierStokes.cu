#include "C06IncompressibleNavierStokes.h"

#include "CascadeIncompSP27_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<C06IncompressibleNavierStokes> C06IncompressibleNavierStokes::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<C06IncompressibleNavierStokes>(new C06IncompressibleNavierStokes(para, level));
}

void C06IncompressibleNavierStokes::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    C06IncompressibleNavierStokes_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_Cascade_Incomp_SP_27 execution failed");
}

C06IncompressibleNavierStokes::C06IncompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	
}

C06IncompressibleNavierStokes::C06IncompressibleNavierStokes()
{
}
