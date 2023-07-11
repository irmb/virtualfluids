#include "K15IncompressibleNavierStokesIsoCheck.h"

#include "K15IncompressibleNavierStokesIsoCheck_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<K15IncompressibleNavierStokesIsoCheck> K15IncompressibleNavierStokesIsoCheck::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K15IncompressibleNavierStokesIsoCheck>(new K15IncompressibleNavierStokesIsoCheck(para, level));
}

void K15IncompressibleNavierStokesIsoCheck::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K15IncompressibleNavierStokesIsoCheck_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->dxxUx,
        para->getParD(level)->dyyUy,
        para->getParD(level)->dzzUz,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_Cum_IsoTest_Incomp_SP_27 execution failed");
}

K15IncompressibleNavierStokesIsoCheck::K15IncompressibleNavierStokesIsoCheck(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	
}

K15IncompressibleNavierStokesIsoCheck::K15IncompressibleNavierStokesIsoCheck()
{
}
