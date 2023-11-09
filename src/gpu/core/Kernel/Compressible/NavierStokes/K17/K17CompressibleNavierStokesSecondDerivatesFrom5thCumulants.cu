#include "K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants.h"

#include "K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda_helper/CudaGrid.h"

std::shared_ptr<K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants> K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants>(new K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants(para, level));
}

void K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants_Device <<< grid.grid, grid.threads >>>(
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

K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants::K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}