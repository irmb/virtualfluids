#include "K17CompressibleNavierStokesChimeraLegacy.h"

#include "Parameter/Parameter.h"
#include "K17CompressibleNavierStokesChimeraLegacy_Device.cuh"
#include "cuda/CudaGrid.h"

std::shared_ptr<K17CompressibleNavierStokesChimeraLegacy> K17CompressibleNavierStokesChimeraLegacy::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K17CompressibleNavierStokesChimeraLegacy>(new K17CompressibleNavierStokesChimeraLegacy(para,level));
}

void K17CompressibleNavierStokesChimeraLegacy::run()
{
	K17CompressibleNavierStokesChimeraLegacy_Device <<< cudaGrid.grid, cudaGrid.threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->typeOfGridNode,
		para->getParD(level)->neighborX,
		para->getParD(level)->neighborY,
		para->getParD(level)->neighborZ,
		para->getParD(level)->distributions.f[0],
		para->getParD(level)->numberOfNodes,
		level,
		para->getIsBodyForce(),
		para->getForcesDev(),
		para->getParD(level)->forceX_SP,
		para->getParD(level)->forceY_SP,
		para->getParD(level)->forceZ_SP,
        para->getQuadricLimitersDev(),
		para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

K17CompressibleNavierStokesChimeraLegacy::K17CompressibleNavierStokesChimeraLegacy(std::shared_ptr<Parameter> para, int level): KernelImp(para, level)
{
	myPreProcessorTypes.push_back(InitCompSP27);
	
	this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
}