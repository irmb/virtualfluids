#include "TurbulentViscosityCumulantK17CompChim.h"
#include "cuda/CudaGrid.h"
#include "Parameter/Parameter.h"
#include "TurbulentViscosityCumulantK17CompChim_Device.cuh"

std::shared_ptr<TurbulentViscosityCumulantK17CompChim> TurbulentViscosityCumulantK17CompChim::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<TurbulentViscosityCumulantK17CompChim>(new TurbulentViscosityCumulantK17CompChim(para,level));
}

void TurbulentViscosityCumulantK17CompChim::run()
{
	vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParH(level)->numberofthreads, para->getParH(level)->numberOfNodes);

	LB_Kernel_TurbulentViscosityCumulantK17CompChim <<< grid.grid, grid.threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->typeOfGridNode,
		para->getParD(level)->neighborX,
		para->getParD(level)->neighborY,
		para->getParD(level)->neighborZ,
		para->getParD(level)->distributions.f[0],
		para->getParD(level)->rho,
		para->getParD(level)->velocityX,
		para->getParD(level)->velocityY,
		para->getParD(level)->velocityZ,
		para->getParD(level)->turbViscosity,
		(unsigned long)para->getParD(level)->numberOfNodes,
		level,
		para->getIsBodyForce(),
		para->getForcesDev(),
		para->getParD(level)->forceX_SP,
		para->getParD(level)->forceY_SP,
		para->getParD(level)->forceZ_SP,
        para->getQuadricLimitersDev(),
		para->getParD(level)->isEvenTimestep);
	getLastCudaError("LB_Kernel_TurbulentViscosityCumulantK17CompChim execution failed");
}

TurbulentViscosityCumulantK17CompChim::TurbulentViscosityCumulantK17CompChim(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}