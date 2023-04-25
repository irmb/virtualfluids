#include "CumulantK17CompChim.h"

#include "Parameter/Parameter.h"
#include "CumulantK17CompChim_Device.cuh"
#include "cuda/CudaGrid.h"

std::shared_ptr<CumulantK17CompChim> CumulantK17CompChim::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK17CompChim>(new CumulantK17CompChim(para,level));
}

void CumulantK17CompChim::run()
{
	LB_Kernel_CumulantK17CompChim <<< cudaGrid.grid, cudaGrid.threads >>>(
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

CumulantK17CompChim::CumulantK17CompChim(std::shared_ptr<Parameter> para, int level): KernelImp(para, level)
{
	myPreProcessorTypes.push_back(InitCompSP27);
	
	this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
}