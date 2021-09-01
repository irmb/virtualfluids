#include "CumulantK17CompChim.h"

#include "Parameter/Parameter.h"
#include "CumulantK17CompChim_Device.cuh"

std::shared_ptr<CumulantK17CompChim> CumulantK17CompChim::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK17CompChim>(new CumulantK17CompChim(para,level));
}

void CumulantK17CompChim::run()
{
    dim3 grid, threads;
    std::tie(grid, threads) =
        *calcGridDimensions(para->getParD(level)->size_Mat_SP, para->getParD(level)->numberofthreads);

	LB_Kernel_CumulantK17CompChim <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->geoSP,
		para->getParD(level)->neighborX_SP,
		para->getParD(level)->neighborY_SP,
		para->getParD(level)->neighborZ_SP,
		para->getParD(level)->d0SP.f[0],
		para->getParD(level)->size_Mat_SP,
		level,
		para->getForcesDev(),
        para->getQuadricLimitersDev(),
		para->getParD(level)->evenOrOdd);
	getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

CumulantK17CompChim::CumulantK17CompChim(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}