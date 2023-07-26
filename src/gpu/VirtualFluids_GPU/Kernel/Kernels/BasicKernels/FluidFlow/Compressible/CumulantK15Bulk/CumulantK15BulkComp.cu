#include "CumulantK15BulkComp.h"

#include "CumulantK15BulkComp_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<CumulantK15BulkComp> CumulantK15BulkComp::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK15BulkComp>(new CumulantK15BulkComp(para, level));
}

void CumulantK15BulkComp::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_CumulantK15BulkComp <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        level,
        para->getForcesDev(),
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_CumulantK15BulkComp execution failed");
}

CumulantK15BulkComp::CumulantK15BulkComp(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

CumulantK15BulkComp::CumulantK15BulkComp()
{
}
