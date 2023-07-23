#include "InitIncompAD7.h"

#include "InitIncompAD7_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<PreProcessorStrategy> InitIncompAD7::getNewInstance(std::shared_ptr<Parameter> para)
{
	return std::shared_ptr<PreProcessorStrategy>(new InitIncompAD7(para));
}

void InitIncompAD7::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Init_Incomp_AD_7 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->concentration,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->distributionsAD.f[0],
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Init_Incomp_AD_7 execution failed");
}

bool InitIncompAD7::checkParameter()
{
	return false;
}

InitIncompAD7::InitIncompAD7(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

InitIncompAD7::InitIncompAD7()
{
}
