#include "InitIncompAD27.h"

#include "InitIncompAD27_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<PreProcessorStrategy> InitIncompAD27::getNewInstance(std::shared_ptr<Parameter> para)
{
	return std::shared_ptr<PreProcessorStrategy>(new InitIncompAD27(para));
}

void InitIncompAD27::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Init_Incomp_AD_27 <<< grid.grid, grid.threads >>>(
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
    getLastCudaError("LB_Init_Incomp_AD_27 execution failed");
}

bool InitIncompAD27::checkParameter()
{
	return false;
}

InitIncompAD27::InitIncompAD27(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

InitIncompAD27::InitIncompAD27()
{
}
