#include "InitF3.h"

#include "InitF3_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<PreProcessorStrategy> InitF3::getNewInstance(std::shared_ptr<Parameter> para)
{
	return std::shared_ptr<PreProcessorStrategy>(new InitF3(para));
}

void InitF3::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Init_F3 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->rho,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->g6.g[0],
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Init_F3 execution failed");
}

bool InitF3::checkParameter()
{
	return false;
}

InitF3::InitF3(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

InitF3::InitF3()
{
}
