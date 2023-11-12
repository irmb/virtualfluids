#include "InitSP27.h"

#include "InitSP27_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<PreProcessorStrategy> InitSP27::getNewInstance(std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<PreProcessorStrategy>(new InitSP27(para));
}

void InitSP27::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Init_SP_27 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->rho,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Init_SP_27 execution failed");
}

bool InitSP27::checkParameter()
{
    return false;
}

InitSP27::InitSP27(std::shared_ptr<Parameter> para)
{
    this->para = para;
}

InitSP27::InitSP27()
{
}
