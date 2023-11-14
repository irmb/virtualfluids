#include "InitCompAD27.h"

#include "InitCompAD27_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<PreProcessorStrategy> InitCompAD27::getNewInstance(std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<PreProcessorStrategy>(new InitCompAD27(para));
}

void InitCompAD27::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Init_Comp_AD_27 <<< grid.grid, grid.threads >>>(
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
    getLastCudaError("LB_Init_Comp_AD_27 execution failed");
}

bool InitCompAD27::checkParameter()
{
    return false;
}

InitCompAD27::InitCompAD27(std::shared_ptr<Parameter> para)
{
    this->para = para;
}

InitCompAD27::InitCompAD27()
{
}
