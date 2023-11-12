#include "InitCompAD7.h"

#include "InitCompAD7_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<InitCompAD7> InitCompAD7::getNewInstance(std::shared_ptr<Parameter> para)
{
    return std::shared_ptr<InitCompAD7>(new InitCompAD7(para));
}

void InitCompAD7::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Init_Comp_AD_7 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->concentration,
        para->getParD(level)->velocityX,
        para->getParD(level)->velocityY,
        para->getParD(level)->velocityZ,
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->distributionsAD7.f[0],
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Init_Comp_AD_7 execution failed");
}

bool InitCompAD7::checkParameter()
{
    return false;
}

InitCompAD7::InitCompAD7(std::shared_ptr<Parameter> para)
{
    this->para = para;
}

InitCompAD7::InitCompAD7()
{
}
