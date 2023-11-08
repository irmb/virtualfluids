#include "InitCompSP27.h"

#include "InitCompSP27_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<PreProcessorStrategy> InitCompSP27::getNewInstance(std::shared_ptr<Parameter> para)
{
	return std::shared_ptr<PreProcessorStrategy>(new InitCompSP27(para));
}

void InitCompSP27::init(int level)
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    if( ! para->getUseInitNeq() )
    {
        LB_Init_Comp_SP_27 <<< grid.grid, grid.threads >>> (
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
        getLastCudaError("LB_Init_Comp_SP_27 execution failed");
    }
    else
    {
        LB_Init_Comp_Neq_SP_27 <<< grid.grid, grid.threads >>> (
            para->getParD(level)->neighborX,
            para->getParD(level)->neighborY,
            para->getParD(level)->neighborZ,
            para->getParD(level)->neighborInverse,
            para->getParD(level)->typeOfGridNode,
            para->getParD(level)->rho,
            para->getParD(level)->velocityX,
            para->getParD(level)->velocityY,
            para->getParD(level)->velocityZ,
            para->getParD(level)->numberOfNodes,
            para->getParD(level)->distributions.f[0],
            para->getParD(level)->omega,
            para->getParD(level)->isEvenTimestep);
        cudaDeviceSynchronize();
        getLastCudaError("LB_Init_Comp_Neq_SP_27 execution failed");
    }



}

bool InitCompSP27::checkParameter()
{
	return false;
}

InitCompSP27::InitCompSP27(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

InitCompSP27::InitCompSP27()
{
}
