#include "B15IncompressibleNavierStokesBGKplus.h"

#include "B15IncompressibleNavierStokesBGKplus_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<B15IncompressibleNavierStokesBGKplus> B15IncompressibleNavierStokesBGKplus::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<B15IncompressibleNavierStokesBGKplus>(new B15IncompressibleNavierStokesBGKplus(para, level));
}

void B15IncompressibleNavierStokesBGKplus::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    B15IncompressibleNavierStokesBGKplus_Device <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_BGK_Plus_Incomp_SP_27 execution failed");
}

B15IncompressibleNavierStokesBGKplus::B15IncompressibleNavierStokesBGKplus(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitSP27);

    
}

B15IncompressibleNavierStokesBGKplus::B15IncompressibleNavierStokesBGKplus()
{
}
