#include "B92CompressibleNavierStokes.h"

#include "B92CompressibleNavierStokes_Device.cuh"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>

std::shared_ptr<B92CompressibleNavierStokes> B92CompressibleNavierStokes::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<B92CompressibleNavierStokes>(new B92CompressibleNavierStokes(para, level));
}

void B92CompressibleNavierStokes::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    B92CompressibleNavierStokes_Device<<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_BGK_Comp_SP_27 execution failed");
}

B92CompressibleNavierStokes::B92CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitCompSP27);
    
}

B92CompressibleNavierStokes::B92CompressibleNavierStokes()
{
}
