#include "Cumulant1hIncompSP27.h"

#include "Cumulant1hIncompSP27_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<Cumulant1hIncompSP27> Cumulant1hIncompSP27::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<Cumulant1hIncompSP27>(new Cumulant1hIncompSP27(para, level));
}

void Cumulant1hIncompSP27::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    LB_Kernel_Cum_1h_Incomp_SP_27 <<< grid.grid, grid.threads >>>(
        para->getParD(level)->omega,
        para->getParD(level)->deltaPhi,
        para->getAngularVelocity(),
        para->getParD(level)->typeOfGridNode,
        para->getParD(level)->neighborX,
        para->getParD(level)->neighborY,
        para->getParD(level)->neighborZ,
        para->getParD(level)->coordinateX,
        para->getParD(level)->coordinateY,
        para->getParD(level)->coordinateZ,
        para->getParD(level)->distributions.f[0],
        para->getParD(level)->numberOfNodes,
        para->getParD(level)->isEvenTimestep);
    getLastCudaError("LB_Kernel_Cum_1h_Incomp_SP_27 execution failed");
}

Cumulant1hIncompSP27::Cumulant1hIncompSP27(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	
}

Cumulant1hIncompSP27::Cumulant1hIncompSP27()
{
}
