#include "K15IncompressibleNavierStokesRotatingVelocityField.h"

#include "K15IncompressibleNavierStokesRotatingVelocityField_Device.cuh"
#include "Parameter/Parameter.h"
#include "cuda/CudaGrid.h"

std::shared_ptr<K15IncompressibleNavierStokesRotatingVelocityField> K15IncompressibleNavierStokesRotatingVelocityField::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K15IncompressibleNavierStokesRotatingVelocityField>(new K15IncompressibleNavierStokesRotatingVelocityField(para, level));
}

void K15IncompressibleNavierStokesRotatingVelocityField::run()
{
    vf::cuda::CudaGrid grid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);

    K15IncompressibleNavierStokesRotatingVelocityField_Device <<< grid.grid, grid.threads >>>(
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
    getLastCudaError("K15IncompressibleNavierStokesRotatingVelocityField_Device execution failed");
}

K15IncompressibleNavierStokesRotatingVelocityField::K15IncompressibleNavierStokesRotatingVelocityField(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	
}

K15IncompressibleNavierStokesRotatingVelocityField::K15IncompressibleNavierStokesRotatingVelocityField()
{
}
