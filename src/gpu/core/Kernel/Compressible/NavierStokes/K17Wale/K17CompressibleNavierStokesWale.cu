#include "K17CompressibleNavierStokesWale.h"

#include "K17CompressibleNavierStokesWale_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<K17CompressibleNavierStokesWale> K17CompressibleNavierStokesWale::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<K17CompressibleNavierStokesWale>(new K17CompressibleNavierStokesWale(para, level));
}

void K17CompressibleNavierStokesWale::run()
{
	int size_Mat = (int)para->getParD(level)->numberOfNodes;
	int numberOfThreads = para->getParD(level)->numberofthreads;

	//int Grid = size_Array / numberOfThreads;
	//dim3 grid(Grid, 1, 1);
	//dim3 threads(numberOfThreads, 1, 1 );

	int Grid = (size_Mat / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid > 512)
	{
		Grid1 = 512;
		Grid2 = (Grid / Grid1) + 1;
	}
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2, 1);
	dim3 threads(numberOfThreads, 1, 1);

	K17CompressibleNavierStokesWale_Device <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->typeOfGridNode,
		para->getParD(level)->neighborX,
		para->getParD(level)->neighborY,
		para->getParD(level)->neighborZ,
		para->getParD(level)->neighborInverse,
		para->getParD(level)->velocityX,
		para->getParD(level)->velocityY,
		para->getParD(level)->velocityZ,
		para->getParD(level)->distributions.f[0],
		para->getParD(level)->turbViscosity,
		para->getParD(level)->numberOfNodes,
		level,
		para->getTimestepOfCoarseLevel(),
		para->getForcesDev(),
        para->getQuadricLimitersDev(),
		para->getParD(level)->isEvenTimestep);
	getLastCudaError("K17CompressibleNavierStokesWale_Device execution failed");
}

K17CompressibleNavierStokesWale::K17CompressibleNavierStokesWale(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	
}

K17CompressibleNavierStokesWale::K17CompressibleNavierStokesWale()
{
}
