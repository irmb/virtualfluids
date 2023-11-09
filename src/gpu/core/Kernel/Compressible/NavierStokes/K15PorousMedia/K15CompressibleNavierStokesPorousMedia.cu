#include "K15CompressibleNavierStokesPorousMedia.h"

#include "K15CompressibleNavierStokesPorousMedia_Device.cuh"
#include "Parameter/Parameter.h"
#include "Calculation/PorousMedia.h"

std::shared_ptr<K15CompressibleNavierStokesPorousMedia> K15CompressibleNavierStokesPorousMedia::getNewInstance(std::shared_ptr<Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level)
{
	return std::shared_ptr<K15CompressibleNavierStokesPorousMedia>(new K15CompressibleNavierStokesPorousMedia(para, pm, level));
}

void K15CompressibleNavierStokesPorousMedia::run()
{
	int size_Mat = (int)para->getParD(level)->numberOfNodes;
	int numberOfThreads = para->getParD(level)->numberofthreads;

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

	for (int i = 0; i < pm.size(); i++) {
		K15CompressibleNavierStokesPorousMedia_Device <<< grid, threads >>>(
			para->getParD(level)->omega,
			para->getParD(level)->neighborX,
			para->getParD(level)->neighborY,
			para->getParD(level)->neighborZ,
			para->getParD(level)->distributions.f[0],
			para->getParD(level)->numberOfNodes,
			level,
			para->getForcesDev(),
			pm[i]->getPorosity(),
			pm[i]->getDarcyLBM(),
			pm[i]->getForchheimerLBM(),
			pm[i]->getSizePM(),
			pm[i]->getHostNodeIDsPM(),
			para->getParD(level)->isEvenTimestep);
		getLastCudaError("K15CompressibleNavierStokesPorousMedia_Device execution failed");
	}
}

K15CompressibleNavierStokesPorousMedia::K15CompressibleNavierStokesPorousMedia(std::shared_ptr<Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level)
{
	this->para = para;
	this->pm = pm;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	
}

K15CompressibleNavierStokesPorousMedia::K15CompressibleNavierStokesPorousMedia()
{
}
