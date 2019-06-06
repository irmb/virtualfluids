#include "PMCumulantOneCompSP27.h"

#include "PMCumulantOneCompSP27_Device.cuh"
#include "Parameter\Parameter.h"
#include "Calculation/PorousMedia.h"

std::shared_ptr<PMCumulantOneCompSP27> PMCumulantOneCompSP27::getNewInstance(std::shared_ptr<Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level)
{
	return std::shared_ptr<PMCumulantOneCompSP27>(new PMCumulantOneCompSP27(para, pm, level));
}

void PMCumulantOneCompSP27::run()
{
	int size_Mat = para->getParD(level)->size_Mat_SP;
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
		LB_Kernel_PM_Cum_One_Comp_SP_27 << < grid, threads >> >(para->getParD(level)->omega,
			para->getParD(level)->neighborX_SP,
			para->getParD(level)->neighborY_SP,
			para->getParD(level)->neighborZ_SP,
			para->getParD(level)->d0SP.f[0],
			para->getParD(level)->size_Mat_SP,
			level,
			para->getForcesDev(),
			pm[i]->getPorosity(),
			pm[i]->getDarcyLBM(),
			pm[i]->getForchheimerLBM(),
			pm[i]->getSizePM(),
			pm[i]->getHostNodeIDsPM(),
			para->getParD(level)->evenOrOdd);
		getLastCudaError("LB_Kernel_PM_Cum_One_Comp_SP_27 execution failed");
	}
}

PMCumulantOneCompSP27::PMCumulantOneCompSP27(std::shared_ptr<Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level)
{
	this->para = para;
	this->pm = pm;
	this->level = level;

	myPreProcessorTypes.push_back(InitSP27);

	myKernelGroup = BasicKernel;
}

PMCumulantOneCompSP27::PMCumulantOneCompSP27()
{
}
