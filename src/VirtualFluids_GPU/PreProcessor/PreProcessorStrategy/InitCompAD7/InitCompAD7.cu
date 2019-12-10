#include "InitCompAD7.h"

#include "InitCompAD7_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<InitCompAD7> InitCompAD7::getNewInstance(std::shared_ptr<Parameter> para)
{
	return std::shared_ptr<InitCompAD7>(new InitCompAD7(para));
}

void InitCompAD7::init(int level)
{
	int numberOfThreads = para->getParD(level)->numberofthreads;
	int size_Mat = para->getParD(level)->size_Mat_SP;

	int Grid = (size_Mat / numberOfThreads) + 1;
	int Grid1, Grid2;
	if (Grid>512)
	{
		Grid1 = 512;
		Grid2 = (Grid / Grid1) + 1;
	}
	else
	{
		Grid1 = 1;
		Grid2 = Grid;
	}
	dim3 grid(Grid1, Grid2);
	dim3 threads(numberOfThreads, 1, 1);

	LB_Init_Comp_AD_7 << < grid, threads >> >(	para->getParD(level)->neighborX_SP,
										para->getParD(level)->neighborY_SP,
										para->getParD(level)->neighborZ_SP,
										para->getParD(level)->geoSP,
										para->getParD(level)->Conc,
										para->getParD(level)->vx_SP,
										para->getParD(level)->vy_SP,
										para->getParD(level)->vz_SP,
										para->getParD(level)->size_Mat_SP,
										para->getParD(level)->d7.f[0],
										para->getParD(level)->evenOrOdd);
	getLastCudaError("LBInitThS7 execution failed");
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
