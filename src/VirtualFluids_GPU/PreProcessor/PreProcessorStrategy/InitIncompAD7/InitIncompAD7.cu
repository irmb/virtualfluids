#include "InitIncompAD7.h"

#include "InitIncompAD7_Device.cuh"
#include "Parameter/Parameter.h"

std::shared_ptr<PreProcessorStrategy> InitIncompAD7::getNewInstance(std::shared_ptr<Parameter> para)
{
	return std::shared_ptr<PreProcessorStrategy>(new InitIncompAD7(para));
}

void InitIncompAD7::init(int level)
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

	LB_Init_Incomp_AD_7 << < grid, threads >> >(	para->getParD(level)->neighborX_SP,
												para->getParD(level)->neighborY_SP,
												para->getParD(level)->neighborZ_SP,
												para->getParD(level)->geoSP,
												para->getParD(level)->Conc,
												para->getParD(level)->vx_SP,
												para->getParD(level)->vy_SP,
												para->getParD(level)->vz_SP,
												para->getParD(level)->size_Mat_SP,
												para->getParD(level)->d27.f[0],
												para->getParD(level)->evenOrOdd);
	getLastCudaError("LBInitIncompAD7 execution failed");
}

bool InitIncompAD7::checkParameter()
{
	return false;
}

InitIncompAD7::InitIncompAD7(std::shared_ptr<Parameter> para)
{
	this->para = para;
}

InitIncompAD7::InitIncompAD7()
{
}
