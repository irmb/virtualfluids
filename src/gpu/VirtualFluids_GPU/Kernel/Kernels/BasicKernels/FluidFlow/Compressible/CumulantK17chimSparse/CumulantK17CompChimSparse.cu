#include "CumulantK17CompChimSparse.h"

#include "Parameter/Parameter.h"
#include "CumulantK17CompChimSparse_Device.cuh"

std::shared_ptr<CumulantK17CompChimSparse> CumulantK17CompChimSparse::getNewInstance(std::shared_ptr<Parameter> para,
                                                                               int level)
{
    return std::shared_ptr<CumulantK17CompChimSparse>(new CumulantK17CompChimSparse(para, level));
}

void CumulantK17CompChimSparse::run()
{
    dim3 grid, threads;
    std::tie(grid, threads) = *calcGridDimensions(para->getParD(level)->numberOfFluidNodes);

	LB_Kernel_CumulantK17CompChimSparse <<< grid, threads >>>(
		para->getParD(level)->omega,
		para->getParD(level)->neighborX_SP,
		para->getParD(level)->neighborY_SP,
		para->getParD(level)->neighborZ_SP,
		para->getParD(level)->d0SP.f[0],
		para->getParD(level)->size_Mat_SP,
		level,
		para->getForcesDev(),
        para->getQuadricLimitersDev(),
		para->getParD(level)->evenOrOdd,
        para->getParD(level)->fluidNodeIndices,
		para->getParD(level)->numberOfFluidNodes);
	getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

void CumulantK17CompChimSparse::runOnIndices(const unsigned int *indices, unsigned int size_indices)
{
    dim3 grid, threads;
    std::tie(grid, threads) = *calcGridDimensions(para->getParD(level)->numberOfFluidNodes);

    LB_Kernel_CumulantK17CompChimSparse<<<grid, threads, 0, para->getStream(0)>>>(
        para->getParD(level)->omega, 
		para->getParD(level)->neighborX_SP, 
		para->getParD(level)->neighborY_SP,
        para->getParD(level)->neighborZ_SP, 
		para->getParD(level)->d0SP.f[0], 
		para->getParD(level)->size_Mat_SP, 
		level,
        para->getForcesDev(), 
		para->getQuadricLimitersDev(), 
		para->getParD(level)->evenOrOdd,
        indices,
		size_indices);
    getLastCudaError("LB_Kernel_CumulantK17CompChim execution failed");
}

CumulantK17CompChimSparse::CumulantK17CompChimSparse(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;
}

std::unique_ptr<std::pair<dim3, dim3>> CumulantK17CompChimSparse::calcGridDimensions(unsigned int size_Mat)
{
    int numberOfThreads = para->getParD(level)->numberofthreads;

    int Grid = (size_Mat / numberOfThreads) + 1;
    int Grid1, Grid2;
    if (Grid > 512) {
        Grid1 = 512;
        Grid2 = (Grid / Grid1) + 1;
    } else {
        Grid1 = 1;
        Grid2 = Grid;
    }
    dim3 grid(Grid1, Grid2);
    dim3 threads(numberOfThreads, 1, 1);
    std::pair<dim3, dim3> dimensions(grid, threads);
    return std::make_unique<std::pair<dim3, dim3>>(dimensions);
}
