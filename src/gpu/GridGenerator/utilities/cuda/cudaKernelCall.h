#ifndef cudaKernelCall_H
#define cudaKernelCall_H

#include "utilities/cuda/cudaDefines.h"
#include "utilities/cuda/CudaErrorCheck.cu"
#include "utilities/cuda/LaunchParameter.cuh"

template<typename Functor, typename... TArgs>
CUDA_HOST float runKernel(Functor kernel, const LaunchParameter& para, TArgs... args)
{
	para.print();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start, 0);
	kernel << < para.blocks, para.threads >> >(args...);
	CudaCheckError();
	cudaDeviceSynchronize();
	cudaEventRecord(stop, 0);

	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return elapsedTime;
}


#endif 
