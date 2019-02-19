#ifndef KERNEL_H
#define KERNEL_H

#include <DataTypes.h>

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

class Parameter;

class Kernel
{
public:
	virtual void run()=0;
	virtual bool checkParameter() = 0;
	
};
#endif