#ifndef KERNEL_H
#define KERNEL_H

#include <DataTypes.h>

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "Kernel/Utilities/KernelGroup.h"
#include "PreProcessor/PreProcessorType.h"

#include <vector>


class Kernel
{
public:
	virtual void run() = 0;

	virtual bool checkParameter() = 0;
	virtual std::vector<PreProcessorType> getPreProcessorTypes() = 0;
	virtual KernelGroup getKernelGroup() = 0;
};
#endif