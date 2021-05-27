#ifndef KERNEL_IMP_H
#define KERNEL_IMP_H

#include "Kernel.h"

#include <memory>

#include "Utilities/CudaGrid.h"

class CheckParameterStrategy;
class Parameter;

class KernelImp : public Kernel
{
public:
	virtual void run() = 0;

	bool checkParameter();
	std::vector<PreProcessorType> getPreProcessorTypes();
	KernelGroup getKernelGroup();

	void setCheckParameterStrategy(std::shared_ptr<CheckParameterStrategy> strategy);

protected:
	KernelImp();

	std::shared_ptr< Parameter> para;
	std::shared_ptr<CheckParameterStrategy> checkStrategy;
	int level;
	std::vector<PreProcessorType> myPreProcessorTypes;
	KernelGroup myKernelGroup;

	vf::gpu::CudaGrid cudaGrid;

};
#endif