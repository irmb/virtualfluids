#include "KernelFactory.h"

#include "Kernel/Advection/Compressible/CumulantOne/CumulantOneCompSP27.h"
#include "Kernel/Advection/Compressible/CumulantAA2016/CumulantAA2016CompSP27.h"
#include "Kernel/Advection/Compressible/CumulantAll4/CumulantAll4CompSP27.h"
#include "Kernel/Advection/Compressible/CumulantF3/CumulantF3CompSP27.h"
#include "Kernel/Advection/Compressible/CumulantF32018/CumulantF32018CompSP27.h"

std::shared_ptr<KernelFactory> KernelFactory::getNewInstance(std::shared_ptr<Parameter> para)
{
	return std::shared_ptr<KernelFactory>(new KernelFactory(para));
}

std::vector<std::shared_ptr<Kernel>> KernelFactory::makeKernels(int maxLevel, std::string kernelName)
{
	kernels.resize(0);
	for (int i = 0; i <= maxLevel; i++)
	{
		kernels.push_back(makeKernel(kernelName, i));
	}
	return kernels;
}

void KernelFactory::setKernelAtLevel(std::vector<std::shared_ptr<Kernel>> kernels, int i, std::string kernelName)
{
	kernels.at(i) = makeKernel(kernelName, i);
}

std::shared_ptr<Kernel> KernelFactory::makeKernel(std::string kernelName, int level)
{
	if (kernelName == "CumulantOneCompSP27")
		return CumulantOneCompSP27::getNewInstance(para, level);
	if (kernelName == "CumulantAA2016CompSP27")
		return CumulantAA2016CompSP27::getNewInstance(para, level);
	if (kernelName == "CumulantAll4CompSP27")
		return CumulantAll4CompSP27::getNewInstance(para, level);
	if (kernelName == "CumulantF3CompSP27")
		return CumulantF3CompSP27::getNewInstance(para, level);
	if (kernelName == "CumulantF32018CompSP27")
		return CumulantF32018CompSP27::getNewInstance(para, level);
}

KernelFactory::KernelFactory(std::shared_ptr<Parameter> para)
{
	this->para = para;
	kernels.resize(0);
}