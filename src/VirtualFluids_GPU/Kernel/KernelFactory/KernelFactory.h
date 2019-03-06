#ifndef KERNEL_FACTORY_H
#define KERNEL_FACTORY_H

#include <memory>
#include <vector>
#include <iostream>

class Kernel;
class Parameter;

enum kernelName { CumulantOneCompSP27, CumulantAA2016CompSP27, CumulantAll4CompSP27, CumulantF3CompSP27, CumulantF32018CompSP27 };

class KernelFactory
{
public:
	static std::shared_ptr< KernelFactory> getNewInstance(std::shared_ptr<Parameter> para);
	std::vector< std::shared_ptr< Kernel>> makeKernels(int maxLevel, std::string kernelName);
	void setKernelAtLevel(std::vector< std::shared_ptr<Kernel>> kernels, int i, std::string kernelName);

private:
	std::shared_ptr< Kernel> makeKernel(std::string kernelName, int level);
	KernelFactory(std::shared_ptr<Parameter> para);
	KernelFactory();

	std::vector< std::shared_ptr< Kernel>> kernels;
	std::shared_ptr< Parameter> para;
};
#endif