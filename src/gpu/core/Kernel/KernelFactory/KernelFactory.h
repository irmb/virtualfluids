#ifndef KERNEL_FACTORY_H
#define KERNEL_FACTORY_H


#include <memory>
#include <string>
#include <vector>



class Kernel;
class AdvectionDiffusionKernel;
class Parameter;
class PorousMedia;

class KernelFactory
{
public:
    virtual ~KernelFactory() = default;
    virtual std::vector< std::shared_ptr< Kernel>> makeKernels(std::shared_ptr<Parameter> para) = 0;
    virtual std::vector< std::shared_ptr< AdvectionDiffusionKernel>> makeAdvDifKernels(std::shared_ptr<Parameter> para) = 0;
    
    virtual void setPorousMedia(std::vector<std::shared_ptr<PorousMedia>> pm) = 0;
};
#endif