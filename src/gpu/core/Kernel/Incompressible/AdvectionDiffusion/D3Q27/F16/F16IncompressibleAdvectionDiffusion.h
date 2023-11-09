#ifndef F16IncompressibleAdvectionDiffusion_H
#define F16IncompressibleAdvectionDiffusion_H

#include "Kernel/AdvectionDiffusionKernel.h"


class F16IncompressibleAdvectionDiffusion : public AdvectionDiffusionKernel
{
public:
	static std::shared_ptr<F16IncompressibleAdvectionDiffusion> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run();

private:
	F16IncompressibleAdvectionDiffusion();
	F16IncompressibleAdvectionDiffusion(std::shared_ptr< Parameter> para, int level);
};
#endif 