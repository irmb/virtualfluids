#ifndef F16CompressibleAdvectionDiffusion_H
#define F16CompressibleAdvectionDiffusion_H

#include "Kernel/ADKernel.h"

class F16CompressibleAdvectionDiffusion : public ADKernel
{
public:
	static std::shared_ptr<F16CompressibleAdvectionDiffusion> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	F16CompressibleAdvectionDiffusion();
	F16CompressibleAdvectionDiffusion(std::shared_ptr< Parameter> para, int level);
};
#endif 