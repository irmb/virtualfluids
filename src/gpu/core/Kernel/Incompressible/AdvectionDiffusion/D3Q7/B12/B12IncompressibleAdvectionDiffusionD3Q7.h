#ifndef B12IncompressibleAdvectionDiffusionD3Q7_H
#define B12IncompressibleAdvectionDiffusionD3Q7_H

#include "Kernel/ADKernel.h"

class B12IncompressibleAdvectionDiffusionD3Q7 : public ADKernel
{
public:
	static std::shared_ptr<B12IncompressibleAdvectionDiffusionD3Q7> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	B12IncompressibleAdvectionDiffusionD3Q7();
	B12IncompressibleAdvectionDiffusionD3Q7(std::shared_ptr< Parameter> para, int level);
};
#endif 