#ifndef M02IncompressibleNavierStokes_H
#define M02IncompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class M02IncompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr<M02IncompressibleNavierStokes> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run();

private:
	M02IncompressibleNavierStokes();
	M02IncompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 