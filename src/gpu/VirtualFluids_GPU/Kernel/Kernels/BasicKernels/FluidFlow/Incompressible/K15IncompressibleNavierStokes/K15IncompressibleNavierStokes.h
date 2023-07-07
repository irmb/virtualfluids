#ifndef K15IncompressibleNavierStokes_H
#define K15IncompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class K15IncompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr<K15IncompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K15IncompressibleNavierStokes();
	K15IncompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 