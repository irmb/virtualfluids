#ifndef K20CompressibleNavierStokes_H
#define K20CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class K20CompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr< K20CompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K20CompressibleNavierStokes();
	K20CompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 