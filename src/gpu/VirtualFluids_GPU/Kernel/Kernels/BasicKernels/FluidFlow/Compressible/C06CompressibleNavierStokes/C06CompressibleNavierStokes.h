#ifndef C06CompressibleNavierStokes_H
#define C06CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class C06CompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr<C06CompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	C06CompressibleNavierStokes();
	C06CompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 