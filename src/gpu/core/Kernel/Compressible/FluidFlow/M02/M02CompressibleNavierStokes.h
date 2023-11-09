#ifndef M02CompressibleNavierStokes_H
#define M02CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"


class M02CompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr<M02CompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	M02CompressibleNavierStokes();
	M02CompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};
#endif 