#ifndef C06IncompressibleNavierStokes_H
#define C06IncompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class C06IncompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr<C06IncompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	C06IncompressibleNavierStokes();
	C06IncompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);

};
#endif 