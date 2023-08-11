#ifndef K08CompressibleNavierStokes_H
#define K08CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class K08CompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr<K08CompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K08CompressibleNavierStokes();
	K08CompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);

};
#endif 