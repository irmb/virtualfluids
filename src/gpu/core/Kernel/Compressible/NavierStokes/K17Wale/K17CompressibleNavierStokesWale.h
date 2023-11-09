#ifndef K17CompressibleNavierStokesWale_H
#define K17CompressibleNavierStokesWale_H

#include "Kernel/KernelImp.h"

class K17CompressibleNavierStokesWale : public KernelImp
{
public:
	static std::shared_ptr<K17CompressibleNavierStokesWale> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K17CompressibleNavierStokesWale();
	K17CompressibleNavierStokesWale(std::shared_ptr< Parameter> para, int level);

};
#endif 