#ifndef K15CompressibleNavierStokesWale_H
#define K15CompressibleNavierStokesWale_H

#include "Kernel/KernelImp.h"

class K15CompressibleNavierStokesWale : public KernelImp
{
public:
	static std::shared_ptr<K15CompressibleNavierStokesWale> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K15CompressibleNavierStokesWale();
	K15CompressibleNavierStokesWale(std::shared_ptr< Parameter> para, int level);

};
#endif 