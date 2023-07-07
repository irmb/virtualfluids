#ifndef K15IncompressibleNavierStokesIsoCheck_H
#define K15IncompressibleNavierStokesIsoCheck_H

#include "Kernel/KernelImp.h"

class K15IncompressibleNavierStokesIsoCheck : public KernelImp
{
public:
	static std::shared_ptr<K15IncompressibleNavierStokesIsoCheck> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K15IncompressibleNavierStokesIsoCheck();
	K15IncompressibleNavierStokesIsoCheck(std::shared_ptr< Parameter> para, int level);

};

#endif 