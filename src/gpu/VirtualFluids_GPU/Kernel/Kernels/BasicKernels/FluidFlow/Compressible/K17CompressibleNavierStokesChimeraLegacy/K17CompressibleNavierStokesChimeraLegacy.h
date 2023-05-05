#ifndef K17CompressibleNavierStokesChimeraLegacy_H
#define K17CompressibleNavierStokesChimeraLegacy_H

#include "Kernel/KernelImp.h"

class K17CompressibleNavierStokesChimeraLegacy : public KernelImp
{
public:
	static std::shared_ptr<K17CompressibleNavierStokesChimeraLegacy> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
    K17CompressibleNavierStokesChimeraLegacy();
    K17CompressibleNavierStokesChimeraLegacy(std::shared_ptr<Parameter> para, int level);
};

#endif 
