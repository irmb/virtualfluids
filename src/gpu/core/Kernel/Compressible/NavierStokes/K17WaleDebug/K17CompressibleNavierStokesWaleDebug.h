#ifndef K17CompressibleNavierStokesWaleDebug_H
#define K17CompressibleNavierStokesWaleDebug_H

#include "Kernel/KernelImp.h"

class K17CompressibleNavierStokesWaleDebug : public KernelImp
{
public:
	static std::shared_ptr<K17CompressibleNavierStokesWaleDebug> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K17CompressibleNavierStokesWaleDebug();
	K17CompressibleNavierStokesWaleDebug(std::shared_ptr< Parameter> para, int level);

};
#endif 