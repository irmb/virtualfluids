#ifndef K15CompressibleNavierStokesSponge_H
#define K15CompressibleNavierStokesSponge_H

#include "Kernel/KernelImp.h"

class K15CompressibleNavierStokesSponge : public KernelImp
{
public:
	static std::shared_ptr<K15CompressibleNavierStokesSponge> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K15CompressibleNavierStokesSponge();
	K15CompressibleNavierStokesSponge(std::shared_ptr< Parameter> para, int level);

};
#endif 