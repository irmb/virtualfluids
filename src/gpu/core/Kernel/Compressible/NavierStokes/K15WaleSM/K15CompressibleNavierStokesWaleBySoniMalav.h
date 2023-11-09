#ifndef K15CompressibleNavierStokesWaleBySoniMalav_H
#define K15CompressibleNavierStokesWaleBySoniMalav_H

#include "Kernel/KernelImp.h"

class K15CompressibleNavierStokesWaleBySoniMalav : public KernelImp
{
public:
	static std::shared_ptr<K15CompressibleNavierStokesWaleBySoniMalav> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K15CompressibleNavierStokesWaleBySoniMalav();
	K15CompressibleNavierStokesWaleBySoniMalav(std::shared_ptr< Parameter> para, int level);

};
#endif 