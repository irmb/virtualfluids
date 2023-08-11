#ifndef K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants_H
#define K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants_H

#include "Kernel/KernelImp.h"

class K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants : public KernelImp
{
public:
	static std::shared_ptr<K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();
	
private:
	K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants();
	K17CompressibleNavierStokesSecondDerivatesFrom5thCumulants(std::shared_ptr< Parameter> para, int level);
};
#endif 