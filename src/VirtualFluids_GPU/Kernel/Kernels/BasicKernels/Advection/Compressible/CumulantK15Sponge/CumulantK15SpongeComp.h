#ifndef CUMULANT_K15_SPONGE_COMP_H
#define CUMULANT_K15_SPONGE_COMP_H

#include "Kernel\KernelImp.h"

class CumulantK15SpongeComp : public KernelImp
{
public:
	static std::shared_ptr<CumulantK15SpongeComp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK15SpongeComp();
	CumulantK15SpongeComp(std::shared_ptr< Parameter> para, int level);

};
#endif 