#ifndef CUMULANT_ONE_COMP_SPONGE_SP27_H
#define CUMULANT_ONE_COMP_SPONGE_SP27_H

#include "Kernel\KernelImp.h"

class CumulantOneCompSpongeSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantOneCompSpongeSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantOneCompSpongeSP27();
	CumulantOneCompSpongeSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 