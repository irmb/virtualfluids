#ifndef CUMULANT_ONE_COMP_SP27_H
#define CUMULANT_ONE_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantOneCompSP27 : public KernelImp
{
public:
	static std::shared_ptr< CumulantOneCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantOneCompSP27();
	CumulantOneCompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 