#ifndef CUMULANT_ONE_INCOMP_SP27_H
#define CUMULANT_ONE_INCOMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantOneIncompSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantOneIncompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantOneIncompSP27();
	CumulantOneIncompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 