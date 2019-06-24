#ifndef WALE_CUMULANT_ONE_COMP_SP27_H
#define WALE_CUMULANT_ONE_COMP_SP27_H

#include "Kernel\KernelImp.h"

class WaleCumulantOneCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<WaleCumulantOneCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleCumulantOneCompSP27();
	WaleCumulantOneCompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 