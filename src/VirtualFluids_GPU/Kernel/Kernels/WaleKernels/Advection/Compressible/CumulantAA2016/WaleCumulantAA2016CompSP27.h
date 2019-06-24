#ifndef WALE_CUMULANT_AA2016_COMP_SP27_H
#define WALE_CUMULANT_AA2016_COMP_SP27_H

#include "Kernel\KernelImp.h"

class WaleCumulantAA2016CompSP27 : public KernelImp
{
public:
	static std::shared_ptr<WaleCumulantAA2016CompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleCumulantAA2016CompSP27();
	WaleCumulantAA2016CompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 