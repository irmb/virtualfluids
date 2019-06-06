#ifndef WALE_CUMULANT_AA2016_DEBUG_COMP_SP27_H
#define WALE_CUMULANT_AA2016_DEBUG_COMP_SP27_H

#include "Kernel\KernelImp.h"

class WaleCumulantAA2016DebugCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<WaleCumulantAA2016DebugCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleCumulantAA2016DebugCompSP27();
	WaleCumulantAA2016DebugCompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 