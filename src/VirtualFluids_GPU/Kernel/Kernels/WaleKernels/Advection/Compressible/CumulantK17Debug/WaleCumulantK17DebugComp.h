#ifndef WALE_CUMULANT_AA2016_DEBUG_COMP_SP27_H
#define WALE_CUMULANT_AA2016_DEBUG_COMP_SP27_H

#include "Kernel\KernelImp.h"

class WaleCumulantK17DebugComp : public KernelImp
{
public:
	static std::shared_ptr<WaleCumulantK17DebugComp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleCumulantK17DebugComp();
	WaleCumulantK17DebugComp(std::shared_ptr< Parameter> para, int level);

};
#endif 