#ifndef WALE_CUMULANT_K15_COMP_H
#define WALE_CUMULANT_K15_COMP_H

#include "Kernel/KernelImp.h"

class WaleCumulantK15Comp : public KernelImp
{
public:
	static std::shared_ptr<WaleCumulantK15Comp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleCumulantK15Comp();
	WaleCumulantK15Comp(std::shared_ptr< Parameter> para, int level);

};
#endif 