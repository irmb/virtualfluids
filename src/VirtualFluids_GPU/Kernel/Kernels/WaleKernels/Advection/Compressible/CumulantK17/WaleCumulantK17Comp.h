#ifndef WALE_CUMULANT_AA2016_COMP_SP27_H
#define WALE_CUMULANT_AA2016_COMP_SP27_H

#include "Kernel\KernelImp.h"

class WaleCumulantK17Comp : public KernelImp
{
public:
	static std::shared_ptr<WaleCumulantK17Comp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	WaleCumulantK17Comp();
	WaleCumulantK17Comp(std::shared_ptr< Parameter> para, int level);

};
#endif 