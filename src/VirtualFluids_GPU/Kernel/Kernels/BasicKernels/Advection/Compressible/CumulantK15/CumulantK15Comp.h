#ifndef CUMULANT_ONE_COMP_SP27_H
#define CUMULANT_ONE_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantK15Comp : public KernelImp
{
public:
	static std::shared_ptr< CumulantK15Comp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK15Comp();
	CumulantK15Comp(std::shared_ptr< Parameter> para, int level);
};
#endif 