#ifndef CUMULANT_F3_COMP_SP27_H
#define CUMULANT_F3_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantK18Comp : public KernelImp
{
public:
	static std::shared_ptr< CumulantK18Comp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK18Comp();
	CumulantK18Comp(std::shared_ptr< Parameter> para, int level);
};

#endif 