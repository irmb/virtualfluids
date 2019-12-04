#ifndef CUMULANT_K18_COMP_H
#define CUMULANT_K18_COMP_H

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