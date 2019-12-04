#ifndef CUMULANT_F32018_COMP_SP27_H
#define CUMULANT_F32018_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantK20Comp : public KernelImp
{
public:
	static std::shared_ptr< CumulantK20Comp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK20Comp();
	CumulantK20Comp(std::shared_ptr< Parameter> para, int level);
};
#endif 