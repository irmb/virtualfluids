#ifndef CUMULANT_AA2016_COMP_SP27_H
#define CUMULANT_AA2016_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantK17Comp : public KernelImp
{
public:
	static std::shared_ptr<CumulantK17Comp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK17Comp();
	CumulantK17Comp(std::shared_ptr< Parameter> para, int level);
};

#endif 
