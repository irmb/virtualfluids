#ifndef CUMULANT_AA2016_COMP_SP27_H
#define CUMULANT_AA2016_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantAA2016CompSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantAA2016CompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantAA2016CompSP27();
	CumulantAA2016CompSP27(std::shared_ptr< Parameter> para, int level);
};

#endif 
