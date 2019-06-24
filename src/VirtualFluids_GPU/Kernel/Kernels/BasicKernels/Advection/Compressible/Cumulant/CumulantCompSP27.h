#ifndef CUMULANT_COMP_SP27_H
#define CUMULANT_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantCompSP27();
	CumulantCompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 