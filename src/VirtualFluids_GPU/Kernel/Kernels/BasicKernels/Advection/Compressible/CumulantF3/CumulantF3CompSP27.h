#ifndef CUMULANT_F3_COMP_SP27_H
#define CUMULANT_F3_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantF3CompSP27 : public KernelImp
{
public:
	static std::shared_ptr< CumulantF3CompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantF3CompSP27();
	CumulantF3CompSP27(std::shared_ptr< Parameter> para, int level);
};

#endif 