#ifndef CUMULANT_F32018_COMP_SP27_H
#define CUMULANT_F32018_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CumulantF32018CompSP27 : public KernelImp
{
public:
	static std::shared_ptr< CumulantF32018CompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantF32018CompSP27();
	CumulantF32018CompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 