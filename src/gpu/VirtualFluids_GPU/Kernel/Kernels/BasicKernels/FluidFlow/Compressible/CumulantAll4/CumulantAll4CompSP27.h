#ifndef CUMULANT_ALL4_COMP_SP27_H
#define CUMULANT_ALL4_COMP_SP27_H

#include "Kernel/KernelImp.h"

class CumulantAll4CompSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantAll4CompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();
	
private:
	CumulantAll4CompSP27();
	CumulantAll4CompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 