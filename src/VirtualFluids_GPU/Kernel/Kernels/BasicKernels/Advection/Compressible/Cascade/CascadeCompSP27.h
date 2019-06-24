#ifndef CASCADE_COMP_SP27_H
#define CASCADE_COMP_SP27_H

#include "Kernel\KernelImp.h"

class CascadeCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<CascadeCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CascadeCompSP27();
	CascadeCompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 