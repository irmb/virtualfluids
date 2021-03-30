#ifndef CASCADE_INCOMP_P27_H
#define CASCADE_INCOMP_P27_H

#include "Kernel/KernelImp.h"

class CascadeIncompSP27 : public KernelImp
{
public:
	static std::shared_ptr<CascadeIncompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CascadeIncompSP27();
	CascadeIncompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 