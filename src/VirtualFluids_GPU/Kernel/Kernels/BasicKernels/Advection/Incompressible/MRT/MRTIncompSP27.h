#ifndef MRT_INCOMP_SP27_H
#define MRT_INCOMP_SP27_H

#include "Kernel\KernelImp.h"

class MRTIncompSP27 : public KernelImp
{
public:
	static std::shared_ptr<MRTIncompSP27> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run();

private:
	MRTIncompSP27();
	MRTIncompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 