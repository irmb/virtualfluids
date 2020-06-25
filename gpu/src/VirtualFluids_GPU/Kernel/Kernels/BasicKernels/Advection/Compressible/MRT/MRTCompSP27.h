#ifndef MRT_COMP_SP27_H
#define MRT_COMP_SP27_H

#include "Kernel/KernelImp.h"


class MRTCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<MRTCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	MRTCompSP27();
	MRTCompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 