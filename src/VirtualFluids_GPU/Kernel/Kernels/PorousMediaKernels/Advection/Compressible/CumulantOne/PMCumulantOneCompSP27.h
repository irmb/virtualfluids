#ifndef PM_CUMULANT_ONE_COMP_SP27_H
#define PM_CUMULANT_ONE_COMP_SP27_H

#include "Kernel\KernelImp.h"

class PorousMedia;

class PMCumulantOneCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<PMCumulantOneCompSP27> getNewInstance(std::shared_ptr< Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level);
	void run();

private:
	PMCumulantOneCompSP27();
	PMCumulantOneCompSP27(std::shared_ptr< Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level);

	std::vector<std::shared_ptr<PorousMedia>> pm;
};

#endif 