#ifndef AD_INCOMP_27_H
#define AD_INCOMP_27_H

#include "Kernel\ADKernel.h"


class ADIncomp27 : public ADKernel
{
public:
	static std::shared_ptr<ADIncomp27> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run();

private:
	ADIncomp27();
	ADIncomp27(std::shared_ptr< Parameter> para, int level);
};
#endif 