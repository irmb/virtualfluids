#ifndef AD_COMP_27_H
#define AD_COMP_27_H

#include "Kernel/ADKernel.h"

class ADComp27 : public ADKernel
{
public:
	static std::shared_ptr<ADComp27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	ADComp27();
	ADComp27(std::shared_ptr< Parameter> para, int level);
};
#endif 