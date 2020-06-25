#ifndef AD_INCOMP_7_H
#define AD_INCOMP_7_H

#include "Kernel/ADKernel.h"

class ADIncomp7 : public ADKernel
{
public:
	static std::shared_ptr<ADIncomp7> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	ADIncomp7();
	ADIncomp7(std::shared_ptr< Parameter> para, int level);
};
#endif 