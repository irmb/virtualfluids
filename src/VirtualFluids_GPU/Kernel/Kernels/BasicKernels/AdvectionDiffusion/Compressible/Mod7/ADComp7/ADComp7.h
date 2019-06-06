#ifndef AD_COMP_7_H
#define AD_COMP_7_H

#include "Kernel\ADKernel.h"

class ADComp7 : public ADKernel
{
public:
	static std::shared_ptr<ADComp7> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	ADComp7();
	ADComp7(std::shared_ptr< Parameter> para, int level);
};
#endif 