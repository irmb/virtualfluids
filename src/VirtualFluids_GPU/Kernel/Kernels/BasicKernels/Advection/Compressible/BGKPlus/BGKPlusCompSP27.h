#ifndef BGK_PLUS_COMP_SP27_H
#define BGK_PLUS_COMP_SP27_H

#include "Kernel\KernelImp.h"

class BGKPlusCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<BGKPlusCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	BGKPlusCompSP27();
	BGKPlusCompSP27(std::shared_ptr< Parameter> para, int level);
};

#endif 