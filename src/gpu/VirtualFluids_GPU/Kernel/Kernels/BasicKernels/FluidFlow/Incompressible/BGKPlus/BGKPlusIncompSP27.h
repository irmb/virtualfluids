#ifndef BGK_PLUS_INCOMP_SP27_H
#define BGK_PLUS_INCOMP_SP27_H

#include "Kernel/KernelImp.h"

class BGKPlusIncompSP27 : public KernelImp
{
public:
	static std::shared_ptr<BGKPlusIncompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	BGKPlusIncompSP27();
	BGKPlusIncompSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 