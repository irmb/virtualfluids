#ifndef BGK_INCOMP_SP27_H
#define BGK_INCOMP_SP27_H

#include "Kernel\KernelImp.h"


class BGKIncompSP27 : public KernelImp
{
public:
	static std::shared_ptr<BGKIncompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	BGKIncompSP27();
	BGKIncompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 