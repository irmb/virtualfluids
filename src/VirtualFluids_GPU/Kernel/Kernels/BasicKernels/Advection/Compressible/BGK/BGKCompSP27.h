#ifndef BGK_COMP_SP27_H
#define BGK_COMP_SP27_H

#include "Kernel\KernelImp.h"

class BGKCompSP27 : public KernelImp
{
public:
	static std::shared_ptr<BGKCompSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();
	
private:
	BGKCompSP27();
	BGKCompSP27(std::shared_ptr< Parameter> para, int level);
};
#endif 