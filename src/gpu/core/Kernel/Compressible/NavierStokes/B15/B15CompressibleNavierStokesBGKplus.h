#ifndef B15CompressibleNavierStokesBGKplus_H
#define B15CompressibleNavierStokesBGKplus_H

#include "Kernel/KernelImp.h"

class B15CompressibleNavierStokesBGKplus : public KernelImp
{
public:
	static std::shared_ptr<B15CompressibleNavierStokesBGKplus> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	B15CompressibleNavierStokesBGKplus();
	B15CompressibleNavierStokesBGKplus(std::shared_ptr< Parameter> para, int level);
};

#endif 