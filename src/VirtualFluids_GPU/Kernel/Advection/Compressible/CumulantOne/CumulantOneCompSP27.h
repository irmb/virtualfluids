#ifndef CUMULANT_ONE_COMP_SP27_H
#define CUMULANT_ONE_COMP_SP27_H

#include "Kernel/Kernel.h"

#include <memory>

class CumulantOneCompSP27 : public Kernel
{
public:
	static std::shared_ptr< Kernel> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();
	bool checkParameter();

private:
	CumulantOneCompSP27();
	CumulantOneCompSP27(std::shared_ptr< Parameter> para, int level);
	std::shared_ptr< Parameter> para;
	int level;
};
#endif 