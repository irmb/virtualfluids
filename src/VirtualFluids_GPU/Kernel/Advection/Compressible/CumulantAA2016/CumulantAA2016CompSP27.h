#ifndef CUMULANT_AA2016_COMP_SP27_H
#define CUMULANT_AA2016_COMP_SP27_H

#include "Kernel/Kernel.h"

#include <memory>

class CumulantAA2016CompSP27 : public Kernel
{
public:
	static std::shared_ptr<Kernel> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();
	bool checkParameter();

	

private:
	CumulantAA2016CompSP27();
	CumulantAA2016CompSP27(std::shared_ptr< Parameter> para, int level);
	std::shared_ptr< Parameter> para;
	int level;
};

#endif 
