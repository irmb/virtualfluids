#ifndef CUMULANT_ALL4_COMP_SP27_H
#define CUMULANT_ALL4_COMP_SP27_H

#include "Kernel/Kernel.h"

#include <memory>

class CumulantAll4CompSP27 : public Kernel
{
public:
	static std::shared_ptr< Kernel> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();
	bool checkParameter();

private:
	CumulantAll4CompSP27(std::shared_ptr< Parameter> para, int level);
	std::shared_ptr< Parameter> para;
	int level;

};

#endif 