#ifndef CUMULANT_ONE_COMP_BULK_SP27_H
#define CUMULANT_ONE_COMP_BULK_SP27_H

#include "Kernel\KernelImp.h"

class CumulantOneCompBulkSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantOneCompBulkSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantOneCompBulkSP27();
	CumulantOneCompBulkSP27(std::shared_ptr< Parameter> para, int level);
};

#endif 