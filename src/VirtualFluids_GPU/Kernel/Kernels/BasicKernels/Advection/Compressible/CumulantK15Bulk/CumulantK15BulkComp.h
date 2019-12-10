#ifndef CUMULANT_K15_BULK_COMP_H
#define CUMULANT_K15_BULK_COMP_H

#include "Kernel/KernelImp.h"

class CumulantK15BulkComp : public KernelImp
{
public:
	static std::shared_ptr<CumulantK15BulkComp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK15BulkComp();
	CumulantK15BulkComp(std::shared_ptr< Parameter> para, int level);
};

#endif 