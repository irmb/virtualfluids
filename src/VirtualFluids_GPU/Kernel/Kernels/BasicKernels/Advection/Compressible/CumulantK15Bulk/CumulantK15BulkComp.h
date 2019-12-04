#ifndef CUMULANT_ONE_COMP_BULK_SP27_H
#define CUMULANT_ONE_COMP_BULK_SP27_H

#include "Kernel\KernelImp.h"

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