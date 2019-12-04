#ifndef CUMULANT_K17_BULK_COMP_H
#define CUMULANT_K17_BULK_COMP_H

#include "Kernel\KernelImp.h"

class CumulantK17BulkComp : public KernelImp
{
public:
	static std::shared_ptr<CumulantK17BulkComp> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantK17BulkComp();
	CumulantK17BulkComp(std::shared_ptr< Parameter> para, int level);

};
#endif 