#ifndef CUMULANT_AA2016_COMP_BULK_SP27_H
#define CUMULANT_AA2016_COMP_BULK_SP27_H

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