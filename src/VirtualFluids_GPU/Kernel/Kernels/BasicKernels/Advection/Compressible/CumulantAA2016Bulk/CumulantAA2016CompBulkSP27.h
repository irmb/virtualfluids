#ifndef CUMULANT_AA2016_COMP_BULK_SP27_H
#define CUMULANT_AA2016_COMP_BULK_SP27_H

#include "Kernel\KernelImp.h"

class CumulantAA2016CompBulkSP27 : public KernelImp
{
public:
	static std::shared_ptr<CumulantAA2016CompBulkSP27> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	CumulantAA2016CompBulkSP27();
	CumulantAA2016CompBulkSP27(std::shared_ptr< Parameter> para, int level);

};
#endif 