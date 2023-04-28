#ifndef CUMULANT_K17_BULK_COMP_H
#define CUMULANT_K17_BULK_COMP_H

#include "Kernel/KernelImp.h"

class K17CompressibleNavierStokesBulkViscosity : public KernelImp
{
public:
	static std::shared_ptr<K17CompressibleNavierStokesBulkViscosity> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K17CompressibleNavierStokesBulkViscosity();
	K17CompressibleNavierStokesBulkViscosity(std::shared_ptr< Parameter> para, int level);

};
#endif 