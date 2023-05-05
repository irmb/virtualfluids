#ifndef K17CompressibleNavierStokesBulkViscosity_H
#define K17CompressibleNavierStokesBulkViscosity_H

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