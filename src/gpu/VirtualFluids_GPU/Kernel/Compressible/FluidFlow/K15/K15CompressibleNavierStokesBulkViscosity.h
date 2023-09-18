#ifndef K15CompressibleNavierStokesBulkViscosity_H
#define K15CompressibleNavierStokesBulkViscosity_H

#include "Kernel/KernelImp.h"

class K15CompressibleNavierStokesBulkViscosity : public KernelImp
{
public:
	static std::shared_ptr<K15CompressibleNavierStokesBulkViscosity> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K15CompressibleNavierStokesBulkViscosity();
	K15CompressibleNavierStokesBulkViscosity(std::shared_ptr< Parameter> para, int level);
};

#endif 