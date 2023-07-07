#ifndef K18CompressibleNavierStokes_H
#define K18CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

class K18CompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr< K18CompressibleNavierStokes> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K18CompressibleNavierStokes();
	K18CompressibleNavierStokes(std::shared_ptr< Parameter> para, int level);
};

#endif 