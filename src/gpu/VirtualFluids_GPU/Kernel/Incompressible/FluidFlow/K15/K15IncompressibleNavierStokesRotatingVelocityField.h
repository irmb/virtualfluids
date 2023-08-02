#ifndef K15IncompressibleNavierStokesRotatingVelocityField_H
#define K15IncompressibleNavierStokesRotatingVelocityField_H

#include "Kernel/KernelImp.h"

class K15IncompressibleNavierStokesRotatingVelocityField : public KernelImp
{
public:
	static std::shared_ptr<K15IncompressibleNavierStokesRotatingVelocityField> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();

private:
	K15IncompressibleNavierStokesRotatingVelocityField();
	K15IncompressibleNavierStokesRotatingVelocityField(std::shared_ptr< Parameter> para, int level);

};

#endif 