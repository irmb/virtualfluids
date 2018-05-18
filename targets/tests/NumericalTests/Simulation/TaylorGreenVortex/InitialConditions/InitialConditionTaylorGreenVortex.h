#ifndef INITIAL_CONDITION_TAYLOR_GREEN_VORTEX_H
#define INITIAL_CONDITION_TAYLOR_GREEN_VORTEX_H

#include "Utilities/InitialCondition/InitialConditionImp.h"

#include <memory>

class InitialConditionTaylorGreen :public InitialConditionImp
{
public:
	static std::shared_ptr< InitialConditionTaylorGreen> getNewInstance(real lx, real lz, real l0, real u0, real amplitude, real rho0);

	real getInitVX(int i, int level);
	real getInitVY(int i, int level);
	real getInitVZ(int i, int level);
	real getInitROH(int i, int level);
	real getInitPRESS(int i, int level);

private:
	InitialConditionTaylorGreen(real lx, real lz, real l0, real u0, real amplitude, real rho0);
	InitialConditionTaylorGreen() {};

	real Amp;
	real rho;
	real L0;
	real Lx, Lz;
	real u0;
};

#endif