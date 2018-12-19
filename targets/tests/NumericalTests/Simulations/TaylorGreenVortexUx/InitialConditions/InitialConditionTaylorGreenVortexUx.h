#ifndef INITIAL_CONDITION_TAYLORGREENVORTEX_UX_H
#define INITIAL_CONDITION_TAYLORGREENVORTEX_UX_H

#include "Utilities/InitialCondition/InitialConditionImp.h"

#include <memory>

class InitialConditionTaylorGreenUx :public InitialConditionImp
{
public:
	static std::shared_ptr< InitialConditionTaylorGreenUx> getNewInstance(real lx, real lz, real l0, real ux, real amplitude, real rho0);

	real getInitVX(int i, int level);
	real getInitVY(int i, int level);
	real getInitVZ(int i, int level);
	real getInitROH(int i, int level);
	real getInitPRESS(int i, int level);

private:
	InitialConditionTaylorGreenUx(real lx, real lz, real l0, real ux, real amplitude, real rho0);
	InitialConditionTaylorGreenUx() {};

	real Amp;
	real rho;
	real L0;
	real Lx, Lz;
	real ux;
};

#endif