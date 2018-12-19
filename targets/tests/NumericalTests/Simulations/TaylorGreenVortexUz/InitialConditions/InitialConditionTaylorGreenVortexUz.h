#ifndef INITIAL_CONDITION_TAYLORGREENVORTEX_UZ_H
#define INITIAL_CONDITION_TAYLORGREENVORTEX_UZ_H

#include "Utilities/InitialCondition/InitialConditionImp.h"

#include <memory>

class InitialConditionTaylorGreenUz :public InitialConditionImp
{
public:
	static std::shared_ptr< InitialConditionTaylorGreenUz> getNewInstance(real lx, real lz, real l0, real uz, real amplitude, real rho0);

	real getInitVX(int i, int level);
	real getInitVY(int i, int level);
	real getInitVZ(int i, int level);
	real getInitROH(int i, int level);
	real getInitPRESS(int i, int level);

private:
	InitialConditionTaylorGreenUz(real lx, real lz, real l0, real u0, real amplitude, real rho0);
	InitialConditionTaylorGreenUz() {};

	real Amp;
	real rho;
	real L0;
	real Lx, Lz;
	real uz;
};

#endif