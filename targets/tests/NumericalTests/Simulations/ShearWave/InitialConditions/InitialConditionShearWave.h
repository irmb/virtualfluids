#ifndef INITIAL_CONDITION_SHEAR_WAVE_H
#define	INITIAL_CONDITION_SHEAR_WAVE_H

#include "Utilities/InitialCondition/InitialConditionImp.h"

class InitialConditionShearWave :public InitialConditionImp
{
public:
	InitialConditionShearWave(real lx, real lz, real l0, real u0, real v0, real rho0);
	real getInitVX(int i, int level);
	real getInitVY(int i, int level);
	real getInitVZ(int i, int level);
	real getInitROH(int i, int level);
	real getInitPRESS(int i, int level);

private:
	InitialConditionShearWave();
	real rho;
	real l0;
	real lx, lz;
	real u0, v0;
};
#endif 