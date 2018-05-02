#include "InitialConditionShearWave.h"

#define _USE_MATH_DEFINES
#include <math.h>


InitialConditionShearWave::InitialConditionShearWave(real lx, real lz, real l0, real u0, real v0, real rho0)
{
	this->l0 = l0;
	this->lx = lx;
	this->lz = lz;
	this->rho = rho0;
	this->u0 = u0;
	this->v0 = v0;
}

real InitialConditionShearWave::getInitVX(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordstopnode) && (y != YCoordstopnode) && (z != ZCoordstopnode))
	{
		real vx = l0 * u0 / lx;
		return vx;
	}
	else
		return (real)0.0;

}

real InitialConditionShearWave::getInitVY(int i, int level)
{
	return (real) 0.0;
}

real InitialConditionShearWave::getInitVZ(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordstopnode) && (y != YCoordstopnode) && (z != ZCoordstopnode))
	{
		real vz = v0 * l0 / lx * cos((real)2.0 * M_PI * z /lz) * sin((real)2.0 * M_PI * x / lx);
		return vz;
	}
	else
		return (real) 0.0;
}

real InitialConditionShearWave::getInitROH(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordstopnode) && (y != YCoordstopnode) && (z != ZCoordstopnode))
	{
		real press = (l0*l0 * v0 * rho * sin(((real)2.0 * M_PI * z) / lz) * ((real)-4.0 * lz * u0 * cos((2 * M_PI * x) / lx) + lx * v0 * sin(((real)2.0 * M_PI * x) / lx)*sin(((real)2.0 * M_PI * x) / lx) * sin(((real)2.0 * M_PI * z) / lz))) / ((real)2.0 * lx*lx*lx);
		return press;
	}
	else
		return (real) 0.0;
}

real InitialConditionShearWave::getInitPRESS(int i, int level)
{
	//nicht benötigt, da Druck aus Dichte berechnet wird
	return (real) 0.0;
}