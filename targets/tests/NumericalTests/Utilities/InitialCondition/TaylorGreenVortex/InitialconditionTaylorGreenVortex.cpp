#include "InitialconditionTaylorGreenVortex.h"

#define _USE_MATH_DEFINES
#include <math.h>

InitialConditionTaylorGreen::InitialConditionTaylorGreen(real lx, real lz, real l0, real u0, real amplitude, real rho0)
{
	this->Amp = amplitude;
	this->L0 = l0;
	this->Lx = lx;
	this->Lz = lz;
	this->rho = rho0;
	this->u0 = u0;
}

real InitialConditionTaylorGreen::getInitVX(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordstopnode) && (y != YCoordstopnode) && (z != ZCoordstopnode))
	{
		real vx = (u0* L0 / Lx + (Amp * L0 * cos((real)2.0 * M_PI * z / Lz) * sin((real)2.0 * M_PI * x / Lx) / Lx));
		return vx;
	}
	else
		return (real)0.0;

}

real InitialConditionTaylorGreen::getInitVY(int i, int level)
{
	return (real) 0.0;
}

real InitialConditionTaylorGreen::getInitVZ(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordstopnode) && (y != YCoordstopnode) && (z != ZCoordstopnode))
	{
		real vz = (-Amp * L0 * Lz * cos((real)2.0 * M_PI * x / Lx) * sin((real)2.0 * M_PI * z / Lz) / (Lx*Lx));
		return vz;
	}
	else
		return (real) 0.0;
}

real InitialConditionTaylorGreen::getInitROH(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordstopnode) && (y != YCoordstopnode) && (z != ZCoordstopnode))
	{
		real press = (Amp*Amp * L0*L0 * rho * ((Lx*Lx * cos((real)4.0 * M_PI * x / Lx)) + (Lz*Lz * cos((real)4.0 * M_PI * z / Lz))) / ((real)4.0 * Lx*Lx*Lx*Lx));
		return press;
	}
	else
		return (real) 0.0;
}

real InitialConditionTaylorGreen::getInitPRESS(int i, int level)
{
	//nicht benötigt, da Druck aus Dichte berechnet wird
	return (real) 0.0;
}