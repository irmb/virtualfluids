#include "InitialConditionTaylorGreenVortexUx.h"

#include "Simulations/TaylorGreenVortexUx/TaylorGreenVortexUxParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

InitialConditionTaylorGreenUx::InitialConditionTaylorGreenUx(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct>  gridInfoStruct)
{
	this->amp = simParaStruct->amplitude;
	this->l0 = simParaStruct->l0;
	this->lx = gridInfoStruct->lx;
	this->lz = gridInfoStruct->lz;
	this->rho = simParaStruct->rho0;
	this->ux = simParaStruct->ux;
}

std::shared_ptr<InitialConditionTaylorGreenUx> InitialConditionTaylorGreenUx::getNewInstance(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct>  gridInfoStruct)
{
	return std::shared_ptr<InitialConditionTaylorGreenUx>(new InitialConditionTaylorGreenUx(simParaStruct, gridInfoStruct));
}

real InitialConditionTaylorGreenUx::getInitVX(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
	{
		real vx = (ux* l0 / lx + (amp * l0 * cos((real)2.0 * M_PI * z / lz) * sin((real)2.0 * M_PI * x / lx) / lx));
		return vx;
	}
	else
		return (real)0.0;

}

real InitialConditionTaylorGreenUx::getInitVY(int i, int level)
{
	return (real) 0.0;
}

real InitialConditionTaylorGreenUx::getInitVZ(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
	{
		real vz = (-amp * l0 * lz * cos((real)2.0 * M_PI * x / lx) * sin((real)2.0 * M_PI * z / lz) / (lx*lx));
		return vz;
	}
	else
		return (real) 0.0;
}

real InitialConditionTaylorGreenUx::getInitROH(int i, int level)
{
	real x = getXCoord(i, level);
	real y = getYCoord(i, level);
	real z = getZCoord(i, level);
	if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
	{
		real press = (amp*pow(l0, (real)2.0)*rho*(amp*pow(lz, (real)2.0)*pow(cos(((real)2.0 * M_PI*z) / lz), (real)2.0) - (real)4.0 * (pow(lx, (real)2.0) - pow(lz, (real)2.0))*ux*cos(((real)2.0 * M_PI*z) / lz)*sin(((real)2.0 * M_PI*x) / lx) - amp*pow(lx, (real)2.0)*pow(sin(((real)2.0 * M_PI*x) / lx), (real)2.0))) / ((real)2.0*pow(lx, (real)4.0));
		real rho = (real)3.0 * press;
		return rho;
	}
	else
		return (real) 0.0;
}

real InitialConditionTaylorGreenUx::getInitPRESS(int i, int level)
{
	//nicht benötigt, da Druck aus Dichte berechnet wird
	return (real) 0.0;
}