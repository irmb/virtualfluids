#include "InitialConditionTaylorGreenVortexUx.h"

#include "Simulations\TaylorGreenVortexUx\TaylorGreenVortexUxParameterStruct.h"
#include "Utilities\Structs\GridInformationStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

InitialConditionTaylorGreenUx::InitialConditionTaylorGreenUx(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct>  gridInfoStruct)
{
	this->Amp = simParaStruct->amplitude;
	this->L0 = simParaStruct->l0;
	this->Lx = gridInfoStruct->lx;
	this->Lz = gridInfoStruct->lz;
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
		real vx = (ux* L0 / Lx + (Amp * L0 * cos((real)2.0 * M_PI * z / Lz) * sin((real)2.0 * M_PI * x / Lx) / Lx));
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
		real vz = (-Amp * L0 * Lz * cos((real)2.0 * M_PI * x / Lx) * sin((real)2.0 * M_PI * z / Lz) / (Lx*Lx));
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
		real press = (Amp*Amp * L0*L0 * rho * ((Lx*Lx * cos((real)4.0 * M_PI * x / Lx)) + (Lz*Lz * cos((real)4.0 * M_PI * z / Lz))) / ((real)4.0 * Lx*Lx*Lx*Lx));
		return press;
	}
	else
		return (real) 0.0;
}

real InitialConditionTaylorGreenUx::getInitPRESS(int i, int level)
{
	//nicht benötigt, da Druck aus Dichte berechnet wird
	return (real) 0.0;
}