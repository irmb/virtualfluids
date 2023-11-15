#include "InitialConditionTaylorGreenVortexUz.h"

#include "Simulations/TaylorGreenVortexUz/TaylorGreenVortexUzParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>

InitialConditionTaylorGreenUz::InitialConditionTaylorGreenUz(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct)
{
    this->amp = simParaStruct->amplitude;
    this->l0 = simParaStruct->l0;
    this->lx = gridInfoStruct->lx;
    this->lz = gridInfoStruct->lz;
    this->rho = simParaStruct->rho0;
    this->uz = simParaStruct->uz;
}

std::shared_ptr<InitialConditionTaylorGreenUz> InitialConditionTaylorGreenUz::getNewInstance(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct)
{
    return std::shared_ptr<InitialConditionTaylorGreenUz>(new InitialConditionTaylorGreenUz(simParaStruct, gridInfoStruct));
}

real InitialConditionTaylorGreenUz::getInitVX(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vx = (amp * l0 * cos((real)2.0 * M_PI * z / lz) * sin((real)2.0 * M_PI * x / lx) / lx);
        return vx;
    }
    else
        return (real)0.0;

}

real InitialConditionTaylorGreenUz::getInitVY(int i, int level)
{
    return (real) 0.0;
}

real InitialConditionTaylorGreenUz::getInitVZ(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vz = (uz* l0 / lz) - (amp * l0 * lz * cos((real)2.0 * M_PI * x / lx) * sin((real)2.0 * M_PI * z / lz) / (lx*lx));
        return vz;
    }
    else
        return (real) 0.0;
}

real InitialConditionTaylorGreenUz::getInitROH(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real press = (amp*pow(l0, (real)2.0)*rho*(amp*pow(lx, (real)2.0)*pow(lz, (real)2.0)*pow(cos(((real)2.0 * M_PI*x) / lx), (real)2.0) - (real)4.0 * pow(lx, (real)2.0)*(pow(lx, (real)2.0) - pow(lz, (real)2.0))*uz*cos(((real)2.0 * M_PI*x) / lx)*sin(((real)2.0 * M_PI*z) / lz) - amp*pow(lz, (real)4.0)*pow(sin(((real)2.0 * M_PI*z) / lz), (real)2.0))) / ((real)2.0*pow(lx, (real)4.0)*pow(lz, (real)2.0));
        real rho = (real)3.0 * press;
        return rho;
    }
    else
        return (real) 0.0;
}

real InitialConditionTaylorGreenUz::getInitPRESS(int i, int level)
{
    //nicht ben�tigt, da Druck aus Dichte berechnet wird
    return (real) 0.0;
}