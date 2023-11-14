#include "InitialConditionShearWave.h"

#include "Simulations/ShearWave/ShearWaveParameterStruct.h"
#include "Utilities/Structs/GridInformationStruct.h"

#define _USE_MATH_DEFINES
#include <math.h>


InitialConditionShearWave::InitialConditionShearWave(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct)
{
    this->l0 = simParaStruct->l0;
    this->lx = gridInfoStruct->lx;
    this->lz = gridInfoStruct->lz;
    this->rho = simParaStruct->rho0;
    this->u0 = simParaStruct->ux;
    this->v0 = simParaStruct->uz;
}

std::shared_ptr<InitialConditionShearWave> InitialConditionShearWave::getNewInstance(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct)
{
    return std::shared_ptr<InitialConditionShearWave>(new InitialConditionShearWave(simParaStruct, gridInfoStruct));
}

real InitialConditionShearWave::getInitVX(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vx = l0 * u0 / lx;
        return vx;
    }
    else
        return (real)0.0;

}

real InitialConditionShearWave::getInitVY(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real vy = v0 * l0 / lx * cos((real)2.0 * M_PI * z / lz) * sin((real)2.0 * M_PI * x / lx);
        return vy;
    }
    else
        return (real) 0.0;
}

real InitialConditionShearWave::getInitVZ(int i, int level)
{
    return (real) 0.0;
}

real InitialConditionShearWave::getInitROH(int i, int level)
{
    real x = getXCoord(i, level);
    real y = getYCoord(i, level);
    real z = getZCoord(i, level);
    if ((i != 0) && (x != XCoordStopNode) && (y != YCoordStopNode) && (z != ZCoordStopNode))
    {
        real press = (l0*l0 * v0 * rho * sin(((real)2.0 * M_PI * z) / lz) * ((real)-4.0 * lz * u0 * cos(((real)2.0 * M_PI * x) / lx) + lx * v0 * sin(((real)2.0 * M_PI * x) / lx)*sin(((real)2.0 * M_PI * x) / lx) * sin(((real)2.0 * M_PI * z) / lz))) / ((real)2.0 * lx*lx*lx);
        real rho = (real)3.0 * press;
        return 0.0;
    }
    else
        return (real) 0.0;
}

real InitialConditionShearWave::getInitPRESS(int i, int level)
{
    //nicht ben�tigt, da Druck aus Dichte berechnet wird
    return (real) 0.0;
}