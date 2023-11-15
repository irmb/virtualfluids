#include "InitialConditionImp.h"

#include "gpu/core/Parameter/Parameter.h"

void InitialConditionImp::setParameter(std::shared_ptr<Parameter> para)
{
    this->para = para;
}

void InitialConditionImp::init(const int level)
{
    XCoordStopNode = para->getGridX().at(level) - 1.0 + 0.5;
    YCoordStopNode = para->getGridY().at(level) - 1.0 + 0.5;
    ZCoordStopNode = para->getGridZ().at(level) - 1.0 + 0.5;
}

real InitialConditionImp::getXCoord(int i, int level)
{
    return para->getParH(level)->coordinateX[i] - 1.0;
}

real InitialConditionImp::getYCoord(int i, int level)
{
    return para->getParH(level)->coordinateY[i] - 1.0;
}

real InitialConditionImp::getZCoord(int i, int level)
{
    return para->getParH(level)->coordinateZ[i] - 1.0;
}