#include "InitialConditionImp.h"

#include "VirtualFluids_GPU/Parameter/Parameter.h"

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
	return (real)(para->getParH(level)->coordX_SP[i] - 1.0);
}

real InitialConditionImp::getYCoord(int i, int level)
{
	return (real)(para->getParH(level)->coordY_SP[i] - 1.0);
}

real InitialConditionImp::getZCoord(int i, int level)
{
	return (real)(para->getParH(level)->coordZ_SP[i] - 1.0);
}